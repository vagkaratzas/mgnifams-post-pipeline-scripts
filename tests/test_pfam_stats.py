import csv
import gzip
import importlib.util
import io
import json
from collections import Counter
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "pfam_stats.py"


def load_module():
    """Import the standalone script without requiring `bin` to be a package."""
    spec = importlib.util.spec_from_file_location("pfam_stats", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_module()


# ---------------------------------------------------------------------------
# fixtures / helpers
# ---------------------------------------------------------------------------
def mkrow(mgyp, seq, full="false", csize=1, hits=()):
    """One CSV line, quoted exactly as the real sequence_explorer_protein.csv is."""
    meta = json.dumps({
        "s": [[112156, [[781358, [[1218588979, 1, 396, 1]]]]]],
        "b": [[120, 1]],
        "p": [list(h) for h in hits],
    })
    buf = io.StringIO()
    csv.writer(buf, lineterminator="").writerow([mgyp, seq, full, csize, meta])
    return buf.getvalue()


def hit(acc, hmm_from=1, hmm_to=60, s_from=1, s_to=60, evalue=1e-8, score=39.0):
    """A "p" entry: [acc, evalue, score, hmm_from, hmm_to, s_from, s_to]."""
    return [acc, evalue, score, hmm_from, hmm_to, s_from, s_to]


def write_csv(path, rows, header=True):
    lines = ["mgyp,sequence,full_length,cluster_size,metadata"] if header else []
    lines.extend(rows)
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(data)
    return data


def scan(path, shards=None):
    """Run the map/reduce over explicit (start, end) shards; default one shard."""
    size = path.stat().st_size
    if shards is None:
        shards = [(0, size)]
    tasks = [(str(path), s, e) for s, e in shards]
    return M.combine([M.worker(t) for t in tasks])


# ---------------------------------------------------------------------------
# byte-offset sharding -- the regression that motivated these tests
# ---------------------------------------------------------------------------
def test_every_split_point_preserves_every_row(tmp_path):
    """Splitting into two shards at ANY byte offset must not gain or lose a row.

    The original code seeked to `start` and dropped one line; when `start` was
    itself a line start that dropped a whole record the previous shard had
    already stopped short of.
    """
    p = tmp_path / "fw.csv"
    rows = [mkrow(f"{i:010d}", "MKVLA", hits=[hit("PF00462")]) for i in range(10)]
    data = write_csv(p, rows)

    baseline = scan(p)
    assert baseline["n_seqs"] == 10

    for split in range(1, len(data)):
        tot = scan(p, [(0, split), (split, len(data))])
        assert tot["n_seqs"] == 10, f"row lost/duplicated splitting at byte {split}"
        assert tot["n_residues"] == baseline["n_residues"]
        assert tot["n_with_pfam"] == baseline["n_with_pfam"]


@pytest.mark.parametrize("nshards", [1, 2, 3, 4, 5, 7, 10])
def test_uniform_shards_match_single_pass(tmp_path, nshards):
    """Fixed-width rows make shard boundaries land exactly on line starts."""
    p = tmp_path / "fw.csv"
    rows = [mkrow(f"{i:010d}", "MKVLA", hits=[hit("PF00462")]) for i in range(20)]
    data = write_csv(p, rows)

    step = len(data) // nshards
    shards = [(i * step, len(data) if i == nshards - 1 else (i + 1) * step)
              for i in range(nshards)]
    assert scan(p, shards)["n_seqs"] == 20


def test_header_skipped_in_whichever_shard_holds_it(tmp_path):
    """Concatenated shard files repeat the header; every copy must be dropped."""
    p = tmp_path / "multi.csv"
    rows = []
    for block in range(3):
        rows.append("mgyp,sequence,full_length,cluster_size,metadata")
        rows.extend(mkrow(f"{block}{i:09d}", "MKVLA") for i in range(5))
    data = write_csv(p, rows)

    assert scan(p)["n_seqs"] == 15
    mid = len(data) // 2
    assert scan(p, [(0, mid), (mid, len(data))])["n_seqs"] == 15


def test_build_tasks_shards_cover_the_file_exactly(tmp_path):
    p = tmp_path / "big.csv"
    p.write_bytes(b"x" * (200 << 20))          # 200 MB, sparse-ish but real enough
    tasks = M.build_tasks([str(p)], jobs=4)

    assert len(tasks) > 1, "a 200 MB file should shard"
    bounds = [(s, e) for _, s, e in tasks]
    assert bounds[0][0] == 0
    assert bounds[-1][1] == p.stat().st_size
    for (_, prev_end), (next_start, _) in zip(bounds, bounds[1:]):
        assert prev_end == next_start, "shards must be contiguous and non-overlapping"


def test_gzip_input_matches_plain_input(tmp_path):
    rows = [mkrow(f"{i:010d}", "MKVLAGG", hits=[hit("PF00462")]) for i in range(12)]
    plain = tmp_path / "a.csv"
    write_csv(plain, rows)
    gz = tmp_path / "a.csv.gz"
    gz.write_bytes(gzip.compress(plain.read_bytes()))

    a = M.worker((str(plain), 0, plain.stat().st_size))
    b = M.worker((str(gz), None, None))
    for k in ("n_seqs", "n_residues", "n_with_pfam", "n_residues_pfam"):
        assert a[k] == b[k]
    assert dict(a["fam_seqs"]) == dict(b["fam_seqs"])


# ---------------------------------------------------------------------------
# interval union -- drives every residue-coverage number in the report
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("intervals,expected", [
    ([], 0),
    ([(1, 10)], 10),
    ([(1, 5), (6, 10)], 10),          # adjacent, must coalesce (inclusive coords)
    ([(1, 5), (7, 10)], 9),           # one-residue gap must NOT coalesce
    ([(1, 10), (3, 5)], 10),          # fully contained
    ([(1, 10), (5, 20)], 20),         # overlapping
    ([(5, 20), (1, 10)], 20),         # unsorted input
    ([(1, 1)], 1),                    # single residue
    ([(10, 20), (1, 5), (30, 40)], 5 + 11 + 11),
])
def test_merge_len(intervals, expected):
    assert M.merge_len(list(intervals)) == expected


# ---------------------------------------------------------------------------
# metadata parsing -- "p" is bracket-matched, "s" is deliberately never decoded
# ---------------------------------------------------------------------------
def test_extract_p_skips_the_s_block():
    meta = json.dumps({
        "s": [[112156, [[781358, [[1218588979, 1, 396, 1]]]]]],
        "b": [[120, 1]],
        "p": [["PF00462", 7.5e-8, 39.0, 2, 58, 53, 113],
              ["PF03960", 1.9e-4, 28.2, 3, 46, 5, 82]],
    })
    hits = M.extract_p(meta)
    assert [h[0] for h in hits] == ["PF00462", "PF03960"]
    assert hits[0][5:7] == [53, 113]


def test_extract_p_absent_or_empty():
    assert M.extract_p(json.dumps({"s": [[1, 2]], "b": [[120, 1]]})) == []
    assert M.extract_p(json.dumps({"s": [[1, 2]], "p": []})) == []
    assert M.extract_p("") == []


def test_extract_p_when_p_is_not_the_last_key():
    meta = '{"p":[["PF00462",1e-8,39.0,1,60,1,60]],"s":[[1,[[2,[[3,1,9,1]]]]]]}'
    assert M.extract_p(meta)[0][0] == "PF00462"


# ---------------------------------------------------------------------------
# accumulation
# ---------------------------------------------------------------------------
def test_consume_basic_counters():
    rows = [
        mkrow("1", "MKVLA", full="true", csize=6, hits=[hit("PF00462", s_from=1, s_to=3)]),
        mkrow("2", "MKV", full="false", csize=2),
        mkrow("3", "MKVLAGGWW", full="TRUE", csize=1, hits=[hit("PF00462", s_from=1, s_to=9)]),
    ]
    acc = M.consume(rows, M.blank())

    assert acc["n_seqs"] == 3
    assert acc["n_residues"] == 5 + 3 + 9
    assert acc["n_full_length"] == 2          # case-insensitive
    assert acc["n_cluster_members"] == 9
    assert acc["n_with_pfam"] == 2            # the unannotated row is excluded
    assert acc["n_residues_pfam"] == 3 + 9
    assert acc["fam_seqs"]["PF00462"] == 2


def test_versioned_accession_is_stripped_and_merged():
    rows = [mkrow("1", "M" * 60, hits=[hit("PF00462.21", s_from=1, s_to=30)]),
            mkrow("2", "M" * 60, hits=[hit("PF00462", s_from=1, s_to=30)])]
    acc = M.consume(rows, M.blank())
    assert dict(acc["fam_seqs"]) == {"PF00462": 2}


def test_reversed_envelope_coordinates_are_swapped():
    """s_to < s_from would otherwise yield a negative covered length."""
    fwd = M.consume([mkrow("1", "M" * 60, hits=[hit("PF00462", s_from=10, s_to=40)])], M.blank())
    rev = M.consume([mkrow("1", "M" * 60, hits=[hit("PF00462", s_from=40, s_to=10)])], M.blank())
    assert rev["n_residues_pfam"] == fwd["n_residues_pfam"] == 31


def test_repeated_domains_counted_once_per_family_but_residues_unioned():
    """Two copies of one family on one sequence: 1 sequence, 2 domains, 1 distinct family."""
    row = mkrow("1", "M" * 100, hits=[hit("PF00462", s_from=1, s_to=20),
                                      hit("PF00462", s_from=15, s_to=40)])
    acc = M.consume([row], M.blank())
    assert acc["fam_seqs"]["PF00462"] == 1
    assert acc["fam_doms"]["PF00462"] == 2
    assert acc["fam_res"]["PF00462"] == 40      # union of 1-20 and 15-40, not 20+26
    assert acc["ndom_hist"][2] == 1
    assert acc["nfam_hist"][1] == 1


def test_overlapping_domains_of_different_families_not_double_counted():
    row = mkrow("1", "M" * 100, hits=[hit("PF00462", s_from=1, s_to=50),
                                      hit("PF03960", s_from=40, s_to=80)])
    acc = M.consume([row], M.blank())
    assert acc["n_residues_pfam"] == 80         # union, not 50 + 41
    assert acc["nfam_hist"][2] == 1
    assert acc["fam_res"]["PF00462"] == 50      # per-family totals stay independent


@pytest.mark.parametrize("s_to,expected_bin", [
    (4, 0),          # inclusive coords: covers 4 of 100 aa -> 0.04 -> bin 0
    (5, 1),          # 0.05 sits on the 0.05-0.10 boundary -> bin 1
    (50, 10),        # 0.50 -> bin 10 (0.50-0.55)
    (100, 19),       # 1.00 must clamp into the last bin, not overflow to 20
])
def test_coverage_histogram_binning(s_to, expected_bin):
    row = mkrow("1", "M" * 100, hits=[hit("PF00462", s_from=1, s_to=s_to)])
    acc = M.consume([row], M.blank())
    assert acc["cov_hist"][expected_bin] == 1
    assert sum(acc["cov_hist"].values()) == 1


# ---------------------------------------------------------------------------
# reduce
# ---------------------------------------------------------------------------
def test_combine_sums_counts_but_maxes_hmm_reach():
    a, b = M.blank(), M.blank()
    a["n_seqs"], b["n_seqs"] = 3, 4
    a["fam_seqs"]["PF00462"], b["fam_seqs"]["PF00462"] = 3, 4
    a["fam_hmm_max"] = {"PF00462": 55, "PF03960": 12}
    b["fam_hmm_max"] = {"PF00462": 20, "PF03960": 99}

    tot = M.combine([a, b])
    assert tot["n_seqs"] == 7
    assert tot["fam_seqs"]["PF00462"] == 7
    assert tot["fam_hmm_max"] == {"PF00462": 55, "PF03960": 99}   # max, not sum


def test_combine_of_one_part_is_that_part():
    acc = M.consume([mkrow("1", "M" * 60, hits=[hit("PF00462")])], M.blank())
    acc["fam_hmm_max"] = dict(acc["fam_hmm_max"])
    tot = M.combine([acc])
    assert tot["n_seqs"] == acc["n_seqs"]
    assert dict(tot["fam_res"]) == dict(acc["fam_res"])


# ---------------------------------------------------------------------------
# quantiles -- the manuscript quotes medians from both code paths
# ---------------------------------------------------------------------------
def test_quantile_paths_agree_on_the_same_data():
    vals = [10] * 10 + [20] * 10 + [30] * 10
    from_values = M.quant_values(vals)
    from_hist = M.quant_hist(Counter(vals))
    assert from_hist[0.5] == from_values[0.5] == 20
    assert from_hist[0.25] == from_values[0.25] == 10
    assert from_hist[0.75] == from_values[0.75] == 30


def test_quant_values_interpolates_and_handles_singletons():
    assert M.quant_values([1, 2, 3, 4])[0.5] == 2.5
    assert M.quant_values([7])[0.5] == 7


def test_quantiles_of_nothing_are_nan():
    assert M.quant_values([])[0.5] != M.quant_values([])[0.5]      # NaN
    assert M.quant_hist(Counter())[0.5] != M.quant_hist(Counter())[0.5]


# ---------------------------------------------------------------------------
# Pfam model lengths
# ---------------------------------------------------------------------------
def test_read_pfam_leng_from_hmm(tmp_path):
    p = tmp_path / "Pfam-A.hmm"
    p.write_text("HMMER3/f\nNAME  Thioredoxin\nACC   PF00085.22\nLENG  103\n"
                 "HMM\n//\nNAME  Glutaredoxin\nACC   PF00462.27\nLENG  60\n//\n")
    assert M.read_pfam_leng(str(p)) == {"PF00085": 103, "PF00462": 60}


def test_read_pfam_leng_from_dat_and_gzip(tmp_path):
    text = "#=GF ID Thioredoxin\n#=GF AC PF00085.22\n#=GF ML 103\n" \
           "#=GF ID Glutaredoxin\n#=GF AC PF00462.27\n#=GF ML 60\n"
    plain = tmp_path / "Pfam-A.hmm.dat"
    plain.write_text(text)
    gz = tmp_path / "Pfam-A.hmm.dat.gz"
    gz.write_bytes(gzip.compress(text.encode()))

    expected = {"PF00085": 103, "PF00462": 60}
    assert M.read_pfam_leng(str(plain)) == expected
    assert M.read_pfam_leng(str(gz)) == expected


# ---------------------------------------------------------------------------
# end to end
# ---------------------------------------------------------------------------
def test_report_writes_family_table_and_replayable_counters(tmp_path, capsys):
    p = tmp_path / "in.csv"
    write_csv(p, [
        mkrow("1", "M" * 100, full="true", csize=6,
              hits=[hit("PF00462", hmm_to=58, s_from=1, s_to=60)]),
        mkrow("2", "M" * 80, hits=[hit("PF00462", hmm_to=60, s_from=1, s_to=40),
                                   hit("PF03960", hmm_to=46, s_from=50, s_to=80)]),
        mkrow("3", "M" * 50),
    ])
    tot = scan(p)
    M.report(tot, str(tmp_path / "out"))
    capsys.readouterr()

    table = (tmp_path / "out_family_sizes.tsv").read_text().splitlines()
    assert table[0] == "pfam_acc\tn_sequences\tn_domains\tn_residues\tmax_hmm_to"
    rows = {ln.split("\t")[0]: ln.split("\t")[1:] for ln in table[1:]}
    assert rows["PF00462"] == ["2", "2", "100", "60"]      # 60 + 40 residues, max hmm_to 60
    assert rows["PF03960"] == ["1", "1", "31", "46"]
    assert table[1].startswith("PF00462"), "table must be sorted by n_sequences desc"

    counters = json.loads((tmp_path / "out_counters.json").read_text())
    assert counters["n_seqs"] == 3
    assert counters["n_with_pfam"] == 2
    assert counters["fam_hmm_max"]["PF00462"] == 60
