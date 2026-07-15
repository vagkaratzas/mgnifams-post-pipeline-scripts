import importlib.util
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "synteny_census.py"


def load_module():
    """Import the standalone script without requiring `bin` to be a package."""
    spec = importlib.util.spec_from_file_location("synteny_census", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_module()


# ---------------------------------------------------------------------------
# pure helpers
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("raw,expected", [
    ("PF01183", "PF01183"),
    ("pf01183", "PF01183"),
    ("PF01183.7", "PF01183"),
    (1183, "PF01183"),
    ("1183", "PF01183"),
    (1183.0, "PF01183"),          # float leaks in via iterrows upcast; must still canonicalise
])
def test_pfam_key_canonicalises(raw, expected):
    assert M.pfam_key(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("+", 1), ("-", -1), ("1", 1), ("-1", -1),
    ("plus", 1), ("minus", -1), ("forward", 1), ("reverse", -1),
    ("?", 0), (None, 0),
])
def test_strand_to_int(raw, expected):
    assert M.strand_to_int(raw) == expected


def test_bp_gap():
    assert M.bp_gap(100, 200, 300, 400) == 100    # b downstream
    assert M.bp_gap(300, 400, 100, 200) == 100    # b upstream
    assert M.bp_gap(100, 300, 200, 400) == 0      # overlap


def test_aa_len():
    assert M.aa_len(1, 300) == 100
    assert M.aa_len(300, 1) == 0                   # guarded against negatives


def test_partner_keys_prefers_pfams():
    assert M.partner_keys({"pfams": {"PF25309", "PF01183"}, "cluster_rep": 900}) \
        == ["PF01183", "PF25309"]                  # sorted
    assert M.partner_keys({"pfams": set(), "cluster_rep": 900}) == ["cluster:900"]


# ---------------------------------------------------------------------------
# anchor_window: signed offsets, AND semantics, edge statuses
# ---------------------------------------------------------------------------
def _gene(rank, start, end, strand=1, pfams=None, rep=1, clen=50_000):
    return {"rank": rank, "start": start, "end": end, "strand_i": strand,
            "pfams": pfams or set(), "cluster_rep": rep, "contig_length": clen,
            "protein_id": str(rank), "contig_id": 1}


def test_window_signed_offset_plus_strand():
    genes = [_gene(0, 15_000, 15_500, 1), _gene(1, 16_000, 16_600, 1),
             _gene(2, 17_000, 17_200, 1)]
    win, status, up, down = M.anchor_window(genes[0], genes, 5, 10_000)
    offs = {g["rank"]: g["offset"] for g in win}
    assert offs == {1: 1, 2: 2}                    # both downstream of anchor
    assert status == "FULL"


def test_window_offset_mirrors_on_minus_strand():
    # same physical layout, anchor on -; downstream in reading frame = lower coords
    genes = [_gene(0, 57_000, 57_200, -1), _gene(1, 58_500, 59_000, -1),
             _gene(2, 60_000, 60_500, -1)]
    win, _, _, _ = M.anchor_window(genes[2], genes, 5, 10_000)   # anchor is rank 2
    offs = {g["rank"]: g["offset"] for g in win}
    assert offs == {1: 1, 0: 2}                    # signed by anchor direction, not coords


def test_window_bp_cap_is_AND_not_or():
    # neighbour is 1 gene away by rank but 76 kb away by bp: must be excluded
    genes = [_gene(0, 1_000, 1_500, 1, clen=200_000),
             _gene(1, 77_000, 77_500, 1, clen=200_000)]
    win, _, _, _ = M.anchor_window(genes[0], genes, 5, 10_000)
    assert win == []


def test_window_status_partial_and_edge():
    # anchor near left edge, plenty of room right -> PARTIAL
    genes = [_gene(0, 200, 700, 1, clen=50_000), _gene(1, 5_000, 5_500, 1, clen=50_000)]
    _, status, up, down = M.anchor_window(genes[0], genes, 5, 10_000)
    assert status == "PARTIAL"
    # short contig, both sides truncated -> EDGE
    genes = [_gene(0, 200, 700, 1, clen=2_000), _gene(1, 1_200, 1_600, 1, clen=2_000)]
    _, status, _, _ = M.anchor_window(genes[0], genes, 5, 10_000)
    assert status == "EDGE"


def test_window_isolated_when_alone():
    genes = [_gene(0, 500, 900, 1, clen=1_200)]
    win, status, _, _ = M.anchor_window(genes[0], genes, 5, 10_000)
    assert win == [] and status == "ISOLATED"


# ---------------------------------------------------------------------------
# end-to-end run() on a synthetic parquet store  (the compute-saver regression)
# ---------------------------------------------------------------------------
def _store(tmp_path, pfam_extra_cols=True):
    """Synthetic metadata + pfam parquets mirroring the real schema.

    contig 1: anchor 100 (+); PF01183 at +1, unannotated cluster 900 at +2.
    contig 4: anchor 400 (-), independent lineage, SAME arrangement mirrored in coords.
    contig 2: anchor 200 (+), FULL window, no PF01183 -> a confident absence.
    contig 3: anchor 300 alone on a 1.2 kb fragment -> ISOLATED.
    """
    meta = pd.DataFrame([
        (1, 100, 15_000, 15_500, 1, 100, 50_000),
        (1, 101, 16_000, 16_600, 1, 101, 50_000),
        (1, 102, 17_000, 17_200, 1, 900, 50_000),
        (2, 200, 20_000, 20_900, 1, 200, 50_000),
        (2, 201, 30_500, 31_000, 1, 201, 50_000),
        (3, 300, 500, 900, 1, 300, 1_200),
        (4, 402, 57_000, 57_200, -1, 900, 90_000),
        (4, 401, 58_500, 59_000, -1, 401, 90_000),
        (4, 400, 60_000, 60_500, -1, 400, 90_000),
    ], columns=["contig_id", "protein_id", "start_position", "end_position",
                "strand", "cluster_rep", "contig_length"])
    meta["contig_name"] = "contig_" + meta["contig_id"].astype(str)
    rows = [(101, 1183), (201, 9999), (401, 1183)]
    if pfam_extra_cols:
        # float column present: reproduces the int/float row upcast that dropped Pfams
        pfam = pd.DataFrame([(p, a, 1e-9) for p, a in rows],
                            columns=["protein_id", "pfam_accession", "i_evalue"])
    else:
        pfam = pd.DataFrame(rows, columns=["protein_id", "pfam_accession"])

    meta_p = tmp_path / "m.parquet"
    pfam_p = tmp_path / "p.parquet"
    ids_p = tmp_path / "ids.txt"
    contigs_p = tmp_path / "c.txt"
    meta.to_parquet(meta_p)
    pfam.to_parquet(pfam_p)
    ids_p.write_text("MGYP100/1-176\nMGYP200\nMGYP300\nMGYP400\n")
    contigs_p.write_text("1\n2\n3\n4\n")
    argv = [str(ids_p), "--metadata", str(meta_p), "--pfam", str(pfam_p),
            "--id-strip-prefix", "MGYP", "--outdir", str(tmp_path / "out")]
    return argv, str(contigs_p)


def test_run_attaches_pfams_despite_float_columns(tmp_path):
    """Regression: a float column in the pfam parquet must not drop every Pfam."""
    argv, _ = _store(tmp_path, pfam_extra_cols=True)
    pdf, _ = M.run(M.parse_args(argv))
    p = pdf.set_index("partner")
    assert "PF01183" in p.index                    # would be absent under the old iterrows path
    assert (pdf["type"] == "pfam").sum() >= 1


def test_run_signed_offsets_and_conservation(tmp_path):
    argv, _ = _store(tmp_path)
    pdf, adf = M.run(M.parse_args(argv))
    p = pdf.set_index("partner")
    # PF01183 recurs on two independent lineages, downstream on both strands, always
    assert p.loc["PF01183", "n_indep"] == 2
    assert p.loc["PF01183", "offset"] == 1
    assert bool(p.loc["PF01183", "same_strand"])
    assert p.loc["PF01183", "frac_conserved"] == 1.0
    assert p.loc["PF01183", "frac_indep"] == round(2 / 3, 3)   # contig 2 = confident absence
    # unannotated neighbour nameable by cluster_rep, conserved at +2
    assert p.loc["cluster:900", "type"] == "unannotated_cluster"
    assert p.loc["cluster:900", "offset"] == 2
    # assessability is target-independent: 3 FULL + 1 ISOLATED
    assert sorted(adf["status"]) == ["FULL", "FULL", "FULL", "ISOLATED"]


def test_run_contigs_flag_matches_full_run(tmp_path):
    argv, contigs_p = _store(tmp_path)
    pdf, _ = M.run(M.parse_args(argv))
    pdf2, _ = M.run(M.parse_args(argv + ["--contigs", contigs_p]))
    assert pdf2.equals(pdf)


def test_run_targets_verdict(tmp_path):
    argv, _ = _store(tmp_path)
    _, adf = M.run(M.parse_args(argv + ["--targets", "PF01183"]))
    assert list(adf[adf.target_hit != ""]["contig_id"]) == [1, 4]


def test_run_pfam_without_extra_columns_still_works(tmp_path):
    """The all-int pfam schema (old self-test shape) must keep working too."""
    argv, _ = _store(tmp_path, pfam_extra_cols=False)
    pdf, _ = M.run(M.parse_args(argv))
    assert "PF01183" in pdf.set_index("partner").index


def test_run_isolated_only_reports_nothing(tmp_path):
    """One gene per contig -> all ISOLATED -> empty partner table, pdf.empty branch."""
    argv, _ = _store(tmp_path)
    # keep only the lone-anchor contig 3
    meta = pd.read_parquet(Path(argv[argv.index("--metadata") + 1]))
    meta[meta.contig_id == 3].to_parquet(argv[argv.index("--metadata") + 1])
    ids_p = Path(argv[0])
    ids_p.write_text("MGYP300\n")
    pdf, _ = M.run(M.parse_args(argv))
    assert pdf.empty


def test_run_threads_and_memory_limit(tmp_path):
    argv, _ = _store(tmp_path)
    pdf, _ = M.run(M.parse_args(argv + ["--threads", "1", "--memory-limit", "512MB"]))
    assert not pdf.empty


def test_run_warns_on_ids_missing_from_metadata(tmp_path, caplog):
    argv, _ = _store(tmp_path)
    Path(argv[0]).write_text("MGYP100\nMGYP200\nMGYP999\n")   # 999 absent from metadata
    with caplog.at_level("WARNING"):
        M.run(M.parse_args(argv))
    assert any("not in metadata" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# error-exit paths
# ---------------------------------------------------------------------------
def test_parse_args_requires_inputs():
    with pytest.raises(SystemExit):
        M.parse_args([])


def test_run_missing_file(tmp_path):
    argv, _ = _store(tmp_path)
    argv[argv.index("--metadata") + 1] = str(tmp_path / "nope.parquet")
    with pytest.raises(SystemExit, match="file not found"):
        M.run(M.parse_args(argv))


def test_run_metadata_missing_required_column(tmp_path):
    argv, _ = _store(tmp_path)
    meta_p = argv[argv.index("--metadata") + 1]
    pd.read_parquet(meta_p).drop(columns=["cluster_rep"]).to_parquet(meta_p)
    with pytest.raises(SystemExit, match="required column"):
        M.run(M.parse_args(argv))


def test_run_pfam_missing_columns(tmp_path):
    argv, _ = _store(tmp_path)
    pfam_p = argv[argv.index("--pfam") + 1]
    pd.DataFrame({"foo": [1], "bar": [2]}).to_parquet(pfam_p)
    with pytest.raises(SystemExit, match="pfam parquet"):
        M.run(M.parse_args(argv))


def test_run_no_usable_ids(tmp_path):
    argv, _ = _store(tmp_path)
    Path(argv[0]).write_text("MGYPabc\n")            # non-numeric, int column -> all dropped
    with pytest.raises(SystemExit, match="no usable protein_ids"):
        M.run(M.parse_args(argv))


def test_run_no_usable_contigs(tmp_path):
    argv, contigs_p = _store(tmp_path)
    Path(contigs_p).write_text("xyz\n")              # non-numeric, int column -> dropped
    with pytest.raises(SystemExit, match="no usable contig_ids"):
        M.run(M.parse_args(argv + ["--contigs", contigs_p]))


def test_run_ids_absent_from_metadata(tmp_path):
    argv, _ = _store(tmp_path)
    Path(argv[0]).write_text("MGYP777\n")            # numeric but not present
    with pytest.raises(SystemExit, match="none of the input ids"):
        M.run(M.parse_args(argv))


def test_run_unknown_contig_has_no_genes(tmp_path):
    argv, contigs_p = _store(tmp_path)
    Path(contigs_p).write_text("999999\n")           # valid int, no genes
    with pytest.raises(SystemExit, match="no genes found"):
        M.run(M.parse_args(argv + ["--contigs", contigs_p]))


def test_run_contigs_without_family_member(tmp_path):
    argv, contigs_p = _store(tmp_path)
    Path(argv[0]).write_text("MGYP100\n")            # anchor on contig 1 only
    Path(contigs_p).write_text("2\n")                # contig 2 carries no family member
    with pytest.raises(SystemExit, match="none of the input protein_ids are on"):
        M.run(M.parse_args(argv + ["--contigs", contigs_p]))


# ---------------------------------------------------------------------------
# small pure branches
# ---------------------------------------------------------------------------
def test_read_ids_skips_blanks_and_comments(tmp_path):
    f = tmp_path / "ids.txt"
    f.write_text("# comment\n\nMGYP1/1-9\nMGYP1\nMGYP2\n")
    assert M.read_ids(f, strip_prefix="MGYP") == ["1", "2"]   # dedup + suffix strip


def test_find_col_returns_none():
    assert M.find_col({"foo": "foo"}, ["bar"]) is None


def test_sanitise_ids_int_column_drops_nonnumeric(caplog):
    with caplog.at_level("WARNING"):
        assert M.sanitise_ids_for_type(["abc", "12"], "BIGINT") == ["12"]
    assert any("non-numeric" in r.message for r in caplog.records)


def test_sanitise_ids_string_column_keeps_all():
    assert M.sanitise_ids_for_type(["abc", "12"], "VARCHAR") == ["abc", "12"]


def test_window_uses_gene_span_when_no_contig_length():
    genes = [_gene(0, 1_000, 1_500, 1, clen=None), _gene(1, 2_000, 2_500, 1, clen=None)]
    win, status, up, down = M.anchor_window(genes[0], genes, 5, 10_000)
    assert [g["offset"] for g in win] == [1]         # neighbour found via span fallback
    assert up == 0 and down == 1_000                 # edges from called-gene span


def test_rank_partners_counts_anchor_once_per_partner():
    """Two neighbours sharing a Pfam must count the anchor once for that partner."""
    win = [
        {"pfams": {"PF01183"}, "cluster_rep": 1, "gap_bp": 100, "start": 1, "end": 300,
         "offset": 1, "same_strand": True},
        {"pfams": {"PF01183"}, "cluster_rep": 2, "gap_bp": 200, "start": 1, "end": 300,
         "offset": 2, "same_strand": True},
    ]
    anchors = [{"rep": "A", "status": "FULL", "contig_id": 1, "window": win}]
    neigh = pd.DataFrame({"pfams": [{"PF01183"}], "cluster_rep": [1]})
    pdf, _ = M.rank_partners(anchors, neigh, set())
    row = pdf.set_index("partner").loc["PF01183"]
    assert row["n_anchors"] == 1                      # deduped, not 2


def test_main_entrypoint(monkeypatch, tmp_path):
    argv, _ = _store(tmp_path)
    monkeypatch.setattr("sys.argv", ["synteny_census.py"] + argv)
    M.main()                                          # covers main() wiring
