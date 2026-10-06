import gzip
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "pfam_version_checks.py"


def run(argv):
    return subprocess.run([sys.executable, str(SCRIPT)] + argv, capture_output=True, text=True)


def write_clans(path, accessions):
    with gzip.open(path, "wt") as fh:
        for acc in accessions:
            fh.write("%s\t\tname\tname\tdescription\n" % acc)


def write_input(path, hits_per_protein):
    """MGnify protein CSV; hits are the "p" entries: [accession, E, score, ...]."""
    with gzip.open(path, "wt") as fh:
        fh.write("mgyp,sequence,full_length,cluster_size,metadata\n")
        for i, accs in enumerate(hits_per_protein):
            meta = json.dumps({"p": [[a, 1e-5, 30.0, 1, 50, 1, 50] for a in accs]})
            fh.write('%d,MKV,false,1,"%s"\n' % (i, meta.replace('"', '""')))


def test_infer_picks_the_release_that_contains_every_accession(tmp_path):
    old, new = tmp_path / "old.tsv.gz", tmp_path / "new.tsv.gz"
    write_clans(old, ["PF00001", "PF00002"])
    write_clans(new, ["PF00001", "PF00002", "PF09999"])
    data = tmp_path / "chunk.csv.gz"
    write_input(data, [["PF00001", "PF09999"], ["PF09999"], []])
    proc = run(["infer", "--input", str(data), "--release", "old=" + str(old),
                "--release", "new=" + str(new)])
    assert proc.returncode == 0, proc.stderr
    lines = [l.split("\t") for l in proc.stdout.splitlines()]
    assert lines[0] == ["release", "families", "accessions_missing", "hits_missing",
                        "distinct_accessions", "hits"]
    rows = {l[0]: dict(zip(lines[0], l)) for l in lines[1:]}
    # the sample holds 2 distinct accessions and 3 hits; PF09999 (2 hits) is new
    assert (rows["old"]["accessions_missing"], rows["old"]["hits_missing"]) == ("1", "2")
    assert (rows["new"]["accessions_missing"], rows["new"]["hits_missing"]) == ("0", "0")
    assert rows["new"]["distinct_accessions"] == "2" and rows["new"]["hits"] == "3"


def test_late_hits_counts_families_whose_pfam_hits_are_all_newer(tmp_path):
    base = tmp_path / "pfam36.tsv.gz"
    write_clans(base, ["PF00001", "PF00002"])
    hits = tmp_path / "mgnifam_pfams.csv"
    hits.write_text(
        "id,pfam,name,e_value\n"
        "1,PF00001,a,1e-5\n"            # family 1: only an old family
        "2,PF09999,b,1e-5\n"            # family 2: only a post-baseline family
        "3,PF00002,c,1e-5\n"
        "3,PF08888,d,1e-5\n"            # family 3: both
    )
    out = tmp_path / "late.tsv"
    proc = run(["late-hits", "--pfams", str(hits), "--baseline", str(base), "--output", str(out)])
    assert proc.returncode == 0, proc.stderr
    assert "families with a Pfam hit: 3" in proc.stdout
    assert "with a hit to a family newer than the baseline: 2" in proc.stdout
    assert "only newer-family hits: 1" in proc.stdout
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    assert rows[0] == ["id", "n_hits", "n_newer", "newer_accessions"]
    assert {r[0]: r[1:] for r in rows[1:]} == {
        "1": ["1", "0", ""], "2": ["1", "1", "PF09999"], "3": ["2", "1", "PF08888"]}
