import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin" / "ida2synteny.py"
spec = importlib.util.spec_from_file_location("ida2synteny", BIN)
i2s = importlib.util.module_from_spec(spec)
# dataclasses resolve their annotations through sys.modules, so register first
sys.modules[spec.name] = i2s
spec.loader.exec_module(i2s)


def domtbl_row(target, hmm_from, hmm_to, qlen=100, score=200.0, evalue=1e-30):
    """A --domtblout line; only the fields the parser reads need to be real."""
    f = ["-"] * 23
    f[0], f[2], f[3] = target, "500", "query"
    f[5] = str(qlen)
    f[6], f[7] = str(evalue), str(score)
    f[15], f[16] = str(hmm_from), str(hmm_to)
    f[22] = "some description with spaces"
    return " ".join(f)


def test_parse_domtbl_merges_domain_envelopes(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text("\n".join([
        "# a comment line",
        domtbl_row("two_domains", 1, 30),      # 30 residues
        domtbl_row("two_domains", 61, 100),    # + 40, disjoint -> 70/100
        domtbl_row("overlapping", 1, 50),
        domtbl_row("overlapping", 40, 60),     # union is 1..60 -> 60/100
        domtbl_row("nested", 1, 80),
        domtbl_row("nested", 20, 40),          # inside the first -> 80/100
    ]) + "\n")

    hits = i2s.parse_domtbl(path)
    assert hits["two_domains"]["coverage"] == 0.7
    assert hits["overlapping"]["coverage"] == 0.6
    assert hits["nested"]["coverage"] == 0.8
    # the full-sequence values repeat per row; they must not be double-counted
    assert hits["nested"]["bitscore"] == 200.0
    assert hits["nested"]["evalue"] == 1e-30


def test_parse_domtbl_skips_short_and_blank_lines(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text("\n\n# header\ntoo few fields here\n" + domtbl_row("ok", 1, 100) + "\n")
    assert list(i2s.parse_domtbl(path)) == ["ok"]


def make_row(rank, strand, anchor_strand, protein_id, hmm_hit=False):
    return i2s.NeighbourRow(
        ida_id="x", ida_string="x", uniprot_acc="x", organism_name="x", tax_id="1",
        proteome="", assembly="", contig_acc="c", contig_len_bp=1000,
        n_genes_on_contig=10, anchor_protein_id="A", anchor_start=1, anchor_end=2,
        anchor_strand=anchor_strand, rank=rank, protein_id=protein_id, locus_tag="",
        product="", start=1, end=2, strand=strand, length_aa=1,
        intergenic_gap_bp=None, same_strand_as_anchor=(strand == anchor_strand),
        hmm_hit=hmm_hit, hmm_evalue=None, hmm_bitscore=None, hmm_coverage=None,
        dist_to_contig_start_bp=0, dist_to_contig_end_bp=0,
        window_truncated=False, fusion_candidate=False,
    )


def test_sketch_is_in_the_anchors_reading_direction():
    # anchor on the minus strand: a neighbour on the minus strand reads forward
    rows = [make_row(-1, -1, -1, "UP"), make_row(0, -1, -1, "ANCHOR"),
            make_row(1, 1, -1, "DOWN", hmm_hit=True)]
    assert i2s.sketch_neighbourhood(rows) == "UP->  [ANCHOR]->  <-*DOWN*"


def test_wgs_master_re():
    assert i2s.WGS_MASTER_RE.match("WPOC01000000")
    assert not i2s.WGS_MASTER_RE.match("WPOC01000026")   # a real contig
    assert not i2s.WGS_MASTER_RE.match("CP146256")       # a complete genome


def test_ncbi_params_picks_up_credentials(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    bare = i2s.ncbi_params("WPOC01000026")
    assert bare["id"] == "WPOC01000026" and bare["rettype"] == "gbwithparts"
    assert "api_key" not in bare and "email" not in bare

    monkeypatch.setenv("NCBI_API_KEY", "KEY")
    monkeypatch.setenv("NCBI_EMAIL", "a@b.c")
    assert i2s.ncbi_params("X")["api_key"] == "KEY"
    assert i2s.ncbi_params("X")["email"] == "a@b.c"


def test_parse_contig_reads_genbank_and_embl(tmp_path):
    """Both flatfile dialects must give the same genes and contig length."""
    gb = tmp_path / "T00000001.gb"
    gb.write_text(
        "LOCUS       T00000001               600 bp    DNA     linear BCT\n"
        "FEATURES             Location/Qualifiers\n"
        "     CDS             10..30\n"
        '                     /protein_id="AAA00001.1"\n'
        '                     /product="first"\n'
        '                     /translation="MKV"\n'
        "     CDS             complement(100..120)\n"
        '                     /protein_id="AAA00002.1"\n'
        '                     /product="second"\n'
        '                     /translation="MTT"\n'
        "ORIGIN\n//\n"
    )
    genes, length = i2s.parse_contig(gb, "T00000001")
    assert length == 600
    assert [g.protein_id for g in genes] == ["AAA00001", "AAA00002"]
    assert [g.start for g in genes] == [10, 100]
    assert [g.strand for g in genes] == [1, -1]
    assert [g.index for g in genes] == [0, 1]
    assert genes[0].length_aa == 3
