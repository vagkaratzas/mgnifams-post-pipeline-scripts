#!/usr/bin/env python3
"""Self-check: union coverage, per-seq annotated flag, and exclusive (subtract) math."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import calculate_stats_from_domtbl as m  # noqa: E402


def _domtbl_line(target, ali_from, ali_to, i_evalue=1e-10):
    # 23-col domtblout: target ... i-Evalue(idx12) ... ali_from(idx17) ali_to(idx18) ... desc
    cols = ["-"] * 23
    cols[0] = target
    cols[12] = str(i_evalue)
    cols[17] = str(ali_from)
    cols[18] = str(ali_to)
    return " ".join(cols) + "\n"


def _write(tmp, name, lines):
    p = tmp / name
    p.write_text("".join(lines))
    return str(p)


def test_stats(tmp_path=Path("/tmp/_stats_selfcheck")):
    tmp_path.mkdir(parents=True, exist_ok=True)
    lengths = {"seqA": 100, "seqB": 50, "seqC": 30}  # seqC never hit

    # pfam: seqA 1-40, seqB 1-50
    pfam = _write(tmp_path, "pfam.domtbl", [_domtbl_line("seqA", 1, 40), _domtbl_line("seqB", 1, 50)])
    # mgnifam: seqA 30-60 (overlaps pfam 30-40), seqB 1-50 (fully covered by pfam already)
    mgnifam = _write(
        tmp_path, "mgnifam.domtbl", [_domtbl_line("seqA", 30, 60), _domtbl_line("seqB", 1, 50)]
    )

    # union (pfam + mgnifam): seqA = 1-60 = 60 residues, seqB = 50 -> 110
    ranges = m.parse_domtbls([pfam, mgnifam])
    union = m.calculate_stats(lengths, ranges, "pfam_mgnifam")
    assert union["total_sequences"] == 3
    assert union["annotated_sequences"] == 2
    assert union["total_amino_acids"] == 180
    assert union["annotated_amino_acids"] == 110, union["annotated_amino_acids"]

    # mgnifam-only with pfam subtracted -> exclusive residues.
    # seqA mgnifam 30-60 (31 res) minus pfam 1-40 -> 41-60 = 20 exclusive
    # seqB mgnifam 1-50 minus pfam 1-50 -> 0 exclusive
    mg_ranges = m.parse_domtbls([mgnifam])
    pf_ranges = m.parse_domtbls([pfam])
    excl = m.calculate_stats(lengths, mg_ranges, "mgnifam", subtract_ranges=pf_ranges)
    assert excl["exclusive_annotated_amino_acids"] == 20, excl["exclusive_annotated_amino_acids"]
    assert excl["exclusive_annotated_sequences"] == 1

    # evalue threshold filters weak hits
    weak = _write(tmp_path, "weak.domtbl", [_domtbl_line("seqC", 1, 30, i_evalue=1.0)])
    filtered = m.parse_domtbls([weak], domain_evalue_threshold=0.001)
    assert "seqC" not in filtered

    print("OK")


if __name__ == "__main__":
    test_stats()
