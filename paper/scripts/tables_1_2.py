"""Tables 1 and 2 of the MGnifams paper from a coverage run's report.

Reads the pipeline's coverage_report.md and library_sizes.tsv (ali coordinates)
and writes the two Drive-ready TSVs, formatted as in the doc.

    python paper/scripts/tables_1_2.py
    python paper/scripts/tables_1_2.py --selftest
"""
import argparse
import re
import sys

DATA = "paper/data"
DBS = ["SwissProt", "TrEMBL", "UniProtKB"]
# report category -> Table 2 row label, in paper order
CATS = [("membrane_a", "membrane-α"), ("membrane_b", "membrane-β"),
        ("tm_any", "transmembrane (any)"), ("disorder", "disordered"),
        ("tm_or_disorder", "transmembrane or disordered"), ("not_tm_disorder", "neither"),
        ("novel_structure_any", "no Pfam or FunFam match"),
        ("novel", "no annotation from any source"), ("any", "all families")]


def num(s):
    return int(s.replace(",", ""))


def md_table(text, heading):
    """Rows (lists of cells) of the first markdown table after a heading starting with `heading`."""
    block = text.split("\n## " + heading, 1)[1].split("\n## ", 1)[0]
    rows = [l.strip().strip("|").split("|") for l in block.splitlines() if l.startswith("|")]
    return [[c.strip() for c in r] for r in rows[2:]]  # drop header and separator


def parse(report):
    t1 = {r[0]: r for r in md_table(report, "Table 1.")}
    t2 = {(r[0], r[1]): r for r in md_table(report, "Table 2.")}
    pfam_seqs = {db: num(n) for db, n in re.findall(
        r"in (\w+), Pfam covered.*?At the sequence level, ([\d,]+) entries", report)}
    mg_seqs = {db: int(n) for db, n in re.findall(
        r"PASS  (\w+)/ali/any exclusive targets <= total  -- \d+ vs (\d+)", report)}
    return t1, t2, pfam_seqs, mg_seqs


def cell(n, total):
    return f"{n:,} ({100 * n / total:.2f})"


def table1(t1, pfam_seqs, mg_seqs):
    out = [["database", "level", "total", "Pfam", "MGnifams", "union", "MGnifams only"]]
    for db in DBS:
        r = t1[db]
        seqs, res, pf_res, mg_res, only_res, only_seqs = (num(r[i]) for i in (1, 2, 3, 5, 7, 10))
        pf_seqs, m_seqs = pfam_seqs[db], mg_seqs[db.lower()]
        out.append([db, "sequence", f"{seqs:,}", cell(pf_seqs, seqs), cell(m_seqs, seqs),
                    cell(pf_seqs + only_seqs, seqs), cell(only_seqs, seqs)])
        out.append([db, "residue", f"{res:,}", cell(pf_res, res), cell(mg_res, res),
                    cell(pf_res + only_res, res), cell(only_res, res)])
    return out


def table2(t2, sizes, t1):
    library = sizes["__library__"]
    db_res = num(t1["UniProtKB"][2])
    excl = num(t2[("UniProtKB", "any")][4])
    out = [["category", "families", "% of library", "with UniProtKB hit", "MGnifams-only residues",
            "% of total", "% of UniProtKB residues", "fold enrichment", "newly covered sequences"]]
    for cat, label in CATS:
        r = t2[("UniProtKB", cat)]
        hit, res, newly = num(r[2]), num(r[4]), num(r[7])
        n = library if cat == "any" else sizes[cat]
        lib, share = 100 * n / library, 100 * res / excl
        if cat == "any":  # shares of themselves
            out.append([label, f"{n:,}", "100", f"{hit:,}", f"{res:,}", "100",
                        f"{100 * res / db_res:.2f}", "—", f"{newly:,}"])
            continue
        # fold from exact counts: the report's Table 3 divides a share already rounded to 2 decimals
        out.append([label, f"{n:,}", f"{lib:.2f}", f"{hit:,}", f"{res:,}", f"{share:.2f}",
                    f"{100 * res / db_res:.2f}", f"{share / lib:.2f}", f"{newly:,}"])
    return out


def read_sizes(path):
    with open(path) as f:
        next(f)
        return {k: int(v) for k, v in (l.split() for l in f if l.strip())}


def build(report_path, sizes_path):
    report = open(report_path).read()
    t1, t2, pfam_seqs, mg_seqs = parse(report)
    return table1(t1, pfam_seqs, mg_seqs), table2(t2, read_sizes(sizes_path), t1)


def write(rows, path):
    with open(path, "w") as f:
        f.writelines("\t".join(r) + "\n" for r in rows)


def selftest():
    a, b = build(f"{DATA}/coverage_report_1.1cats.md", f"{DATA}/library_sizes_1.1cats.tsv")
    # Table 1 as in the doc (published run, 2026-10-07)
    assert a[5] == ["UniProtKB", "sequence", "149,810,139", "117,233,257 (78.25)", "47,921,275 (31.99)",
                    "121,115,785 (80.85)", "3,882,528 (2.59)"], a[5]
    assert a[1][4] == "164,295 (28.55)" and a[4][5] == "34,660,567,341 (59.71)"
    # Table 2 on release 1.1 categories
    assert b[2] == ["membrane-β", "366", "1.03", "352", "113,195,167", "2.13", "0.19", "2.07", "84,456"], b[2]
    assert b[6][5] == "86.80" and b[7][7] == "0.73" and b[8][7] == "0.83"
    assert b[9] == ["all families", "35,459", "100", "34,144", "5,306,060,877", "100", "9.11", "—", "3,882,528"]
    print("selftest ok")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--report", default=f"{DATA}/coverage_report_1.1cats.md")
    p.add_argument("--sizes", default=f"{DATA}/library_sizes_1.1cats.tsv")
    p.add_argument("--out1", default="tables/table-1_uniprot-coverage.tsv")
    p.add_argument("--out2", default="tables/table-2_mgnifam-only-by-category.tsv")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest:
        return selftest()
    t1, t2 = build(a.report, a.sizes)
    write(t1, a.out1)
    write(t2, a.out2)
    print(a.out1, a.out2, sep="\n")


if __name__ == "__main__":
    sys.exit(main())
