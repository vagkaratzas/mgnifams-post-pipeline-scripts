"""Annotation evidence per family (Results catalogue paragraph 6, Supplementary Table S1), release 1.0 and 1.1.

Usage: python3 annotation_evidence.py --old 1.0/output --new 1.1/output --reann 1.1/output_annotate/annotate_not_updated \
           --ids 1.1/output/update_families/families.tsv.gz --mask10 Pfam36.0.clans.tsv.gz --mask11 Pfam38.1.clans.tsv.gz

A family "has" a resource when its table_data CSV has a row for it (the DB's has_* flags, finalize_db.sql):
mgnifam_pfams (sequence level), mgnifam_model_pfams (hhsearch, profile level), mgnifam_funfams, mgnifam_folds
(Foldseek; AlphaFoldDB when the target starts with "AF-", else PDB). The 1.0 block must reproduce the published
17,576 / 7,169 / 9,570 / 4,290 / 1,174 / 13,399 / 25,988 / 26,416 / 7,362 / 7,266 / 14,628; a grid of hhsearch
probability and Foldseek E-value cutoffs shows which rule did, if the unfiltered one doesn't.
1.1: updated families from --new, not-updated ones from --reann, model Pfams from --old for all (seed MSAs and
HMMs unchanged, ANNOTATE_MODELS not rerun). "Late" Pfam = every hit accession absent from the masking release
(1.0 and the 1.1 not-updated: --mask10; 1.1 updated: --mask11). "DUF only" (Discussion paragraph 1) = every Pfam
accession matched at sequence or profile level is named DUF* or UPF* (names from --names, first file wins; default
--mask11; 1.1 sequence hits are Pfam-A 38.2, so pass the 38.2 list first), and no FunFam;
reported with and without a structure match. Read-only.
"""
import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)
PROBS = (None, 50.0, 90.0, 95.0)
EVALUES = (None, 0.01, 0.001)


def rows(path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh)


def accession(raw):
    return raw.strip().split(".")[0]


def read_mask(path):
    """Pfam accessions of a release: Pfam-A.clans.tsv.gz (first column) or Pfam-A.hmm.dat.gz (#=GF AC)."""
    out = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                out.add(accession(line.split()[2]))
            elif line.startswith("PF"):
                out.add(accession(line.split("\t")[0]))
    return out


def read_names(path):
    """Pfam accession -> short name (pfamA_id, column 4 of Pfam-A.clans.tsv.gz)."""
    with gzip.open(path, "rt") as fh:
        return {accession(f[0]): f[3] for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) > 3}


def unknown(name):
    return name.startswith(("DUF", "UPF"))


def fold_db(name):
    return "afdb" if name.startswith("AF-") else "pdb"


def fold_shape(name):
    return re.sub(r"[0-9]", "9", re.sub(r"[A-Za-z]", "a", name))[:12]


def load(table_dir, keep):
    """Per-family hit sets of one table_data folder, restricted to the family IDs in keep."""
    d = defaultdict(set)
    acc = defaultdict(set)
    for r in rows(f"{table_dir}/mgnifam_pfams.csv"):
        if r["id"] in keep:
            d["pfam"].add(r["id"])
            acc[r["id"]].add(accession(r["pfam"]))
    for r in rows(f"{table_dir}/mgnifam_funfams.csv"):
        if r["id"] in keep:
            d["funfam"].add(r["id"])
    folds = defaultdict(lambda: defaultdict(lambda: float("inf")))  # id -> db -> best E-value
    shapes = Counter()
    for r in rows(f"{table_dir}/mgnifam_folds.csv"):
        if r["id"] in keep:
            db = fold_db(r["fold"])
            shapes[(db, fold_shape(r["fold"]))] += 1
            folds[r["id"]][db] = min(folds[r["id"]][db], float(r["e_value"]))
    return d, acc, folds, shapes


def load_models(table_dir):
    best = defaultdict(float)  # id -> best hhsearch probability
    acc = defaultdict(set)
    for r in rows(f"{table_dir}/mgnifam_model_pfams.csv"):
        best[r["id"]] = max(best[r["id"]], float(r["prob"]))
        acc[r["id"]].add(accession(r["pfam"]))
    return best, acc


def evidence(ids, d, models, folds, prob=None, evalue=None):
    """Per-family evidence sets under an hhsearch probability floor and a Foldseek E-value cap."""
    model = {i for i in ids if i in models and (prob is None or models[i] >= prob)}
    struct = {db: {i for i in ids if db in folds.get(i, {}) and (evalue is None or folds[i][db] <= evalue)}
              for db in ("pdb", "afdb")}
    pfam, funfam = d["pfam"] & ids, d["funfam"] & ids
    family = pfam | model | funfam
    return dict(pfam=pfam, model=model, funfam=funfam, family=family, pdb=struct["pdb"], afdb=struct["afdb"],
                structure=struct["pdb"] | struct["afdb"])


def pct(x, n):
    return f"{x} ({100 * x / n:.1f}%)"


def report(name, ids, e, late):
    n = len(ids)
    any_ = e["family"] | e["structure"]
    afdb_only_struct = e["afdb"] - e["pdb"]
    afdb_only = e["afdb"] - e["pdb"] - e["family"]
    print(f"== {name}: {n} families")
    print(f"family resource {pct(len(e['family']), n)}: Pfam {len(e['pfam'])}, profile only "
          f"{pct(len(e['model'] - e['pfam']), n)}, FunFam {len(e['funfam'])}, FunFam only {len(e['funfam'] - e['pfam'] - e['model'])}")
    for label, (with_pfam, only_late) in late.items():
        print(f"  late Pfam ({label}): {only_late} of {with_pfam} with Pfam hit only families absent from the masking release"
              f" ({100 * only_late / with_pfam:.1f}%)" if with_pfam else f"  late Pfam ({label}): no Pfam hits")
    print(f"structure: PDB {len(e['pdb'])}, AlphaFoldDB {len(e['afdb'])}, either {pct(len(e['structure']), n)}; "
          f"AlphaFoldDB not PDB {len(afdb_only_struct)}; AlphaFoldDB only evidence {len(afdb_only)}")
    print(f"no match at all {pct(n - len(any_), n)}; no Pfam (seq or profile), FunFam or PDB "
          f"{pct(n - len(e['family'] | e['pdb']), n)}; no Pfam (seq), FunFam or PDB {pct(n - len(e['pfam'] | e['funfam'] | e['pdb']), n)}")


def grid(name, ids, d, models, folds):
    print(f"-- {name} grid: prob floor / E cap -> profile only, PDB, AlphaFoldDB, either, AFDB only evidence, none")
    for p in PROBS:
        for ev in EVALUES:
            e = evidence(ids, d, models, folds, p, ev)
            none = len(ids) - len(e["family"] | e["structure"])
            print(f"  prob {p} E {ev}: {len(e['model'] - e['pfam'])} {len(e['pdb'])} {len(e['afdb'])} "
                  f"{len(e['structure'])} {len(e['afdb'] - e['pdb'] - e['family'])} {none}")


def duf_only(e, acc, macc, names):
    """Families whose every Pfam match (sequence or profile) is DUF/UPF, without FunFam; accessions without a name."""
    out = set()
    for i in e["pfam"] | e["model"]:
        pf = acc.get(i, set()) | macc.get(i, set())
        if pf and all(unknown(names.get(a, "")) for a in pf) and i not in e["funfam"]:
            out.add(i)
    unnamed = {i: {a for a in acc.get(i, set()) | macc.get(i, set()) if a not in names} for i in e["pfam"] | e["model"]}
    return out, {i: u for i, u in unnamed.items() if u}


def duf_report(name, ids, e, acc, macc, names):
    n = len(ids)
    duf, unnamed = duf_only(e, acc, macc, names)
    print(f"DUF only ({name}): {pct(len(duf), n)} with only DUF/UPF Pfam matches and no FunFam; of them "
          f"no structure match {pct(len(duf - e['structure']), n)}, PDB {len(duf & e['pdb'])}, AlphaFoldDB {len(duf & e['afdb'])}; "
          f"seq-level DUF only {len({i for i in duf if acc.get(i)})}; matched accessions without a name "
          f"{len(set().union(*unnamed.values()))}, on {len(unnamed)} families")


def late_counts(ids, acc, mask):
    with_pfam = [i for i in ids if acc.get(i)]
    return len(with_pfam), sum(1 for i in with_pfam if not (acc[i] & mask))


def selftest():
    assert fold_db("AF-A0A1-F1-model_v4") == "afdb" and fold_db("2k4x.cif.gz_MODEL_17_A") == "pdb"
    assert accession("PF00005.30") == "PF00005"
    ids = {"1", "2", "3"}
    d = {"pfam": {"1"}, "funfam": set()}
    folds = {"2": {"afdb": 0.5}, "3": {"pdb": 1e-5, "afdb": 1e-9}}
    e = evidence(ids, d, {"3": 40.0}, folds)
    assert e["family"] == {"1", "3"} and e["afdb"] == {"2", "3"} and e["pdb"] == {"3"}
    e = evidence(ids, d, {"3": 40.0}, folds, prob=50.0, evalue=0.01)
    assert e["family"] == {"1"} and e["afdb"] == {"3"}
    assert late_counts(ids, {"1": {"PF1", "PF2"}, "2": {"PF3"}}, {"PF1"}) == (2, 1)
    e = dict(pfam={"1", "2"}, model={"3"}, funfam={"2"})
    names = {"PF1": "DUF285", "PF2": "ABC_tran", "PF3": "UPF0118"}
    assert duf_only(e, {"1": {"PF1"}, "2": {"PF1"}}, {"3": {"PF3"}}, names) == ({"1", "3"}, {})
    assert duf_only(e, {"1": {"PF1", "PF2"}}, {"3": {"PF9"}}, names) == (set(), {"3": {"PF9"}})


def main():
    selftest()
    p = argparse.ArgumentParser()
    for k in ("old", "new", "reann", "ids", "mask10", "mask11"):
        p.add_argument(f"--{k}", required=True)
    p.add_argument("--names", nargs="+", help="Pfam-A.clans.tsv.gz files for DUF names, first wins (default --mask11)")
    a = p.parse_args()
    with gzip.open(a.ids, "rt") as fh:
        status = {r["family_id"]: r["status"] for r in csv.DictReader(fh, delimiter="\t")}
    all_ids = set(status)
    upd = {i for i, s in status.items() if s == "updated"}
    nup = all_ids - upd
    mask10, mask11 = read_mask(a.mask10), read_mask(a.mask11)
    print(f"families {len(all_ids)}: updated {len(upd)}, not updated {len(nup)}; mask10 {len(mask10)}, mask11 {len(mask11)} accessions")
    models, macc = load_models(f"{a.old}/table_data")
    names = {}
    for path in reversed(a.names or [a.mask11]):
        names.update(read_names(path))

    d10, acc10, f10, s10 = load(f"{a.old}/table_data", all_ids)
    e10 = evidence(all_ids, d10, models, f10)
    report("1.0", all_ids, e10, {"36.0": late_counts(all_ids, acc10, mask10)})
    duf_report("1.0", all_ids, e10, acc10, macc, names)
    grid("1.0", all_ids, d10, models, f10)

    du, accu, fu, su = load(f"{a.new}/table_data", upd)
    dn, accn, fn, sn = load(f"{a.reann}/table_data", nup)
    d11 = {k: du[k] | dn[k] for k in ("pfam", "funfam")}
    f11 = {**fu, **fn}
    late = {"updated vs 38.1": late_counts(upd, accu, mask11), "not updated vs 36.0": late_counts(nup, accn, mask10)}
    e11 = evidence(all_ids, d11, models, f11)
    report("1.1", all_ids, e11, late)
    duf_report("1.1", all_ids, e11, {**accu, **accn}, macc, names)
    grid("1.1", all_ids, d11, models, f11)
    for name, s in (("1.0", s10), ("1.1 updated", su), ("1.1 not updated", sn)):
        print(f"fold name shapes {name}: " + ", ".join(f"{db}:{sh} {c}" for (db, sh), c in s.most_common(6)))


if __name__ == "__main__":
    main()
