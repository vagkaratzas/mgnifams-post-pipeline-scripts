#!/usr/bin/env python3
"""Do family representatives resolve in the MGnify Proteins API?

Usage: python -I rep_resolve.py FAMILIES.tsv.gz OUT.tsv [--sample N] [--seed 1]   (or --selftest)
Checks every not_updated family's representative (they keep their release 1.0 MGYP) and a random sample of N
updated ones. Per family: status ok (sequence returned, rep_region inside it) / short (region beyond the
sequence) / missing (HTTP 404) / error. Read-only GETs, 4 at a time.
"""
import csv
import gzip
import json
import random
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

API = "https://www.ebi.ac.uk/metagenomics/proteins/api/v1/protein/"


def region_end(region):
    return 0 if region == "-" else int(region.split("-")[1])


def classify(seq, region, http=None):
    if http == 404:
        return "missing"
    if http is not None or seq is None:
        return "error"
    return "ok" if region_end(region) <= len(seq) else "short"


def get(mgyp, tries=4):
    """(sequence, HTTP status or None)."""
    for i in range(tries):
        try:
            with urllib.request.urlopen(API + mgyp, timeout=60) as r:
                return json.load(r).get("sequence"), None
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, 404
            err = e.code
        except (urllib.error.URLError, TimeoutError):
            err = -1
        time.sleep(5 * (i + 1))
    return None, err


def check(row):
    seq, http = get(row["rep_id"])
    return {**row, "length": len(seq) if seq else "", "status": classify(seq, row["rep_region"], http)}


def main(families, out, sample=200, seed=1):
    with gzip.open(families, "rt") as fh:
        rows = [dict(family_id=r["family_id"], group=r["status"], rep_id=r["rep_id"], rep_region=r["rep_region"])
                for r in csv.DictReader(fh, delimiter="\t")]
    nup = [r for r in rows if r["group"] == "not_updated"]
    upd = random.Random(seed).sample([r for r in rows if r["group"] == "updated"], sample)
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(check, nup + upd))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0]), delimiter="\t")
        w.writeheader()
        w.writerows(res)
    for grp in ("not_updated", "updated"):
        c = Counter(r["status"] for r in res if r["group"] == grp)
        print(f"{grp}: {sum(c.values())} checked: {dict(sorted(c.items()))}")
    print(f"distinct representative MGYPs checked: {len({r['rep_id'] for r in res})}; wrote {out}")


def selftest():
    assert classify("A" * 10, "3-10") == "ok" and classify("A" * 10, "3-11") == "short"
    assert classify("A" * 10, "-") == "ok" and classify(None, "1-5", 404) == "missing"
    assert classify(None, "1-5", 500) == "error"
    print("selftest ok")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["--selftest"]:
        selftest()
    else:
        kw = {"sample": int(a[a.index("--sample") + 1]) if "--sample" in a else 200,
              "seed": int(a[a.index("--seed") + 1]) if "--seed" in a else 1}
        main(a[0], a[1], **kw)
