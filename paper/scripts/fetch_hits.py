#!/usr/bin/env python3
"""Fetch the UniProtKB hit proteins of uniprot_hit_taxa.tsv (column 3, sp|ACC|… / tr|ACC|…) as FASTA (stdlib).

Usage: fetch_hits.py uniprot_hit_taxa.tsv --out hits.fasta
stdout: the UniProt release served (the hits came from 2026_02), and any accession not returned.
"""

import argparse
import csv
import json
import urllib.request

URL = "https://rest.uniprot.org/uniprotkb/accessions?format=fasta&accessions="
UNISAVE = "https://rest.uniprot.org/unisave/"


def unisave(acc):
    """Latest archived entry of an accession no longer in UniProtKB (FASTA, header gains its UniSave release)."""
    with urllib.request.urlopen(f"{UNISAVE}{acc}?format=json", timeout=120) as r:
        latest = json.load(r)["results"][0]
    with urllib.request.urlopen(f"{UNISAVE}{acc}?format=fasta&versions={latest['entryVersion']}", timeout=120) as r:
        return r.read().decode(), latest["lastRelease"]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("taxa")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    accs = sorted({row[2].split("|")[1] for row in csv.reader(open(a.taxa), delimiter="\t")})
    fasta, release = [], set()
    for i in range(0, len(accs), 100):
        with urllib.request.urlopen(URL + ",".join(accs[i:i + 100]), timeout=120) as r:
            release.add(r.headers.get("X-UniProt-Release"))
            fasta.append(r.read().decode())
    text = "".join(fasta)
    got = {line.split("|")[1] for line in text.splitlines() if line.startswith(">")}
    print(f"# {len(accs)} accessions, {len(got)} returned, UniProt release {sorted(release)}")
    for acc in sorted(set(accs) - got):  # retired since the 2026_02 search
        record, last = unisave(acc)
        text += record if record.endswith("\n") else record + "\n"
        print(f"UNISAVE\t{acc}\tlast in {last}")
    open(a.out, "w").write(text)
    print(f"# {text.count('>')} sequences -> {a.out}")


if __name__ == "__main__":
    main()
