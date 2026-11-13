#!/usr/bin/env python3
"""Fetch the members' source contigs from the MGnify API (stdlib only; workstation, needs internet).

Per assembly (ERZ): its analysis (MGYA) of the member's pipeline version, then
- v5.0: api/v1/analyses/<MGYA>/contigs/<contig_name>, one call per contig;
- older, or a v5.0 contig the API does not find: the analysis's "Processed contigs with pCDS" FASTA, streamed once.
Each sequence must be as long as `contig_length` of the lookup, or it is reported as a mismatch and not written.
Resumable: contigs already in --out are skipped.

Usage: fetch_contigs.py members.tsv --out contigs.fasta [--threads 8]
FASTA header: >contig_id assembly_accession contig_name
"""

import argparse
import csv
import gzip
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

API = "https://www.ebi.ac.uk/metagenomics/api/v1"


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            err = e
        except (urllib.error.URLError, TimeoutError) as e:
            err = e
        time.sleep(5 * (i + 1))
    raise RuntimeError(f"{url}: {err}")


def analysis(erz, version):
    """MGYA of this assembly for the pipeline version, else None."""
    data = json.loads(get(f"{API}/assemblies/{erz}/analyses") or b'{"data": []}')["data"]
    hits = [a["id"] for a in data if a["attributes"].get("pipeline-version") == version]
    return hits[0] if hits else None


def pcds_fasta(mgya, wanted):
    """{contig_name: sequence} for the wanted names, streaming the pCDS FASTA."""
    data = json.loads(get(f"{API}/analyses/{mgya}/downloads?page_size=250"))["data"]
    urls = [d["links"]["self"] for d in data
            if d["attributes"]["description"]["label"] == "Processed contigs with pCDS"]
    if not urls:
        return {}
    found, name, seq = {}, None, []
    with urllib.request.urlopen(urls[0], timeout=600) as r, gzip.open(r, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name in wanted:
                    found[name] = "".join(seq)
                    if len(found) == len(wanted):
                        return found
                name, seq = line[1:].split()[0], []
            elif name in wanted:
                seq.append(line.strip())
    if name in wanted:
        found[name] = "".join(seq)
    return found


def fetch_assembly(erz, version, contigs):
    """contigs: {contig_name: (contig_id, length)} -> (records, problems)."""
    mgya = analysis(erz, version)
    if mgya is None:
        return [], [(erz, name, "no analysis of version " + version) for name in contigs]
    seqs = {}
    if version == "5.0":
        for name in contigs:
            body = get(f"{API}/analyses/{mgya}/contigs/{urllib.parse.quote(name, safe='')}")
            if body:
                seqs[name] = "".join(body.decode().splitlines()[1:])
    rest = {n for n in contigs if n not in seqs}
    if rest:
        seqs.update(pcds_fasta(mgya, rest))
    records, problems = [], []
    for name, (cid, length) in contigs.items():
        s = seqs.get(name)
        if s is None:
            problems.append((erz, name, f"not found ({mgya})"))
        elif len(s) != length:
            problems.append((erz, name, f"length {len(s)} != {length} ({mgya})"))
        else:
            records.append((cid, erz, name, s))
    return records, problems


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("members", help="contig_lookup.py members.tsv")
    p.add_argument("--out", required=True, help="contigs FASTA (appended; resumable)")
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args(argv)
    done = set()
    if Path(a.out).exists():
        done = {l[1:].split()[0] for l in open(a.out) if l.startswith(">")}
    todo = defaultdict(dict)
    for r in csv.DictReader(open(a.members), delimiter="\t"):
        if r["contig_id"] not in done:
            version = f"{float(r['pipeline_version']):.1f}"  # parquet float32: 4.0999… -> 4.1
            todo[(r["assembly_accession"], version)][r["contig_name"]] = (r["contig_id"], int(r["contig_length"]))
    print(f"# {len(done)} contigs already fetched; {sum(map(len, todo.values()))} to fetch "
          f"from {len(todo)} assemblies", flush=True)
    n_ok, problems = 0, []
    with ThreadPoolExecutor(a.threads) as pool, open(a.out, "a") as out:
        jobs = {pool.submit(fetch_assembly, erz, v, c): erz for (erz, v), c in todo.items()}
        for i, job in enumerate(as_completed(jobs), 1):
            try:
                records, probs = job.result()
            except Exception as e:  # one failing assembly must not lose the others' contigs
                records, probs = [], [(jobs[job], "*", f"error: {e}")]
            for cid, erz, name, s in records:
                out.write(f">{cid} {erz} {name}\n{s}\n")
            out.flush()
            n_ok += len(records)
            problems += probs
            if i % 25 == 0 or i == len(jobs):
                print(f"#   {i}/{len(jobs)} assemblies, {n_ok} contigs, {len(problems)} problems", flush=True)
    for erz, name, why in problems:
        print(f"PROBLEM\t{erz}\t{name}\t{why}")
    print(f"# fetched {n_ok}, problems {len(problems)}; total in {a.out}: {len(done) + n_ok}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
