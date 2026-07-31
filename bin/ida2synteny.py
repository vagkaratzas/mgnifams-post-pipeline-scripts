#!/usr/bin/env python3
"""
ida2synteny.py
==============

Use a Pfam domain architecture as a pre-computed proxy for a homology search,
then recover the genomic context of every matching protein in isolate genomes
and test a query HMM (e.g. an MGnifam) against the flanking genes.

Rationale
---------
Running an HMM against every isolate genome does not scale. InterPro has
already run Pfam across all of UniProtKB, so querying by domain architecture
(IDA) is a homology search someone else has already paid for. That gives a
small, high-precision anchor set. Only then is it cheap to fetch contigs and
run the expensive HMM against the handful of neighbouring proteins.

Note that an IDA is a property of a *protein*, not of a contig. Domains that
sit on separate genes will never co-occur in one architecture -- except in
gene-fusion cases, which this script flags explicitly. Recovering gene-level
synteny is what stages 3-6 below are for.

Pipeline
--------
  1. InterPro  : Pfam list          -> domain architectures (ida_id)
  2. InterPro  : ida_id             -> UniProtKB accessions
  3. UniProt   : accession          -> ENA contig accession + protein_id
  4. ENA       : contig accession   -> EMBL flatfile with CDS features
  5. local     : extract +/-N gene neighbourhood around each anchor
  6. HMMER     : hmmsearch query HMM vs all neighbour proteins
  7. report    : joined TSVs + human-readable neighbourhood sketch

Usage
-----
  ./ida2synteny.py --pfams pfams.txt --hmm MGYF243.hmm -o results/

  # only architectures containing exactly these domains, in this order
  ./ida2synteny.py --pfams pfams.txt --hmm MGYF243.hmm -o results/ \
      --ordered --exact --window 5

Requirements
------------
  python >= 3.8, requests, biopython, and hmmsearch (HMMER3) on $PATH.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import requests
from Bio import SeqIO

# --------------------------------------------------------------------------
# Endpoints
# --------------------------------------------------------------------------

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
UNIPROT_API = "https://rest.uniprot.org"
ENA_API = "https://www.ebi.ac.uk/ena/browser/api"

USER_AGENT = "ida2synteny/1.0 (https://github.com/; bioinformatics pipeline)"

# WGS master records (contig number all zeros) carry no CDS features -- they
# only describe the set. Fetching one is a common and silent failure mode.
WGS_MASTER_RE = re.compile(r"^[A-Z]{4,6}\d{2}0{6,}$")


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------


@dataclass
class Anchor:
    """A UniProt protein whose Pfam architecture matched the query."""

    uniprot_acc: str
    ida_id: str
    ida_string: str
    organism_name: str = ""
    tax_id: str = ""
    proteome: str = ""
    contig_acc: str = ""
    protein_id: str = ""
    n_domains: int = 0
    is_fusion_candidate: bool = False


@dataclass
class Gene:
    """A CDS feature on a contig."""

    contig_acc: str
    protein_id: str
    locus_tag: str
    product: str
    start: int
    end: int
    strand: int
    translation: str
    index: int = -1  # positional index along the contig

    @property
    def length_aa(self) -> int:
        return len(self.translation)


@dataclass
class NeighbourRow:
    """One gene in the neighbourhood of one anchor -- a row in the report."""

    ida_id: str
    ida_string: str
    uniprot_acc: str
    organism_name: str
    tax_id: str
    proteome: str
    assembly: str
    contig_acc: str
    contig_len_bp: int
    n_genes_on_contig: int
    anchor_protein_id: str
    anchor_start: int
    anchor_end: int
    anchor_strand: int
    rank: int  # 0 = anchor, -1 = one gene upstream, etc.
    protein_id: str
    locus_tag: str
    product: str
    start: int
    end: int
    strand: int
    length_aa: int
    intergenic_gap_bp: Optional[int]
    same_strand_as_anchor: bool
    hmm_hit: bool
    hmm_evalue: Optional[float]
    hmm_bitscore: Optional[float]
    hmm_coverage: Optional[float]
    dist_to_contig_start_bp: int
    dist_to_contig_end_bp: int
    window_truncated: bool
    fusion_candidate: bool


# --------------------------------------------------------------------------
# HTTP helpers
# --------------------------------------------------------------------------


class Fetcher:
    """Polite HTTP client with on-disk caching and retry/backoff."""

    def __init__(self, cache_dir: Path, sleep: float = 0.25, retries: int = 4):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.sleep = sleep
        self.retries = retries
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": USER_AGENT})

    def _cache_path(self, key: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", key)[:180]
        return self.cache_dir / safe

    def get_text(self, url: str, params: Optional[dict] = None,
                 cache_key: Optional[str] = None) -> Optional[str]:
        key = cache_key or (url + "?" + json.dumps(params or {}, sort_keys=True))
        cached = self._cache_path(key)
        if cached.exists():
            return cached.read_text()

        delay = 1.0
        for attempt in range(self.retries):
            try:
                r = self.session.get(url, params=params, timeout=90)
            except requests.RequestException as exc:
                sys.stderr.write(f"  [warn] {exc} -- retrying\n")
                time.sleep(delay)
                delay *= 2
                continue

            if r.status_code == 200:
                cached.write_text(r.text)
                time.sleep(self.sleep)
                return r.text
            if r.status_code == 204:
                # InterPro returns 204 for an empty result set
                cached.write_text("")
                return ""
            if r.status_code == 404:
                sys.stderr.write(f"  [warn] 404 {r.url}\n")
                return None
            if r.status_code in (429, 500, 502, 503, 504):
                sys.stderr.write(
                    f"  [warn] HTTP {r.status_code} -- backing off {delay:.0f}s\n"
                )
                time.sleep(delay)
                delay *= 2
                continue
            sys.stderr.write(f"  [warn] HTTP {r.status_code} for {r.url}\n")
            return None

        sys.stderr.write(f"  [error] gave up on {url}\n")
        return None

    def get_json(self, url: str, params: Optional[dict] = None,
                 cache_key: Optional[str] = None) -> Optional[dict]:
        text = self.get_text(url, params=params, cache_key=cache_key)
        if not text:
            return None
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            sys.stderr.write(f"  [warn] non-JSON response from {url}\n")
            return None


# --------------------------------------------------------------------------
# Stage 1: Pfam list -> domain architectures
# --------------------------------------------------------------------------


def read_pfam_list(path: Path) -> List[str]:
    """One Pfam accession per line; blank lines and # comments ignored."""
    pfams = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        acc = line.upper()
        if not re.fullmatch(r"PF\d{5}", acc):
            sys.stderr.write(f"  [warn] skipping malformed Pfam accession: {line}\n")
            continue
        pfams.append(acc)
    if not pfams:
        sys.exit("No valid Pfam accessions found in the input file.")
    return pfams


def find_architectures(fetcher: Fetcher, pfams: Sequence[str], ordered: bool,
                       exact: bool, ignore: Sequence[str]) -> List[dict]:
    """Query the InterPro IDA endpoint for architectures containing the Pfams."""
    params: Dict[str, object] = {"ida_search": ",".join(pfams)}
    if ordered:
        params["ordered"] = ""
    if exact:
        params["exact"] = ""
    if ignore:
        params["ida_ignore"] = ",".join(ignore)

    data = fetcher.get_json(f"{INTERPRO_API}/entry", params=params)
    if not data:
        return []

    architectures = []
    for res in data.get("results", []):
        # Field naming has shifted across InterPro releases; accept either.
        ida_id = res.get("ida_id") or res.get("ida_accession") or ""
        ida_string = res.get("ida") or ""
        count = res.get("counts") or res.get("unique_proteins") or 0
        if not ida_id:
            continue
        n_domains = len([d for d in ida_string.split("-") if d]) if ida_string else 0
        architectures.append(
            {
                "ida_id": ida_id,
                "ida": ida_string,
                "count": count,
                "n_domains": n_domains,
            }
        )
    return architectures


# --------------------------------------------------------------------------
# Stage 2: architecture -> UniProt accessions
# --------------------------------------------------------------------------


def proteins_for_ida(fetcher: Fetcher, ida_id: str, reviewed_only: bool,
                     max_proteins: int) -> List[str]:
    """Page through the InterPro protein endpoint for one architecture."""
    source = "reviewed" if reviewed_only else "uniprot"
    url = f"{INTERPRO_API}/protein/{source}"
    params: Optional[dict] = {"ida": ida_id, "page_size": 200}
    accessions: List[str] = []
    page = 0

    while url and len(accessions) < max_proteins:
        data = fetcher.get_json(url, params=params,
                                cache_key=f"ida_prot_{ida_id}_{page}")
        params = None  # the `next` URL already carries the query string
        if not data:
            break
        for res in data.get("results", []):
            meta = res.get("metadata", {})
            acc = meta.get("accession")
            if acc:
                accessions.append(acc.upper())
        url = data.get("next")
        page += 1

    return accessions[:max_proteins]


# --------------------------------------------------------------------------
# Stage 3: UniProt -> contig accession + protein_id
# --------------------------------------------------------------------------


def uniprot_records(fetcher: Fetcher, accessions: Sequence[str],
                    batch_size: int = 100) -> Dict[str, dict]:
    """Batch-fetch UniProt entries and pull out ENA cross-references."""
    out: Dict[str, dict] = {}

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        query = " OR ".join(f"accession:{a}" for a in batch)
        data = fetcher.get_json(
            f"{UNIPROT_API}/uniprotkb/stream",
            params={
                "query": f"({query})",
                "format": "json",
                "fields": ("accession,xref_embl,xref_proteomes,"
                           "organism_name,organism_id,length"),
            },
            cache_key=f"uniprot_batch_{i}_{batch[0]}_{batch[-1]}",
        )
        if not data:
            continue

        for entry in data.get("results", []):
            acc = entry.get("primaryAccession")
            if not acc:
                continue
            organism = entry.get("organism", {}) or {}
            rec = {
                "accession": acc,
                "organism_name": organism.get("scientificName", ""),
                "tax_id": str(organism.get("taxonId", "")),
                "length": (entry.get("sequence", {}) or {}).get("length", 0),
                "embl": [],       # list of (contig_acc, protein_id)
                "proteomes": [],
            }
            for xref in entry.get("uniProtKBCrossReferences", []) or []:
                db = xref.get("database")
                props = {p.get("key"): p.get("value")
                         for p in (xref.get("properties") or [])}
                if db == "EMBL":
                    contig = (xref.get("id") or "").strip().rstrip(";")
                    protein_id = (props.get("ProteinId") or "").strip()
                    if contig:
                        rec["embl"].append((contig, protein_id))
                elif db == "Proteomes":
                    pid = (xref.get("id") or "").strip()
                    if pid:
                        rec["proteomes"].append(pid)
            out[acc] = rec

    return out


def proteome_to_assembly(fetcher: Fetcher, proteome_id: str) -> str:
    """Resolve a UniProt proteome ID to a GCA/GCF assembly accession."""
    if not proteome_id:
        return ""
    data = fetcher.get_json(f"{UNIPROT_API}/proteomes/{proteome_id}.json",
                            cache_key=f"proteome_{proteome_id}")
    if not data:
        return ""
    return data.get("genomeAssembly", {}).get("assemblyId", "") or ""


# --------------------------------------------------------------------------
# Stage 4/5: ENA contig -> genes -> neighbourhood
# --------------------------------------------------------------------------


def fetch_contig(fetcher: Fetcher, contig_acc: str, contig_dir: Path) -> Optional[Path]:
    """Download an EMBL flatfile for one contig, skipping WGS master records."""
    if WGS_MASTER_RE.match(contig_acc):
        sys.stderr.write(
            f"  [warn] {contig_acc} looks like a WGS master record "
            "(no CDS features) -- skipping\n"
        )
        return None

    path = contig_dir / f"{contig_acc}.embl"
    if path.exists() and path.stat().st_size > 0:
        return path

    text = fetcher.get_text(f"{ENA_API}/embl/{contig_acc}",
                            cache_key=f"ena_embl_{contig_acc}")
    if not text or "FT   CDS" not in text:
        sys.stderr.write(f"  [warn] no CDS features retrieved for {contig_acc}\n")
        return None

    path.write_text(text)
    return path


def parse_contig(path: Path, contig_acc: str) -> tuple:
    """Parse an EMBL flatfile into an ordered list of Gene objects."""
    genes: List[Gene] = []
    contig_len = 0

    # The declared length on the ID line is authoritative; len(record.seq) can
    # disagree if the flatfile is truncated or the sequence block is absent.
    with path.open() as fh:
        for line in fh:
            if line.startswith("ID   "):
                m = re.search(r"(\d+)\s+BP\.", line)
                if m:
                    contig_len = int(m.group(1))
                break

    for record in SeqIO.parse(str(path), "embl"):
        contig_len = max(contig_len, len(record.seq))
        for feat in record.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            translation = (q.get("translation") or [""])[0]
            if not translation:
                continue  # pseudogenes carry no translation
            genes.append(
                Gene(
                    contig_acc=contig_acc,
                    protein_id=(q.get("protein_id") or [""])[0].split(".")[0],
                    locus_tag=(q.get("locus_tag") or [""])[0],
                    product=(q.get("product") or [""])[0],
                    start=int(feat.location.start) + 1,  # 1-based, inclusive
                    end=int(feat.location.end),
                    strand=feat.location.strand or 0,
                    translation=translation,
                )
            )

    genes.sort(key=lambda g: (g.start, g.end))
    for i, g in enumerate(genes):
        g.index = i
    return genes, contig_len


def locate_anchor(genes: Sequence[Gene], protein_id: str,
                  uniprot_len: int) -> Optional[int]:
    """Find the anchor gene by protein_id, falling back to sequence length."""
    if protein_id:
        stem = protein_id.split(".")[0]
        for g in genes:
            if g.protein_id == stem:
                return g.index

    # Fallback: unique length match. Ambiguous matches are rejected rather
    # than guessed, so a wrong anchor never silently enters the report.
    if uniprot_len:
        hits = [g.index for g in genes if abs(g.length_aa - uniprot_len) <= 1]
        if len(hits) == 1:
            return hits[0]

    return None


# --------------------------------------------------------------------------
# Stage 6: HMMER
# --------------------------------------------------------------------------


def run_hmmsearch(hmm: Path, faa: Path, out_dir: Path, evalue: float,
                  cpus: int) -> Dict[str, dict]:
    """Run hmmsearch and parse --domtblout into best-hit-per-target."""
    domtbl = out_dir / "neighbours.domtbl"
    cmd = [
        "hmmsearch",
        "--domtblout", str(domtbl),
        "-E", str(evalue),
        "--cpu", str(cpus),
        "-o", str(out_dir / "hmmsearch.log"),
        str(hmm), str(faa),
    ]
    sys.stderr.write(f"  running: {' '.join(cmd)}\n")
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        sys.exit("hmmsearch not found on $PATH. Install HMMER3 and retry.")
    except subprocess.CalledProcessError as exc:
        sys.exit(f"hmmsearch failed with exit code {exc.returncode}.")

    hits: Dict[str, dict] = {}
    for line in domtbl.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split(None, 22)
        if len(f) < 22:
            continue
        target, tlen = f[0], int(f[2])
        qlen = int(f[5])
        full_evalue, full_score = float(f[6]), float(f[7])
        hmm_from, hmm_to = int(f[15]), int(f[16])
        coverage = (hmm_to - hmm_from + 1) / qlen if qlen else 0.0

        prev = hits.get(target)
        if prev is None or full_score > prev["bitscore"]:
            hits[target] = {
                "evalue": full_evalue,
                "bitscore": full_score,
                "coverage": round(coverage, 3),
                "target_len": tlen,
            }
    return hits


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def sketch_neighbourhood(rows: Sequence[NeighbourRow]) -> str:
    """Render one locus as a single-line arrow diagram."""
    parts = []
    for r in sorted(rows, key=lambda x: x.rank):
        label = r.protein_id or r.locus_tag or "?"
        if r.hmm_hit:
            label = f"*{label}*"
        if r.rank == 0:
            label = f"[{label}]"
        parts.append(f"{label}->" if r.strand >= 0 else f"<-{label}")
    return "  ".join(parts)


def write_reports(rows: List[NeighbourRow], out_dir: Path,
                  architectures: List[dict], stats: dict) -> None:
    # 1. Full per-gene table
    genes_tsv = out_dir / "neighbourhood_genes.tsv"
    with genes_tsv.open("w", newline="") as fh:
        if rows:
            writer = csv.DictWriter(fh, fieldnames=list(asdict(rows[0]).keys()),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for r in rows:
                writer.writerow(asdict(r))

    # 2. Per-locus summary, one row per anchor
    loci: Dict[tuple, List[NeighbourRow]] = {}
    for r in rows:
        loci.setdefault((r.contig_acc, r.anchor_protein_id, r.uniprot_acc), []).append(r)

    summary_tsv = out_dir / "loci_summary.tsv"
    with summary_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "uniprot_acc", "organism_name", "tax_id", "assembly", "contig_acc",
            "anchor_protein_id", "anchor_start", "anchor_end", "anchor_strand",
            "n_neighbours", "n_hmm_hits", "best_hmm_bitscore", "best_hmm_evalue",
            "hmm_hit_ranks", "fusion_candidate", "window_truncated",
            "anchor_dist_to_contig_start_bp", "anchor_dist_to_contig_end_bp",
            "ida_string", "sketch",
        ])
        for (contig, anchor_pid, acc), group in sorted(loci.items()):
            hits = [g for g in group if g.hmm_hit]
            anchor = next((g for g in group if g.rank == 0), group[0])
            best = max(hits, key=lambda g: g.hmm_bitscore or 0) if hits else None
            writer.writerow([
                acc, anchor.organism_name, anchor.tax_id, anchor.assembly, contig,
                anchor_pid, anchor.anchor_start, anchor.anchor_end,
                anchor.anchor_strand,
                len(group) - 1, len(hits),
                best.hmm_bitscore if best else "",
                best.hmm_evalue if best else "",
                ",".join(str(g.rank) for g in sorted(hits, key=lambda x: x.rank)),
                anchor.fusion_candidate,
                any(g.window_truncated for g in group),
                anchor.dist_to_contig_start_bp, anchor.dist_to_contig_end_bp,
                anchor.ida_string,
                sketch_neighbourhood(group),
            ])

    # 3. Human-readable report
    report = out_dir / "report.txt"
    with report.open("w") as fh:
        fh.write("ida2synteny report\n")
        fh.write("=" * 70 + "\n\n")
        fh.write(f"Query Pfams        : {stats['pfams']}\n")
        fh.write(f"Query HMM          : {stats['hmm']}\n")
        fh.write(f"Architectures found: {len(architectures)}\n")
        fh.write(f"Anchor proteins    : {stats['n_anchors']}\n")
        fh.write(f"Contigs retrieved  : {stats['n_contigs']}\n")
        fh.write(f"Loci reconstructed : {len(loci)}\n")
        fh.write(f"Loci with HMM hit  : "
                 f"{sum(1 for g in loci.values() if any(x.hmm_hit for x in g))}\n")
        fh.write(f"Fusion candidates  : {stats['n_fusion']}\n\n")

        fh.write("Architectures\n" + "-" * 70 + "\n")
        for a in architectures:
            fh.write(f"  {a['ida_id'][:12]}  n={a['count']:<7} "
                     f"domains={a['n_domains']}  {a['ida']}\n")
        fh.write("\n")

        fh.write("Loci\n" + "-" * 70 + "\n")
        fh.write("Legend: [anchor]  *HMM hit*  > forward  < reverse\n\n")
        for (contig, anchor_pid, acc), group in sorted(loci.items()):
            anchor = next((g for g in group if g.rank == 0), group[0])
            flags = []
            if anchor.fusion_candidate:
                flags.append("FUSION-CANDIDATE")
            if any(g.window_truncated for g in group):
                flags.append("WINDOW-TRUNCATED")
            if anchor.dist_to_contig_start_bp < 5000 or \
               anchor.dist_to_contig_end_bp < 5000:
                flags.append("NEAR-CONTIG-EDGE")
            hits = [g for g in group if g.hmm_hit]

            fh.write(f"{acc}  {anchor.organism_name} (taxid {anchor.tax_id})\n")
            fh.write(f"  contig {contig} ({anchor.contig_len_bp} bp, "
                     f"{anchor.n_genes_on_contig} CDS)  assembly "
                     f"{anchor.assembly or 'n/a'}\n")
            fh.write(f"  anchor {anchor_pid} at {anchor.anchor_start}"
                     f"..{anchor.anchor_end} ({'+' if anchor.anchor_strand >= 0 else '-'})\n")
            fh.write(f"  {sketch_neighbourhood(group)}\n")
            if hits:
                for h in sorted(hits, key=lambda x: x.rank):
                    fh.write(f"    HMM hit rank {h.rank:+d}: {h.protein_id} "
                             f"bitscore={h.bitscore_str()} "
                             f"E={h.hmm_evalue:.2g} cov={h.hmm_coverage} "
                             f"| {h.product}\n")
            else:
                fh.write("    no HMM hits in window\n")
            if flags:
                fh.write(f"    flags: {', '.join(flags)}\n")
            fh.write("\n")

    sys.stderr.write(f"\nWrote:\n  {genes_tsv}\n  {summary_tsv}\n  {report}\n")


# small helper attached post-hoc to keep the dataclass declaration clean
def _bitscore_str(self: NeighbourRow) -> str:
    return f"{self.hmm_bitscore:.1f}" if self.hmm_bitscore is not None else "-"


NeighbourRow.bitscore_str = _bitscore_str  # type: ignore[attr-defined]


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(
        description="Pfam architecture -> isolate genome synteny -> HMM scan.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--pfams", required=True, type=Path,
                   help="Text file of Pfam accessions, one per line.")
    p.add_argument("--hmm", required=True, type=Path,
                   help="Query HMM file (e.g. an MGnifam) to scan neighbours with.")
    p.add_argument("-o", "--outdir", required=True, type=Path,
                   help="Output directory.")
    p.add_argument("--window", type=int, default=5,
                   help="Genes to report either side of the anchor (default: 5).")
    p.add_argument("--ordered", action="store_true",
                   help="Require the Pfams in the given order.")
    p.add_argument("--exact", action="store_true",
                   help="Only architectures containing exactly these Pfams.")
    p.add_argument("--ida-ignore", default="",
                   help="Comma-separated Pfams to exclude from architectures.")
    p.add_argument("--reviewed-only", action="store_true",
                   help="Restrict to Swiss-Prot rather than all of UniProtKB.")
    p.add_argument("--max-proteins", type=int, default=500,
                   help="Cap on anchor proteins per architecture (default: 500).")
    p.add_argument("--evalue", type=float, default=1e-5,
                   help="hmmsearch full-sequence E-value threshold (default: 1e-5).")
    p.add_argument("--cpus", type=int, default=4, help="Threads for hmmsearch.")
    p.add_argument("--cache", type=Path, default=None,
                   help="Cache directory (default: <outdir>/cache).")
    args = p.parse_args()

    if not args.hmm.exists():
        sys.exit(f"HMM file not found: {args.hmm}")

    out_dir = args.outdir
    out_dir.mkdir(parents=True, exist_ok=True)
    contig_dir = out_dir / "contigs"
    contig_dir.mkdir(exist_ok=True)
    fetcher = Fetcher(args.cache or (out_dir / "cache"))

    pfams = read_pfam_list(args.pfams)
    ignore = [x.strip().upper() for x in args.ida_ignore.split(",") if x.strip()]

    # -- Stage 1 ------------------------------------------------------------
    sys.stderr.write(f"[1/6] InterPro IDA search: {','.join(pfams)}\n")
    architectures = find_architectures(fetcher, pfams, args.ordered,
                                       args.exact, ignore)
    if not architectures:
        sys.exit("No domain architectures matched. Check the Pfam accessions, "
                 "or relax --ordered/--exact.")
    total = sum(a["count"] for a in architectures)
    sys.stderr.write(f"      {len(architectures)} architectures, "
                     f"{total} proteins total\n")

    # -- Stage 2 ------------------------------------------------------------
    sys.stderr.write("[2/6] Resolving architectures to UniProt accessions\n")
    anchors: List[Anchor] = []
    seen = set()
    for arch in architectures:
        accs = proteins_for_ida(fetcher, arch["ida_id"], args.reviewed_only,
                                args.max_proteins)
        for acc in accs:
            if acc in seen:
                continue
            seen.add(acc)
            anchors.append(Anchor(
                uniprot_acc=acc,
                ida_id=arch["ida_id"],
                ida_string=arch["ida"],
                n_domains=arch["n_domains"],
                # An architecture longer than the query is a candidate fusion:
                # domains that normally sit on separate genes merged into one.
                is_fusion_candidate=arch["n_domains"] > len(pfams),
            ))
    sys.stderr.write(f"      {len(anchors)} unique anchor proteins\n")
    if not anchors:
        sys.exit("No proteins retrieved for these architectures.")

    # -- Stage 3 ------------------------------------------------------------
    sys.stderr.write("[3/6] Fetching UniProt cross-references\n")
    records = uniprot_records(fetcher, [a.uniprot_acc for a in anchors])
    assembly_cache: Dict[str, str] = {}
    resolved: List[Anchor] = []
    for a in anchors:
        rec = records.get(a.uniprot_acc)
        if not rec:
            continue
        a.organism_name = rec["organism_name"]
        a.tax_id = rec["tax_id"]
        a.proteome = rec["proteomes"][0] if rec["proteomes"] else ""
        if not rec["embl"]:
            sys.stderr.write(f"  [warn] {a.uniprot_acc} has no ENA cross-reference "
                             "-- no genomic context available\n")
            continue
        # One UniProt entry can map to many identical records across assemblies.
        # Take the first; add --all-xrefs here if you want every instance.
        a.contig_acc, a.protein_id = rec["embl"][0]
        resolved.append(a)
    sys.stderr.write(f"      {len(resolved)}/{len(anchors)} anchors have a contig\n")

    # -- Stage 4/5 ----------------------------------------------------------
    sys.stderr.write("[4/6] Retrieving contigs and extracting neighbourhoods\n")
    contig_cache: Dict[str, tuple] = {}
    pending: List[tuple] = []   # (anchor, genes, window, contig_len)
    faa_records: List[tuple] = []

    for a in resolved:
        if a.contig_acc not in contig_cache:
            path = fetch_contig(fetcher, a.contig_acc, contig_dir)
            if path is None:
                contig_cache[a.contig_acc] = ([], 0)
            else:
                contig_cache[a.contig_acc] = parse_contig(path, a.contig_acc)
        genes, contig_len = contig_cache[a.contig_acc]
        if not genes:
            continue

        idx = locate_anchor(genes, a.protein_id,
                            records.get(a.uniprot_acc, {}).get("length", 0))
        if idx is None:
            sys.stderr.write(f"  [warn] could not locate {a.protein_id or a.uniprot_acc} "
                             f"on {a.contig_acc}\n")
            continue

        lo, hi = idx - args.window, idx + args.window
        truncated = lo < 0 or hi >= len(genes)
        window = genes[max(0, lo):min(len(genes), hi + 1)]
        pending.append((a, window, genes[idx], contig_len, len(genes), truncated))

        for g in window:
            key = f"{a.contig_acc}|{g.index}|{g.protein_id or g.locus_tag or 'na'}"
            faa_records.append((key, g.translation))

    n_contigs = sum(1 for v in contig_cache.values() if v[0])
    sys.stderr.write(f"      {n_contigs} contigs parsed, {len(pending)} loci, "
                     f"{len(faa_records)} neighbour proteins\n")
    if not pending:
        sys.exit("No loci reconstructed -- nothing to scan.")

    # -- Stage 6 ------------------------------------------------------------
    sys.stderr.write("[5/6] Running hmmsearch against neighbour proteins\n")
    faa = out_dir / "neighbours.faa"
    with faa.open("w") as fh:
        written = set()
        for key, seq in faa_records:
            if key in written:
                continue
            written.add(key)
            fh.write(f">{key}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    hmm_hits = run_hmmsearch(args.hmm, faa, out_dir, args.evalue, args.cpus)
    sys.stderr.write(f"      {len(hmm_hits)} proteins hit the HMM\n")

    # -- Assemble rows ------------------------------------------------------
    sys.stderr.write("[6/6] Building report\n")
    rows: List[NeighbourRow] = []
    for a, window, anchor_gene, contig_len, n_genes, truncated in pending:
        if a.proteome and a.proteome not in assembly_cache:
            assembly_cache[a.proteome] = proteome_to_assembly(fetcher, a.proteome)
        assembly = assembly_cache.get(a.proteome, "")

        prev_end = None
        for g in window:
            key = f"{a.contig_acc}|{g.index}|{g.protein_id or g.locus_tag or 'na'}"
            hit = hmm_hits.get(key)
            gap = (g.start - prev_end - 1) if prev_end is not None else None
            prev_end = g.end
            rows.append(NeighbourRow(
                ida_id=a.ida_id,
                ida_string=a.ida_string,
                uniprot_acc=a.uniprot_acc,
                organism_name=a.organism_name,
                tax_id=a.tax_id,
                proteome=a.proteome,
                assembly=assembly,
                contig_acc=a.contig_acc,
                contig_len_bp=contig_len,
                n_genes_on_contig=n_genes,
                anchor_protein_id=anchor_gene.protein_id,
                anchor_start=anchor_gene.start,
                anchor_end=anchor_gene.end,
                anchor_strand=anchor_gene.strand,
                rank=g.index - anchor_gene.index,
                protein_id=g.protein_id,
                locus_tag=g.locus_tag,
                product=g.product,
                start=g.start,
                end=g.end,
                strand=g.strand,
                length_aa=g.length_aa,
                intergenic_gap_bp=gap,
                same_strand_as_anchor=(g.strand == anchor_gene.strand),
                hmm_hit=hit is not None,
                hmm_evalue=hit["evalue"] if hit else None,
                hmm_bitscore=hit["bitscore"] if hit else None,
                hmm_coverage=hit["coverage"] if hit else None,
                dist_to_contig_start_bp=anchor_gene.start - 1,
                dist_to_contig_end_bp=max(0, contig_len - anchor_gene.end),
                window_truncated=truncated,
                fusion_candidate=a.is_fusion_candidate,
            ))

    stats = {
        "pfams": ",".join(pfams),
        "hmm": str(args.hmm),
        "n_anchors": len(resolved),
        "n_contigs": n_contigs,
        "n_fusion": sum(1 for a in resolved if a.is_fusion_candidate),
    }
    write_reports(rows, out_dir, architectures, stats)


if __name__ == "__main__":
    main()
