#!/usr/bin/env python3

"""Build per-family domain architecture JSONs from the MGnify re-annotation CSV.

See SPEC.md for the frozen behaviour. Replaces the pipeline's parse_domains.py, which read
per-family TSVs from the old MGnify protein database and assumed one MGnifam hit per sequence.
"""

import argparse
import csv
import gzip
import io
import json
import logging
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict, namedtuple
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

# Field order is the sort order: ali_from, then ali_to (shortest first on a tie), then id.
Chip = namedtuple("Chip", "start end id name link")

# Families absent from the clan file are reported once each rather than per sequence.
_UNKNOWN_FAMILIES = set()


def load_clan_membership(clan_membership_file):
    """Return (family_id -> clan_id, clan_id -> representative family_id)."""
    family_to_clan = {}
    clan_to_rep = {}

    with open(clan_membership_file, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            clan = row["Cluster Id"]
            clan_to_rep[clan] = row["Family Rep Id"]
            for family_id in row["Family Ids"].split(";"):
                family_id = family_id.strip()
                if family_id:
                    family_to_clan[family_id] = clan

    log.info(f"Loaded {len(clan_to_rep)} clans covering {len(family_to_clan)} families")

    return family_to_clan, clan_to_rep


def load_pfam_mapping(pfam_mapping_file):
    """Return accession -> human readable name.

    Accepts the two-column `accession<TAB>name` form and the wider
    `pfam_id<TAB>name<TAB>clan_id` form; anything past the name is ignored, as is a header row.
    """
    pfam_mapping = {}
    skipped = 0

    with open(pfam_mapping_file) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0].startswith("PF"):
                skipped += 1
                continue
            pfam_mapping[fields[0]] = fields[1]

    log.info(f"Loaded {len(pfam_mapping)} pfam entries")
    if skipped:
        log.info(f"Skipped {skipped} non-accession lines in {pfam_mapping_file}")

    return pfam_mapping


def overlaps(a, b, overlap_fraction):
    """True when the two (start, end) spans share more than overlap_fraction of the shorter one."""
    shared = min(a[1], b[1]) - max(a[0], b[0]) + 1
    if shared <= 0:
        return False

    shorter = min(a[1] - a[0] + 1, b[1] - b[0] + 1)

    return shared > overlap_fraction * shorter


def cluster_hits(hits, family_to_clan, overlap_fraction):
    """Single-linkage cluster (family_id, start, end) hits, merging only within a clan.

    A family missing from the clan file gets a private pseudo-clan, so it can only ever merge
    with itself. Returns groups ordered by span.
    """
    by_clan = {}
    for hit in hits:
        clan = family_to_clan.get(hit[0])
        if clan is None:
            clan = f"__unknown__{hit[0]}"
            if hit[0] not in _UNKNOWN_FAMILIES:
                _UNKNOWN_FAMILIES.add(hit[0])
                log.warning(f"Family {hit[0]} is absent from the clan file, never merging it")
        by_clan.setdefault(clan, []).append(hit)

    groups = []
    for clan_hits in by_clan.values():
        # ponytail: O(n^2) single linkage. n is the hit count on one sequence (tens at worst);
        # switch to union-find if a pathological sequence ever makes this the bottleneck.
        clan_groups = []
        for hit in clan_hits:
            merged = [hit]
            kept = []
            for group in clan_groups:
                if any(overlaps(hit[1:], other[1:], overlap_fraction) for other in group):
                    merged.extend(group)
                else:
                    kept.append(group)
            clan_groups = kept + [merged]
        groups.extend(clan_groups)

    return sorted(groups, key=lambda group: (min(h[1] for h in group), max(h[2] for h in group)))


def construct_name(family_id):
    return "MGYF" + f"{int(family_id):010d}"


def details_link(family_id, base_url):
    return f"{base_url.rstrip('/')}/{construct_name(family_id)}"


def mgnifam_chip(group, family_to_clan, clan_to_rep, base_url):
    """One merged group of MGnifam hits becomes one chip, clan-labelled only if it spans >1 family."""
    start = min(hit[1] for hit in group)
    end = max(hit[2] for hit in group)
    families = {hit[0] for hit in group}

    if len(families) > 1:
        clan = family_to_clan[next(iter(families))]
        return Chip(start, end, clan, f"MGnifam clan {clan.split('_', 1)[1]}",
                    details_link(clan_to_rep[clan], base_url))

    family_id = families.pop()

    return Chip(start, end, family_id, f"MGnifam{family_id}", details_link(family_id, base_url))


def pfam_chip(hit, pfam_mapping):
    """One "p" entry becomes one chip, positioned by its alignment (not HMM) coordinates."""
    accession = hit[0]

    return Chip(hit[5], hit[6], accession, pfam_mapping.get(accession, accession), None)


def build_chips(metadata, family_to_clan, clan_to_rep, pfam_mapping, base_url, overlap_fraction):
    """Every chip on one sequence, ordered by (ali_from, ali_to, id)."""
    chips = [pfam_chip(hit, pfam_mapping) for hit in metadata.get("p") or []]

    hits = [(hit[0], hit[3], hit[4]) for hit in metadata.get("m") or []]
    for group in cluster_hits(hits, family_to_clan, overlap_fraction):
        chips.append(mgnifam_chip(group, family_to_clan, clan_to_rep, base_url))

    return sorted(chips)


def architecture_key(chips):
    return "\t".join(chip.id for chip in chips)


def string_to_hex_color(s):
    hash_val = 0

    for char in s:
        hash_val = ord(char) + ((hash_val << 5) - hash_val)

    color = '#'

    for i in range(3):
        value = (hash_val >> (i * 8)) & 0xFF
        color += ('00' + format(value, 'x'))[-2:]

    return color


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')

    return (int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16))


def calculate_luminosity(rgb):
    def linearize(color):
        c = color / 255.0
        if c <= 0.03928:
            return c / 12.92
        return ((c + 0.055) / 1.055) ** 2.4

    r, g, b = rgb

    return 0.2126 * linearize(r) + 0.7152 * linearize(g) + 0.0722 * linearize(b)


def decide_font_color(hex_color):
    return 'black' if calculate_luminosity(hex_to_rgb(hex_color)) > 0.2 else 'white'


def raise_csv_field_limit():
    """Protein metadata fields can exceed Python's default 131072-byte csv limit."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = int(limit / 10)


def _metadata_column(proteins_file):
    """Index of the metadata column, read from the header alone."""
    opener = gzip.open if str(proteins_file).endswith(".gz") else open
    with opener(proteins_file, "rt", newline="") as handle:
        header = next(csv.reader(handle))

    return header.index("metadata")


def iter_metadata(proteins_file, use_prefilter=True):
    """Yield the parsed metadata of every row carrying an "m" annotation.

    With use_prefilter, decompression and the discard of "m"-less rows are pushed into
    `zcat | grep`, so only survivors reach json.loads. The "m" check is repeated here, so both
    paths yield exactly the same rows.
    """
    raise_csv_field_limit()
    column = _metadata_column(proteins_file)
    prefiltered = use_prefilter and str(proteins_file).endswith(".gz")

    if prefiltered:
        # The raw CSV doubles the quotes inside the metadata field, hence '""m"":'.
        zcat = subprocess.Popen(["zcat", str(proteins_file)], stdout=subprocess.PIPE)
        grep = subprocess.Popen(["grep", "-F", '""m"":'], stdin=zcat.stdout,
                                stdout=subprocess.PIPE, env={**os.environ, "LC_ALL": "C"})
        zcat.stdout.close()
        handle = io.TextIOWrapper(grep.stdout, newline="")
    else:
        opener = gzip.open if str(proteins_file).endswith(".gz") else open
        handle = opener(proteins_file, "rt", newline="")

    try:
        reader = csv.reader(handle)
        if not prefiltered:
            next(reader, None)  # grep already dropped the header on the prefiltered path
        for row in reader:
            if len(row) <= column:
                continue
            metadata = json.loads(row[column] or "{}")
            if metadata.get("m"):
                yield metadata
    finally:
        handle.close()
        if prefiltered:
            zcat.wait()
            grep.wait()
            if grep.returncode not in (0, 1):  # 1 is grep's "no lines matched"
                raise RuntimeError(f"prefilter failed with exit code {grep.returncode}")


def count_architectures(metadata_rows, family_to_clan, clan_to_rep, pfam_mapping, base_url,
                        overlap_fraction, log_every=1000000):
    """Tally architecture keys per family. Every family hit on a row is credited exactly once."""
    # ponytail: the whole tally (~35K families) lives in RAM for the pass. If it ever stops
    # fitting, spill `family\tkey` lines to disk and finish with `sort | uniq -c`.
    counters = defaultdict(Counter)
    kept = 0

    for kept, metadata in enumerate(metadata_rows, start=1):
        chips = build_chips(metadata, family_to_clan, clan_to_rep, pfam_mapping, base_url,
                            overlap_fraction)
        key = sys.intern(architecture_key(chips))
        for family_id in {hit[0] for hit in metadata["m"]}:
            counters[family_id][key] += 1

        if kept % log_every == 0:
            log.info(f"Counted {kept} annotated sequences across {len(counters)} families")

    log.info(f"Counted {kept} annotated sequences across {len(counters)} families")

    return counters


def resolve_chip(chip_id, clan_to_rep, pfam_mapping, base_url):
    """Turn a chip id back into the viewer's domain object."""
    if chip_id in clan_to_rep:
        name = f"MGnifam clan {chip_id.split('_', 1)[1]}"
        link = details_link(clan_to_rep[chip_id], base_url)
    elif chip_id.isdigit():
        name = f"MGnifam{chip_id}"
        link = details_link(chip_id, base_url)
    else:
        name = pfam_mapping.get(chip_id, chip_id)
        link = f"https://www.ebi.ac.uk/interpro/entry/pfam/{chip_id}"

    color = string_to_hex_color(name)

    return {"id": chip_id, "color": color, "link": link, "name": name,
            "font_color": decide_font_color(color)}


def architecture_json(counter, top, clan_to_rep, pfam_mapping, base_url):
    containers = []

    for key, count in Counter(counter).most_common(top):
        domains = [resolve_chip(chip_id, clan_to_rep, pfam_mapping, base_url)
                   for chip_id in key.split("\t") if chip_id]
        containers.append({"architecture_text": str(count), "domains": domains})

    return {"architecture_containers": containers}


def _family_sort_key(family_id):
    return (0, int(family_id), "") if family_id.isdigit() else (1, 0, family_id)


def write_outputs(counters, expected_families, output_dir, top, clan_to_rep, pfam_mapping,
                  base_url):
    """Write one JSON per family and return the expected families that got no hits.

    Families with no hits still get a valid empty file, so updating the database cannot fail
    on a missing input.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for family_id, counter in counters.items():
        payload = architecture_json(counter, top, clan_to_rep, pfam_mapping, base_url)
        (output_dir / f"{family_id}.json").write_text(json.dumps(payload, indent=4))

    missing = sorted(set(expected_families) - set(counters), key=_family_sort_key)
    for family_id in missing:
        (output_dir / f"{family_id}.json").write_text(
            json.dumps({"architecture_containers": []}, indent=4))

    (output_dir / "missing_families.txt").write_text(
        "".join(f"{family_id}\n" for family_id in missing))

    log.info(f"Wrote {len(counters)} families with hits to {output_dir}")
    if missing:
        log.warning(f"{len(missing)} expected families had no annotated sequence; wrote empty "
                    f"architectures and listed them in {output_dir / 'missing_families.txt'}")

    unexpected = sorted(set(counters) - set(expected_families), key=_family_sort_key)
    if unexpected:
        log.warning(f"{len(unexpected)} families were annotated but absent from the clan file, "
                    f"e.g. {unexpected[:5]}")

    return missing


def main():
    parser = argparse.ArgumentParser(
        description="Build per-family domain architecture JSONs from the MGnify re-annotation CSV.")
    parser.add_argument("--proteins", required=True,
                        help="Re-annotated proteins CSV(.gz) with mgyp and metadata columns")
    parser.add_argument("--clan-membership", required=True,
                        help="clan_membership.csv with Cluster Id, Family Rep Id and Family Ids")
    parser.add_argument("--pfam-mapping", required=True,
                        help="TSV mapping Pfam accessions to names")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for the per-family JSONs")
    parser.add_argument("--overlap-fraction", type=float, default=0.5,
                        help="Same-clan hits merge when they share more than this fraction of the "
                             "shorter hit")
    parser.add_argument("--top", type=int, default=15,
                        help="Number of architectures kept per family")
    parser.add_argument("--base-url", default="http://mgnifams-demo.mgnify.org/details/",
                        help="Base URL for MGnifam detail pages")
    parser.add_argument("--log-every", type=int, default=1000000,
                        help="Log progress every N annotated sequences")
    parser.add_argument("--no-prefilter", action="store_true",
                        help="Skip the zcat|grep prefilter and read the CSV directly")

    args = parser.parse_args()
    started = time.time()
    log.info("Starting parse_domain_architectures")

    family_to_clan, clan_to_rep = load_clan_membership(args.clan_membership)
    pfam_mapping = load_pfam_mapping(args.pfam_mapping)

    counters = count_architectures(
        iter_metadata(args.proteins, use_prefilter=not args.no_prefilter),
        family_to_clan, clan_to_rep, pfam_mapping, args.base_url, args.overlap_fraction,
        args.log_every)

    write_outputs(counters, set(family_to_clan), args.output_dir, args.top, clan_to_rep,
                  pfam_mapping, args.base_url)

    log.info(f"parse_domain_architectures complete in {time.time() - started:.1f}s")


if __name__ == "__main__":
    main()
