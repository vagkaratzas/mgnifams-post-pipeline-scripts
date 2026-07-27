#!/usr/bin/env python3

"""Build per-family domain architecture JSONs from the MGnify re-annotation CSV.

See SPEC.md for the frozen behaviour. Replaces the pipeline's parse_domains.py, which read
per-family TSVs from the old MGnify protein database and assumed one MGnifam hit per sequence.
"""

import csv
import logging
from collections import namedtuple

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
    """Return accession -> human readable name."""
    pfam_mapping = {}

    with open(pfam_mapping_file) as handle:
        for line in handle:
            accession, name = line.rstrip("\n").split("\t", 1)
            pfam_mapping[accession] = name

    log.info(f"Loaded {len(pfam_mapping)} pfam entries")

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
