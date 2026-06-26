#!/usr/bin/env python3
"""Rank profile HMMs from a HMMER `hmmstat` report by a composite quality score.

`hmmstat` emits one row per profile with these columns:

    idx  name  accession  nseq  eff_nseq  M  relent  info  p relE  compKL

This tool reads that report, computes a composite score that ranks the profiles
from "best" to "worst", and prints the original columns plus the score, sorted
best-first.

What makes a profile HMM "good"?
--------------------------------
There is no single ground-truth quality number in `hmmstat`, so the composite
score combines the columns that genuinely reflect how trustworthy and useful a
profile is for searching:

  * eff_nseq  (REWARD, heaviest weight)
        The effective number of sequences after HMMER's entropy weighting.
        This is the best available proxy for statistical support and family
        diversity. A profile built from many effectively-independent sequences
        generalises; one built from 1-2 near-identical sequences (eff_nseq < 1)
        is fragile.

        Note: eff_nseq is used *raw*, not divided by nseq. The ratio eff_nseq/nseq
        would measure training-data diversity ("what fraction of sequences were
        independent?"), not profile quality. A profile with eff_nseq=7 (from 100
        sequences) has far better statistical support than one with eff_nseq=1.8
        (from 2 sequences); the ratio 0.07 vs 0.9 would invert that ranking.
        eff_nseq is already the entropy-corrected measure HMMER uses for model
        calibration — it is the right metric as-is.

  * p relE    (REWARD)
        Mean positional relative entropy per match state in bits, including
        transition (insertion/deletion) probabilities. The HMMER manual describes
        this as "a fancier version of the per-match-state relative entropy" that
        "may be a more accurate estimation of the average score contributed per
        model consensus position." It is preferred over `relent` (which is pinned
        to the entropy-weighting target ~0.59 for default models and therefore
        has low variance) and over `info` (which the HMMER manual flags as
        "probably not useful").

  * M         (REWARD, light weight)
        Model length in consensus match states. Longer models accumulate score
        over more positions and therefore tend to yield more statistically
        significant hits, all else being equal.

  * compKL    (PENALTY)
        KL divergence of the model's residue composition from the background.
        The HMMER manual notes that highly biased profiles can slow the
        acceleration pipeline by causing too many non-homologous sequences to
        pass the filters, so it is penalised.

Note on relent / info: both are excluded. relent is the entropy-weighting target
itself (~0.59 bits for default HMMER3 models), giving it very low variance and
little ranking power. info is explicitly described as "probably not useful" in
the HMMER manual. p relE supersedes both.

The formula
-----------
Each contributing column is min-max normalised to [0, 1] across the profiles in
the input, so columns on very different scales (eff_nseq ~0-7, M ~40-520,
p_relE ~0.5-1.2, compKL ~0.01-0.22) contribute comparably. The composite is a
weighted sum of the rewards minus the weighted penalty:

    score =  w_eff   * norm(eff_nseq)
           + w_prelE * norm(p_relE)
           + w_M     * norm(M)
           - w_kl    * norm(compKL)

Default weights (rewards sum to 1.0): w_eff=0.45, w_prelE=0.30, w_M=0.15,
w_kl=0.10. All four are overridable on the command line.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field

# Order of the data columns in an hmmstat report.
COLUMNS = ["idx", "name", "accession", "nseq", "eff_nseq", "M",
           "relent", "info", "p_relE", "compKL"]


@dataclass
class Profile:
    idx: int
    name: str
    accession: str
    nseq: int
    eff_nseq: float
    M: int
    relent: float
    info: float
    p_relE: float
    compKL: float
    score: float = field(default=0.0)


def parse_hmmstat(handle) -> list[Profile]:
    """Parse an hmmstat report into a list of Profile records.

    Comment/header lines (starting with '#') and blank lines are skipped. Each
    data row is expected to have 10 whitespace-separated fields.
    """
    profiles: list[Profile] = []
    for lineno, raw in enumerate(handle, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != len(COLUMNS):
            raise ValueError(
                f"line {lineno}: expected {len(COLUMNS)} columns, "
                f"got {len(fields)}: {line!r}"
            )
        profiles.append(Profile(
            idx=int(fields[0]),
            name=fields[1],
            accession=fields[2],
            nseq=int(fields[3]),
            eff_nseq=float(fields[4]),
            M=int(fields[5]),
            relent=float(fields[6]),
            info=float(fields[7]),
            p_relE=float(fields[8]),
            compKL=float(fields[9]),
        ))
    return profiles


def _normaliser(values: list[float]):
    """Return a min-max normaliser mapping the given values onto [0, 1].

    If every value is identical (max == min) the column carries no ranking
    information, so the normaliser returns 0.0 for all inputs (neutral).
    """
    lo, hi = min(values), max(values)
    span = hi - lo
    if span == 0:
        return lambda _x: 0.0
    return lambda x: (x - lo) / span


def score_profiles(profiles: list[Profile], w_eff: float, w_prelE: float,
                   w_M: float, w_kl: float) -> None:
    """Compute the composite score for each profile, in place."""
    if not profiles:
        return
    norm_eff = _normaliser([p.eff_nseq for p in profiles])
    norm_prelE = _normaliser([p.p_relE for p in profiles])
    norm_M = _normaliser([float(p.M) for p in profiles])
    norm_kl = _normaliser([p.compKL for p in profiles])

    for p in profiles:
        p.score = (
            w_eff * norm_eff(p.eff_nseq)
            + w_prelE * norm_prelE(p.p_relE)
            + w_M * norm_M(float(p.M))
            - w_kl * norm_kl(p.compKL)
        )


def write_table(profiles: list[Profile], handle, tsv: bool = False) -> None:
    """Write ranked profiles (best first) with a rank column and the score."""
    header = ["rank"] + COLUMNS + ["score"]

    def fmt(p: Profile, rank: int) -> list[str]:
        return [
            str(rank), str(p.idx), p.name, p.accession, str(p.nseq),
            f"{p.eff_nseq:.2f}", str(p.M), f"{p.relent:.2f}", f"{p.info:.2f}",
            f"{p.p_relE:.2f}", f"{p.compKL:.2f}", f"{p.score:.4f}",
        ]

    rows = [fmt(p, rank) for rank, p in enumerate(profiles, start=1)]

    if tsv:
        print("\t".join(header), file=handle)
        for row in rows:
            print("\t".join(row), file=handle)
        return

    widths = [len(h) for h in header]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def render(cells: list[str]) -> str:
        # Left-justify the text identifier columns, right-justify the numbers.
        left = {1, 2, 3}  # idx, name, accession positions within `header`
        out = []
        for i, cell in enumerate(cells):
            out.append(cell.ljust(widths[i]) if i in left
                       else cell.rjust(widths[i]))
        return "  ".join(out)

    print(render(header), file=handle)
    print("  ".join("-" * w for w in widths), file=handle)
    for row in rows:
        print(render(row), file=handle)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Rank profile HMMs from a HMMER hmmstat report by a "
                    "composite quality score.")
    parser.add_argument(
        "input", nargs="?", default="-",
        help="hmmstat report file (default: stdin)")
    parser.add_argument(
        "-o", "--output", default="-",
        help="output file (default: stdout)")
    parser.add_argument(
        "--tsv", action="store_true",
        help="emit tab-separated output instead of an aligned table")
    parser.add_argument("--w-eff", type=float, default=0.45,
                        help="weight for eff_nseq (default: 0.45)")
    parser.add_argument("--w-prelE", type=float, default=0.30,
                        help="weight for p relE (default: 0.30)")
    parser.add_argument("--w-m", type=float, default=0.15,
                        help="weight for M / model length (default: 0.15)")
    parser.add_argument("--w-kl", type=float, default=0.10,
                        help="penalty weight for compKL (default: 0.10)")
    args = parser.parse_args(argv)

    in_handle = sys.stdin if args.input == "-" else open(args.input)
    try:
        profiles = parse_hmmstat(in_handle)
    finally:
        if in_handle is not sys.stdin:
            in_handle.close()

    if not profiles:
        print("error: no data rows found in input", file=sys.stderr)
        return 1

    score_profiles(profiles, args.w_eff, args.w_prelE, args.w_m, args.w_kl)
    profiles.sort(key=lambda p: p.score, reverse=True)

    out_handle = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        write_table(profiles, out_handle, tsv=args.tsv)
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
