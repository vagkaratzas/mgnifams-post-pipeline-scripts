#!/usr/bin/env python3
"""Rebuild the released MGnifams full MSAs with the released HMMs.

Why this exists
---------------
Legacy `bin/generate_families.py` exports a family HMM built with `hmmbuild --hand`
from the final seed alignment (`:547`), but recruits (`:549`) and aligns (`:574`) the
full MSA with `hmm` -- the `fast`-architecture model of the previous round. The shipped
HMM is therefore not the model that produced the shipped full MSA. On the workspace's
own 9-family output, 4 full MSAs have a match-column count that differs from their
HMM's `LENG`.

This re-recruits and re-aligns every family with the *shipped* HMM, which is the
artefact users and the website actually have.

Writes, and nothing else -- it never touches a database. `update_mgnifam_db_full_msa_bug.py` applies
the CSV this produces.

What it deliberately reuses
---------------------------
Search, hit filtering, alignment and renumbering come from `mgnifam.generate_families`
(the fixed standalone package), so the regenerated MSAs cannot drift from the pipeline's
own semantics. The exit-branch recruitment rule is reproduced exactly: E-value cutoff
only, envelope-length filter waived (`exit_flag=True`), envelopes masked out of their
record, rows in HMMER ranking order. The representative columns are derived exactly as
`renumber_sto_msa` derived them: row 0, ungapped and upper-cased.

What it does NOT do
-------------------
* It does not re-run the seed/redundancy loop. Seed MSAs, RF lines, HMMs and therefore
  the `consensus` and `converged` columns are unchanged.
* It does not re-apply the membership check (`--discard_min_starting_membership`): the
  original cluster members are not recoverable from an HMM. Families that the corrected
  model would have discarded stay released; the report's `rep_length` lets outliers be
  spotted.
* Family membership changes, so the `biome_blob` and `domain_blob` columns -- built by
  `update_mgnifams_db` from `refined_families.tsv` -- go stale for any family whose
  `full_size` moved. Regenerating those is a separate pipeline run.

Outputs
-------
    <outdir>/full_msa/<id>.sto[.gz]   rebuilt alignments
    <outdir>/mgnifam_updates.csv      DB-ready values, plus diagnostics

Usage
-----
    /home/vangelis/Desktop/Projects/mgnifam/.venv/bin/python bin/rebuild_full_msas.py \
        --hmm_lib mgnifams_hmm.lib.gz \
        --fasta output/mgnifams_v2.fa \
        --outdir full_msa_fixed \
        --old_full_msa_dir output/generate_families/families/full_msa_sto \
        --cpus 8

`--fasta` must be the pipeline's own search database (`SETUP_CLUSTERS.out.mgnifams_input_fa`,
e.g. `output/mgnifams_v2.fa`), uncompressed, not the raw MGnify protein FASTA: recruitment
is only comparable against the database the release was built from.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import importlib.util
import itertools
import logging
import re
import sys
import time
from importlib.metadata import PackageNotFoundError
from pathlib import Path

# `mgnifam.generate_families` is 3.13-only (`itertools.batched(strict=)`, bare
# `Generator[T]`), and so is this script by extension.
if sys.version_info < (3, 13):
    raise SystemExit(
        "python >= 3.13 required; run with mgnifam's interpreter, e.g.\n"
        "  /home/vangelis/Desktop/Projects/mgnifam/.venv/bin/python bin/rebuild_full_msas.py ..."
    )

DEFAULT_MGNIFAM_SRC = Path("/home/vangelis/Desktop/Projects/mgnifam/src")

# The `mgnifam` table columns this rebuild invalidates. `consensus` and `converged` are
# not here: both come from the HMM and the seed loop, neither of which is re-run.
DB_COLUMNS = ("full_size", "protein_rep", "rep_region", "rep_length", "rep_sequence")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--hmm_lib", type=Path, help="released HMM library (plain or .gz)")
    parser.add_argument("--fasta", type=Path, help="uncompressed pipeline search database")
    parser.add_argument("--outdir", type=Path, help="output directory")
    parser.add_argument(
        "--old_full_msa_dir", type=Path, default=None, help="released full MSAs, for the diff columns"
    )
    parser.add_argument("--ssi", type=Path, default=None, help="SSI index path (default: <fasta>.ssi)")
    parser.add_argument("--recruit_evalue_cutoff", type=float, default=0.001)
    parser.add_argument("--batch_size", type=int, default=1000, help="HMMs per database pass")
    parser.add_argument("--cpus", type=int, default=0, help="0 = all available")
    parser.add_argument("--prefetch_targets", action="store_true", help="hold the database in RAM")
    parser.add_argument("--gzip_output", action="store_true", help="write .sto.gz instead of .sto")
    parser.add_argument("--dry_run", action="store_true", help="report only, write no alignments")
    parser.add_argument("--mgnifam_src", type=Path, default=DEFAULT_MGNIFAM_SRC, help="path to mgnifam's src/")
    parser.add_argument("--self_test", action="store_true", help="run the built-in check and exit")
    return parser.parse_args(argv)


def load_generate_families(src: Path):
    """Import `mgnifam.generate_families`, falling back to loading the file directly.

    The fallback bypasses `mgnifam/__init__.py`, which reads installed package metadata
    and so fails on a plain `sys.path` insert of an uninstalled checkout.
    """
    try:
        from mgnifam import generate_families
    except (ModuleNotFoundError, PackageNotFoundError):
        module_path = src / "mgnifam" / "generate_families.py"
        if not module_path.exists():
            raise SystemExit(f"{module_path} not found; pass --mgnifam_src")
        spec = importlib.util.spec_from_file_location("mgnifam_generate_families", module_path)
        generate_families = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(generate_families)
    return generate_families


def text(value: str | bytes) -> str:
    """Normalise an Easel name. pyhmmer hands back `str` or `bytes` by version."""
    return value.decode() if isinstance(value, bytes) else value


def read_msa_rows(path: Path) -> tuple[list[str], str]:
    """Return (row names, RF line) of a Stockholm file, reading .gz transparently."""
    opener = gzip.open if path.suffix == ".gz" else open
    names, reference = [], []
    with opener(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if not fields or line.startswith("//"):
                continue
            if line.startswith("#=GC RF"):
                reference.append(fields[2])
            elif not line.startswith("#"):
                names.append(fields[0])
    return names, "".join(reference)


def write_msa(msa, path: Path, gzip_output: bool) -> None:
    opener = gzip.open if gzip_output else open
    with opener(path, "wb") as handle:
        msa.write(handle, format="pfam")


def representative_columns(msa) -> dict[str, object]:
    """Derive the `mgnifam` rep_* values from row 0, as `renumber_sto_msa` did.

    Row 0 is the family representative: hits arrive in HMMER's ranking order, so the
    best-scoring sequence leads the alignment. A row spanning a whole protein carries no
    `/start-end` suffix, and the legacy convention records its region as `-`.
    """
    name = text(msa.names[0])
    residues = re.sub(r"[.\-~]", "", text(msa.alignment[0])).upper()
    protein, _, region = name.partition("/")
    return {
        "protein_rep": protein,
        "rep_region": region or "-",
        "rep_length": len(residues),
        "rep_sequence": residues,
    }


def rebuild_family(gf, hmm, hits, indexed, cpus: int):
    """Re-recruit and re-align one family with its own shipped HMM.

    Mirrors the legacy exit branch: `exit_flag=True` waives the envelope-length filter,
    so the length percentage passed here never applies and is fixed at 0.0.
    """
    records = gf.extract_records(hits)
    sequences = gf.filter_hits(records, hits.query.M, True, 0.0, indexed)
    if not sequences:
        return None, "no hits above the E-value cutoff"
    msa = gf.renumber_msa(gf.run_hmmalign(hmm, sequences, cpus), text(hmm.name), indexed)

    # Duplicate row names are illegal Stockholm. The fixed `parse_protein_name` resolves
    # repeat domains to distinct coordinates, so this should never fire -- it is kept
    # because the legacy renumbering had to drop duplicates and the counts must be
    # comparable if it ever does.
    seen: set[str] = set()
    keep = [text(name) not in seen and not seen.add(text(name)) for name in msa.names]
    if not all(keep):
        import pyhmmer

        deduplicated = pyhmmer.easel.TextMSA(
            name=msa.name,
            sequences=[row for row, wanted in zip(msa.sequences, keep) if wanted],
        )
        deduplicated.reference = msa.reference
        msa = deduplicated
    return msa, None


def compare_with_release(old_dir: Path, family: str, msa, hmm) -> dict[str, object]:
    """Diff the rebuilt MSA against the released one, if it is there."""
    for candidate in (old_dir / f"{family}.sto", old_dir / f"{family}.sto.gz"):
        if candidate.exists():
            break
    else:
        return {}
    old_names, old_reference = read_msa_rows(candidate)
    return {
        "old_full_size": len(old_names),
        "old_match_columns": old_reference.count("x"),
        "old_hmm_mismatch": "yes" if old_reference.count("x") != hmm.M else "no",
        "old_rep": old_names[0] if old_names else "",
        "size_changed": "yes" if len(old_names) != len(msa.names) else "no",
        "rep_changed": "yes" if (old_names[0] if old_names else None) != text(msa.names[0]) else "no",
    }


def write_report(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        fieldnames.extend(key for key in row if key not in fieldnames)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    options = parse_args(argv)
    if options.self_test:
        return self_test(options)
    for required in ("hmm_lib", "fasta", "outdir"):
        if getattr(options, required) is None:
            raise SystemExit(f"--{required} is required")

    gf = load_generate_families(options.mgnifam_src)
    import pyhmmer

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger("rebuild_full_msas")

    msa_dir = options.outdir / "full_msa"
    msa_dir.mkdir(parents=True, exist_ok=True)
    updates_path = options.outdir / "mgnifam_updates.csv"

    ssi_path = options.ssi or Path(f"{options.fasta}.ssi")
    if not ssi_path.exists():
        logger.info("building SSI index at %s", ssi_path)
        gf.build_ssi_index(options.fasta, ssi_path)

    started = time.monotonic()
    rows: list[dict[str, object]] = []
    mismatched_after = 0

    with contextlib.ExitStack() as stack:
        index_reader = stack.enter_context(pyhmmer.easel.SSIReader(ssi_path))
        indexed = gf.IndexedSequences(
            stack.enter_context(
                pyhmmer.easel.SequenceFile(options.fasta, digital=False, index=index_reader)
            )
        )
        target_file = stack.enter_context(
            pyhmmer.easel.SequenceFile(options.fasta, digital=True, alphabet=gf.ALPHABET)
        )
        targets = target_file.read_block() if options.prefetch_targets else target_file
        hmm_file = stack.enter_context(pyhmmer.plan7.HMMFile(options.hmm_lib))

        for batch_number, batch in enumerate(itertools.batched(hmm_file, options.batch_size), 1):
            hmms = list(batch)
            with gf.search(
                hmms,
                targets,
                cpus=options.cpus,
                evalue=options.recruit_evalue_cutoff,
                logger=logger,
                batch_number=batch_number,
                round_number=0,
            ) as searched:
                for hmm, hits in zip(hmms, searched, strict=True):
                    family = text(hmm.name)
                    msa, failure = rebuild_family(gf, hmm, hits, indexed, options.cpus)
                    row: dict[str, object] = {"id": family, "status": "rebuilt", "hmm_length": hmm.M}
                    if msa is None:
                        row.update(status="skipped", note=failure)
                        rows.append(row)
                        logger.warning("family %s: %s", family, failure)
                        continue

                    reference = msa.reference or ""
                    row["full_size"] = len(msa.names)
                    row.update(representative_columns(msa))
                    row["match_columns"] = reference.count("x")
                    row["note"] = ""
                    # The whole point of the rebuild: the shipped HMM must be the model the
                    # shipped alignment came from, so its match states must be the alignment's
                    # RF columns. A failure here means the rebuild did not do what it claims.
                    if row["match_columns"] != hmm.M:
                        mismatched_after += 1
                        row["note"] = "RF columns still differ from HMM length"
                        logger.error("family %s: RF %s != LENG %s", family, row["match_columns"], hmm.M)

                    if options.old_full_msa_dir is not None:
                        row.update(compare_with_release(options.old_full_msa_dir, family, msa, hmm))
                    if not options.dry_run:
                        suffix = ".sto.gz" if options.gzip_output else ".sto"
                        write_msa(msa, msa_dir / f"{family}{suffix}", options.gzip_output)
                    rows.append(row)

    write_report(updates_path, rows)
    rebuilt = [row for row in rows if row["status"] == "rebuilt"]
    logger.info(
        "rebuilt=%d skipped=%d still_mismatched=%d elapsed=%.1fs updates=%s",
        len(rebuilt),
        len(rows) - len(rebuilt),
        mismatched_after,
        time.monotonic() - started,
        updates_path,
    )
    if options.old_full_msa_dir is not None:
        logger.info(
            "size_changed=%d rep_changed=%d",
            sum(1 for row in rebuilt if row.get("size_changed") == "yes"),
            sum(1 for row in rebuilt if row.get("rep_changed") == "yes"),
        )
    logger.info("apply with: update_mgnifam_db_full_msa_bug.py --csv %s --sqlite <db> --apply", updates_path)
    return 1 if mismatched_after else 0


def self_test(options: argparse.Namespace) -> int:
    """Rebuild a toy family end to end and assert the invariant the fix is about.

    Builds an HMM by hand from a small alignment, plants its sequences in a FASTA, then
    runs the same rebuild path and checks the alignment's RF columns match the HMM and
    that the representative columns come out in `mgnifam`'s shape.
    """
    import tempfile

    gf = load_generate_families(options.mgnifam_src)
    import pyhmmer

    records = [
        ("p1_1_39", "MKVLAAGIVGLNLGGSTVAAMKVLAAGIVGLNLGGSTV"),
        ("p2_1_39", "MKVLAAGIVGLNLAGSTVAAMKVLAAGIVGLNLGGSTV"),
        ("p3_1_39", "MKVLAAGIVGLNLGGSTVAAMKVLAAGIVGLNLGGSTA"),
    ]
    with tempfile.TemporaryDirectory() as workdir:
        work = Path(workdir)
        fasta = work / "db.fa"
        fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))

        msa = pyhmmer.easel.TextMSA(
            name=b"1",
            sequences=[
                pyhmmer.easel.TextSequence(name=name.encode(), sequence=sequence)
                for name, sequence in records
            ],
        ).digitize(gf.ALPHABET)
        hmm = gf.run_hmmbuild(msa, "1")

        ssi = work / "db.fa.ssi"
        gf.build_ssi_index(fasta, ssi)
        with contextlib.ExitStack() as stack:
            reader = stack.enter_context(pyhmmer.easel.SSIReader(ssi))
            indexed = gf.IndexedSequences(
                stack.enter_context(pyhmmer.easel.SequenceFile(fasta, digital=False, index=reader))
            )
            targets = stack.enter_context(
                pyhmmer.easel.SequenceFile(fasta, digital=True, alphabet=gf.ALPHABET)
            ).read_block()
            with gf.search(
                [hmm], targets, cpus=1, evalue=0.001,
                logger=logging.getLogger("self_test"), batch_number=1, round_number=0,
            ) as searched:
                hits = next(iter(searched))
                rebuilt, failure = rebuild_family(gf, hmm, hits, indexed, 1)

        assert failure is None, failure
        assert rebuilt is not None
        assert (rebuilt.reference or "").count("x") == hmm.M, "RF columns must equal the HMM length"
        assert len(rebuilt.names) == len(records), "every planted sequence must be recruited"
        assert len({text(name) for name in rebuilt.names}) == len(rebuilt.names), "names must be unique"

        columns = representative_columns(rebuilt)
        assert set(columns) == set(DB_COLUMNS) - {"full_size"}
        assert columns["rep_length"] == len(columns["rep_sequence"])
        assert columns["rep_sequence"].isupper()
        assert re.fullmatch(r"\d+-\d+|-", str(columns["rep_region"])), columns["rep_region"]

        updates = work / "mgnifam_updates.csv"
        write_report(updates, [{"id": "1", "status": "rebuilt", "full_size": len(rebuilt.names), **columns}])
        written = list(csv.DictReader(updates.open()))
        assert written[0]["rep_sequence"] == columns["rep_sequence"]

    print("self_test: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
