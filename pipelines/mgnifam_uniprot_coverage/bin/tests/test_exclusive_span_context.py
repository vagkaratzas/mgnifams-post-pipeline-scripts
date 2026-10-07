import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "exclusive_span_context.py"
spec = importlib.util.spec_from_file_location("ctx", SCRIPT)
ctx = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctx)


def test_classify_partitions_exclusive_residues():
    mgnifam = {
        "a": [(1, 50), (40, 120)],  # one block 1-120 over Pfam 61-100: 60 + 20 extension
        "b": [(101, 110)],          # abuts Pfam 1-100: extension
        "c": [(300, 400)],          # Pfam 1-100 elsewhere: separate
        "d": [(1, 80)],             # no Pfam hit: no_pfam
        "e": [(10, 20)],            # inside Pfam: nothing exclusive
    }
    pfam = {"a": [(61, 100)], "b": [(1, 100)], "c": [(1, 100)], "e": [(1, 100)]}
    assert dict(ctx.classify(mgnifam, pfam)) == {
        ("extension", ">=75"): [80, 1],
        ("extension", "<30"): [10, 1],
        ("separate", ">=75"): [101, 1],
        ("no_pfam", ">=75"): [80, 1],
    }
