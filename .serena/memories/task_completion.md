# Task Completion

- For script changes, run the focused pytest file first, e.g. `python -m pytest tests/test_biome_analysis.py -q`.
- Also run the changed script's `--help` if CLI arguments changed.
- Inspect `git diff -- <changed files>` before final response and mention any checks that could not be run.