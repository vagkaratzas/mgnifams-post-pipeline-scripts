# Conventions

- Root `bin/*.py` scripts usually expose small pure functions plus `parse_args`, `output_paths` when needed, and `main()` guarded by `if __name__ == "__main__"`.
- Prefer narrow, script-local helpers over shared abstractions for one-off post-processing tools.
- Keep CLI defaults in module constants near the top; write outputs under `output/` unless script-specific docs say otherwise.