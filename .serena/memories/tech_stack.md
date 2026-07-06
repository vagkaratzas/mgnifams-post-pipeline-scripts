# Tech Stack

- Python scripts target the local system Python and mostly use stdlib plus common scientific stack where needed (`pandas`, `numpy`, `matplotlib`, `biopython`).
- Tests are pytest tests under `tests/` and import standalone scripts via `importlib.util.spec_from_file_location`.
- SQLite data is accessed with stdlib `sqlite3`; avoid loading large DB-derived blobs into pandas inside tight loops unless a vectorized table operation is genuinely needed.