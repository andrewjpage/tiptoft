# Tiptoft

## Quick Reference
- **Language:** Python 3.8+
- **Test runner:** `python -m pytest tiptoft/tests/ -v` (54 tests)
- **Dev environment:** Package not installed — use `PYTHONPATH=. python scripts/tiptoft` to run scripts
- **Check for deprecation warnings:** `python -W error::DeprecationWarning -c "from tiptoft.Kmers import Kmers"`

## Architecture
- `tiptoft/` — main package (TipToft orchestrator, Kmers, Fasta, Fastq, Gene, Blocks, Read)
- `scripts/` — CLI entry points (`tiptoft`, `tiptoft_database_downloader`)
- `homopolymer_compression.py` — pure Python fallback for optional Cython extension (.pyx/.c)
- `tiptoft/data/` — bundled plasmid database (plasmid_data.fa)

## Code Style
- No type annotations in existing code — don't add them
- Docstrings use Google-style format (Args/Returns/Raises)
- Class naming: PascalCase; method naming: snake_case

## Key Patterns
- Bundled data accessed via `importlib.resources.files('tiptoft').joinpath('data', ...)`
- Version from `importlib.metadata.version("tiptoft")` with `PackageNotFoundError` fallback
- Homopolymer compression: import prefers compiled `.so` over `.py` fallback automatically
- setup.py uses Cython conditionally (`USE_CYTHON` flag)
