# AGENTS.md

Orientation and guardrails for agents working in the SURFA repository.

## What this project is

SURFA (Small Upstream Reading Frame Annotator) is a Python 3.13+ CLI and library that
scans Ensembl **GTF + FASTA** inputs and writes a SQLite database of **5' UTR upstream
open reading frames (uORFs)**. A second command exports one transcript's annotations as
JSON.

Key facts:

- CLI entry point: `surfa` (defined in `pyproject.toml`, implemented in `src/surfa/cli.py`).
- Runtime dependency on a locally installed **`bedtools`** binary. The Python dependency
  set is intentionally small (`numpy`, `pandas`).
- Output database: `<output-dir>/uorfs.db`.
- Primary package code lives in `src/surfa/`; tests in `tests/`.
- Docs are published with **Zensical** (a MkDocs-style tool) from `docs/`.
- Current package version is `0.1.0` (see version caveat below).

## Repository map

```
README.md                  User-facing install/usage; mirrors docs/index.md
docs/
  index.md                 Installation, build/query usage, seqid mapping
  annotations.md           Biological model: 5' UTR, CDS, N-terminal extensions,
                           missing CDS, near-cognate start codons
  schema.md                Human-readable database schema (may lag schema.json)
  api.md                   mkdocstrings API pages for three modules
src/surfa/
  __init__.py              __version__
  __main__.py              Empty (see gotchas)
  cli.py                   argparse root parser, logging setup, subcommand dispatch
  commands/
    build.py               `surfa build` parser + run()
    query.py               `surfa query` parser + run()
  uorf_utils.py            Core pipeline: GTF parsing -> UTR extraction -> uORF calling
                           -> SQLite export (gtf_to_uorf_db)
  fasta_utils.py           BED creation, bedtools getfasta calls, strand handling
  db_utils.py              schema loading, DataFrame validation/coercion, SQLite writes
  json_converter.py        Query DB tables and assemble/export transcript JSON
  schema.json              Source of truth for SQLite tables/types/nullability
tests/
  unit/                    Fast unittest suite for check_identity/fasta_codon_search
  integration/             DB build script + sequence-identity tests (Docker-oriented)
  data/                    Mini GTF/FASTA, NRG1/MEF2C reference sequences
data/seqid_map.csv         Example human chromosome/GenBank/RefSeq mapping table
Dockerfile                 Ubuntu + uv + bedtools image used by CI
.github/workflows/main.yml CI: build image, build test DB, run unittest discovery
.github/workflows/docs.yml Deploy docs on pushes to main/master
```

## Setup

### Option A: local development

```bash
brew install bedtools          # or use your platform's package manager
uv sync
source .venv/bin/activate
```

`uv sync` installs the project plus the `dev` dependency group (ruff, zensical,
mkdocstrings, mdformat, pip-audit). The `surfa` console script is then available.

### Option B: Docker (closest to CI)

```bash
docker build -t surfa-test .
```

The container installs bedtools, `uv sync --frozen --no-dev`, and working directory is
`/app`.

## Common commands

```bash
# Build a database. --output-dir MUST already exist.
surfa build \
  --gtf tests/data/mini.gtf.gz \
  --fasta tests/data/minifasta.fa \
  --output-dir /tmp/surfa-build \
  --ensembl-source ensembl_havana

# Export one transcript. Transcript ID must include the version.
surfa query \
  --db /tmp/surfa-build/uorfs.db \
  --transcript ENST00000504921.7 \
  --output /tmp/uorfs.json
```

Global flags available on every command: `-l/--log-level` and `--log-file`
(default `surfa.log` in the current directory). Subcommand help:
`surfa build --help`, `surfa query --help`.

### Testing

CI runs integration and unit tests **inside the Docker image**:

```bash
docker run --rm surfa-test bash -c \
  "tests/integration/build_test_db.sh && python -m unittest discover -s tests -p '*test*.py' -v"
```

Local unit tests (no bedtools or fixture DB required):

```bash
uv run python -m unittest discover -s tests/unit -p '*test*.py' -v
```

Local integration tests are possible only after producing `/app/tests/uorfs.db` and the
`/app/tests/data/...` fixtures expected by the hardcoded paths, or by adapting those
paths. Prefer Docker when touching sequence/strand behavior. See "Gotchas".

Tests use **`unittest`**, not pytest. There is no `conftest.py`.

### Linting

Run Ruff before proposing a change (the PR template requires it). There is no checked-in
Ruff config, so defaults apply:

```bash
uv run ruff check .
uv run ruff format .          # only on files you changed; avoid repo-wide churn
```

### Docs

```bash
uv run zensical build --clean   # outputs site/
uv run zensical serve           # local preview
```

API docs are generated from Google-style docstrings via mkdocstrings (`paths = ["src"]`).

## Architecture and data flow

`surfa build` -> `commands/build.py:run` -> `uorf_utils.gtf_to_uorf_db`:

1. **Read GTF** with pandas (tab-separated, comment `#`, named GTF columns).
2. **Filter to `--ensembl-source`** (`source` column), e.g. `ensembl_havana` or `ensembl`.
3. **CDS**: take the first CDS exon per transcript (lowest exon number) as the coding
   reference. Its genomic `start` defines the CDS start boundary.
4. **5' UTR extraction**:
   - Newer GTFs use `feature == "five_prime_utr"`; exon numbers are assigned by
     transcript order (ascending start on `+`, descending on `-`).
   - Older GTFs use `feature == "UTR"` and `exon_number`; 5'/3' identity is inferred by
     `check_identity(...)` relative to CDS start.
   - Transcripts without a mapped CDS start are dropped (logged).
5. **Frame state**: `sum(UTR exon lengths) % 3` per transcript, stored as
   `utr.frame_state`.
6. **Sequence retrieval**: UTR and first-CDS intervals become BED records
   (`fasta_utils.make_bed`, converting 1-based inclusive GTF to 0-based BED), then
   `bedtools getfasta` is invoked via subprocess. Optional seqid mapping remaps GTF
   chromosome names to FASTA names before the call.
7. **Transcript FASTA**: per-transcript exon sequences are concatenated in exon order;
   negative-strand transcripts are reverse-complemented (`get_transcript_FASTA`).
8. **uORF calling** (`get_uorfs`): each transcript is copied into frames 0/1/2; DNA is
   converted to RNA (`T` -> `U`); codons are extracted; start and stop codons are matched
   to the next downstream stop; sequences and annotations are generated.
9. **IDs and counts**: sort and construct `uorf_id = <transcript>.<start_codon>.<n>`;
   count uORFs per transcript.
10. **SQLite export**: `db_utils.write_to_db` loads `schema.json`, validates/coerces each
    DataFrame (`validate_and_coerce`), creates tables, and inserts rows. Writes are
    transactional; foreign keys are enabled.

`surfa query` -> `json_converter.assemble_json_from_transcript` reads the `transcripts`,
`utr`, `cds`, and `uorfs` tables for one transcript and emits a JSON structure with
transcript and uORF region blocks.

### Module responsibilities

| Module | Responsibility |
|---|---|
| `cli.py` | Parser construction, logging, dispatch. Commands expose `register()` and `run()`. |
| `commands/build.py` | CLI argument surface for database builds. |
| `commands/query.py` | CLI argument surface for JSON exports. |
| `uorf_utils.py` | GTF parsing, UTR/CDS logic, codon matching, uORF annotations, build orchestration. |
| `fasta_utils.py` | BED conversion, bedtools subprocess integration, seqid remapping, strand handling. |
| `db_utils.py` | Schema-driven validation, type coercion, table creation, SQLite writes, md5 metadata. |
| `json_converter.py` | DB queries and JSON formatting/export. |
| `schema.json` | Authoritative table/column definitions; packaged as package data. |

## Domain invariants (do not change casually)

These encode biology and coordinate semantics. Changes require a clear rationale,
updated docs, and test coverage.

- **Only 5' UTR uORFs are called.** 3' UTR and other regions are out of scope.
- **Start codons** are the canonical `AUG` plus near-cognate codons:
  `CUG, GUG, UUG, ACG, AAG, AGG, AUC, AUU, AUA`. **Stop codons** are `UAA, UAG, UGA`.
  The docs justify near-cognate inclusion; removing codons changes output semantics.
- **Sequences stored in the DB are mixed alphabet**: `utr.FASTA` and `cds.FASTA` are DNA
  (`T`); `transcripts.transcript_FASTA` comes from bedtools and is DNA; uORF sequences
  (`uorfs.uORF_FASTA`, and the working transcript sequence) are RNA (`U`). Preserve this
  convention or update all consumers, including tests and JSON output.
- **Coordinate systems**:
  - GTF `start`/`end`: 1-based, inclusive.
  - BED `chromStart`: 0-based, half-open (`make_bed` subtracts 1).
  - `rel_start`/`rel_stop`, `rel_start_pos`/`rel_stop_pos`: 0-based indices into the
    relevant FASTA. uORF `rel_stop_pos` is the exclusive end of the stop codon.
- **Negative-strand handling** relies on per-transcript exon concatenation followed by
  reverse complementation. Integration tests (`tests/integration/test_sequence_identity.py`)
  check NRG1 (positive) and MEF2C (negative) against reference sequences.
- **`NO_UTR_STOP`** means no in-UTR stop was found and the uORF reads into the CDS. For
  those rows, `rel_stop_pos` and `uorf_length` are set to the transcript length.
- **`CDS_frame_status`** is computed as `(utr_len - frame) % 3 == 0 -> "in-frame"`, else
  `"out-frame"`; `uorf_location` is `"UTR_only"` or `"CDS_overlap"`.
- **Known documented limitations**: transcripts without a defined CDS are dropped;
  N-terminal extensions of the CDS are not called (assumed represented by alternate
  transcripts). See `docs/annotations.md`.

## Coding conventions

- Python 3.13+, `from __future__ import annotations` in modules that use modern hints.
- Keep the existing style: module-level `logger = logging.getLogger(__name__)`, small
  focused helpers, `__all__` for public module APIs.
- Docstrings are **Google style** (mkdocstrings renders them with `docstring_style = "google"`).
  New public functions should have parameter/return documentation.
- Arguments are passed as explicit positional parameters through the pipeline
  (`gtf_to_uorf_db(...)`, `gtf_to_sequence(...)`), not via globals.
- pandas is used heavily for table operations; be careful with dtype changes, index
  preservation, and in-place mutations (`inplace=True` is common).
- New SQLite tables/columns require updating `src/surfa/schema.json` **and**
  `docs/schema.md`, plus any DataFrame construction in `uorf_utils.py`.
- Add fast logic tests under `tests/unit/` (unittest). Integration/sequence changes
  should be validated through the Docker workflow.

## Git and PR conventions

- Branch from `main` using a descriptive `feature/...` or fix branch.
- CI runs on pushes to `main`/`develop` and on PRs to `main`.
- Keep changes focused; avoid unrelated formatting or refactors.
- Run Ruff checks and, where feasible, the Docker test workflow before opening a PR.
- Follow `.github/pull_request_template.md`: issue/description, description of changes,
  passing checks, Ruff run, branch up to date.

## Gotchas / known issues

These are current repo facts, not necessarily desired end states. Confirm intent before
"fixing" them.

1. **Integration tests hardcode `/app` paths.** `tests/integration/build_test_db.sh`
   writes `/app/tests/uorfs.db`; `test_sequence_identity.py` reads `/app/tests/...`.
   Local runs outside Docker need path edits or an equivalent mount.
2. **`logger.warnings` typo** in `uorf_utils.py` (around the check for unmapped UTR
   sequences). `logging.Logger` has `warning`, not `warnings`; this path raises
   `AttributeError` when sequences are dropped.
3. **`score_kozak` is a stub.** It always returns `1` regardless of codon and is not
   called by the build pipeline. Do not treat it as functional scoring.
4. **Schema type quirks.** `uorfs.uorf_length`, `rel_start_pos`, and `rel_stop_pos` are
   declared `TEXT` in `schema.json` even though the pipeline computes integers. Values
   round-trip as text under SQLite/validation. Changing types affects docs, JSON output,
   and downstream consumers.
5. **Docs drift.** `docs/schema.md` says the database has "five tables" but `schema.json`
   defines six (`transcripts`, `uorfs`, `utr`, `cds`, `exons`, `build_data`). It also
   contains a stray `## icon: lucide/folder-kanban` heading. `docs/schema.md` should be
   reconciled with `schema.json` when schema work happens.
6. **Version mismatch.** `pyproject.toml` and `__init__.py` declare `0.1.0`; `cli.py`
   hardcodes `%(prog)s 0.1.0`; `CITATION.cff` declares version `1.0.0`. Update all of
   these together when releasing.
7. **`src/surfa/commands/` has no `__init__.py`.** Imports work with source/editable
   installs, but verify wheel packaging if you add command modules or build a
   distribution.
8. **`src/surfa/__main__.py` is empty.** `python -m surfa` does nothing; use the `surfa`
   console script or call `surfa.cli.main` directly.
9. **`--output-dir` must already exist**, and is used for `uorfs.db` and the temporary
   `tmp.bed`. Concurrent builds sharing an output directory will collide on `tmp.bed`.
10. **Gene metadata parsing is positional and Ensembl-specific.** `gene_id` is taken
    from the second space-separated token and `gene_name` from a fixed attribute
    position. GTFs with different attribute ordering may produce wrong metadata.
11. **`query_uorf_db` swallows SQLite errors** (prints and returns `[]`), and
    `assemble_json_from_transcript` assumes exactly one transcript row and at least one
    uORF row. Missing transcripts/uORFs can surface as unrelated downstream exceptions.
12. **Test discovery pattern is `*test*.py`**, which includes both unit and integration
    tests. Running discovery at the repo root locally will pick up the `/app`-dependent
    integration tests.

## Verification checklist for changes

- `uv run ruff check .` (or at minimum, the files you touched).
- Unit tests: `uv run python -m unittest discover -s tests/unit -p '*test*.py' -v`.
- If touching GTF parsing, strand/sequence handling, codon logic, or DB schema:
  run the Docker build/test workflow and, if possible, compare a rebuilt DB against a
  known-good one for unchanged inputs.
- If touching docs or public functions: `uv run zensical build --clean`.
