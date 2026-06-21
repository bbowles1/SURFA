---
icon: lucide/flask-conical
---

# SURF-A

SURF-A is a Python utility that allows users to annotate small upstream open reading frames (uORFs) in mRNA sequences. The tool provides a way to build a basic database of uORF sequences from any input set of Ensembl GTF and FASTA files.

## Installation

This work is in development and full PyPI hosting is coming soon. This package depends on a locally available version of Bedtools.

For development, you can install this tool locally using the following steps:

1. Install bedtools (https://bedtools.readthedocs.io/en/stable/).
2. Install uv: (https://docs.astral.sh/uv/).
3. `uv sync` to install all required dependencies.
4. `source .venv/bin/activate` to activate the environment.

## Inputs

The required inputs for the SURF-A database build are:

- An Ensembl GTF file for your organism of choice.
- A matching FASTA file (should use the same reference genome as the GTF).

SURF-A wraps efficient queries using Bedtools and uses the resulting sequences to call uORF regions.


## Commands

#### Build 

```
surfa build --gtf "/Users/bbowles/Documents/Code/refdata/ensembl/Homo_sapiens.GRCh38.115.gtf.gz" \
    --fasta  '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa' \
    --output-dir "/Users/bbowles/Documents/Code/GitHub/Upstream-Display/data/" \
    --ensembl-source "ensembl_havana" \
    --log-level DEBUG
```


Query









## Usage

```
surfa <command> [options]
```

### Global Options

| Option | Description |
|--------|-------------|
| -V, --version	| Show program version and exit |
| -l, --log-level LEVEL| Set logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO) |
| --log-file FILE | Write logs to specified file (default: surfa.log) |

### Commands

#### build

Build a SQLite database of uORF calls from an input FASTA and GTF sequence.

Note: Specific build options are available via `surfa build --help.`

```
surfa build --gtf <gtf-file> --fasta <fasta-file> --output-dir <directory> [options]
```

**Required Arguments**

| Argument | Description |
|----------|-------------|
| --gtf <FILE>	Path to Ensembl-format GTF file
| --fasta <FILE>	Path to input FASTA file
| --output-dir <DIR>	Output directory for generated files

**Optional Arguments**

| Argument | Description |
|----------|-------------|
| --ensembl-source <SOURCE>	| Which Ensembl GTF data source to use (default: ensembl_havana). |
| --seqid-map <FILE>	| Dictionary file mapping GenBank/RefSeq identifiers |
| --seqid-key <COLUMN>	| Column key in seqid_map for GTF chromosome identifiers |
| --seqid-value <COLUMN> | Column in seqid_map containing remapped values |

**Example Build**

```
surfa build \
    --gtf genes.gtf \
    --fasta sequences.fasta \
    --output-dir ./uorf_database \
    --ensembl-source ensembl
```


#### query

Query the uORF database for specific annotations.

Note: Specific query options are available via `surfa query --help.`

Examples

Basic uORF Database Construction

| Argument | Description |
|----------|-------------|
| --db <FILE> | Path to uorfs.db (created using surfa build command). |
| --transcript <STR> | Target transcript (including version number). |
| --output | Output file name (incl JSON extension). Defaults to `query.json` |
