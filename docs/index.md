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


## Usage

The surfa cli tool has several options such as `build` and `query` which can be called using:

```
surfa <command> [options]
```

Each command has the following set of global options:

| Option | Description |
|--------|-------------|
| -V, --version	| Show program version and exit |
| -l, --log-level LEVEL| Set logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO) |
| --log-file FILE | Write logs to specified file (default: surfa.log) |


## surfa build

`surfa build` creates a SQLite database of uORF calls from an input FASTA and GTF sequence. Specific build options are available via `surfa build --help.`

```
surfa build \
  --gtf <gtf-file> \
  --fasta <fasta-file> \
  --output-dir <directory> \
  [options]
```

**Required Arguments**

| Argument | Description |
|----------|-------------|
| --gtf <FILE> | Path to Ensembl-format GTF file |
| --fasta <FILE> | Path to input FASTA file |
| --output-dir <DIR> | Output directory for generated files |

**Optional Arguments**

| Argument | Description |
|----------|-------------|
| --ensembl-source <SOURCE>	| Which Ensembl GTF data source to use (default: ensembl_havana). |
| --seqid-map <FILE>	| Dictionary file mapping GenBank/RefSeq identifiers |
| --seqid-key <COLUMN>	| Column key in seqid_map for GTF chromosome identifiers |
| --seqid-value <COLUMN> | Column in seqid_map containing remapped values |


**Seqid Mapping**

Occasionally, seqids in the Ensembl GTF will not match the names in your FASTA input file, as can happen when you switch between RefSeq and GenBank identifiers. You can provide a .csv file mapping between these values using the `--seqid-map` argument, where `--seqid-key` is the set of GTF seqid values while `--seqid-value` is the set of corresponding FASTA identifiers. An example input csv is below:


| molecule_name | number | chr_abbreviation | chrom_abbreviation | genbank_sequence | refseq_sequence |
|---------------|--------|------------------|--------------------|------------------|-----------------|
| Chromosome 1 | 1 | chr1 | chrom1 | CM000663.1 | NC_000001.10 |
| Chromosome 2 | 2 | chr2 | chrom2 | CM000664.1 | NC_000002.11 |

Your seqid-map.csv can contain an arbitrary set of columns so long as you refer to the key-value pairs using your `--seqid-key` and `--seqid-value` arguments.


**Example Build**

```
surfa build \
    --gtf genes.gtf \
    --fasta sequences.fasta \
    --output-dir ./uorf_database \
    --ensembl-source ensembl
```

**Example Build Using Seqid Mappings**

```
surfa build \ 
  --gtf Homo_sapiens.GRCh38.115.gtf.gz \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output-dir /data \
  --ensembl-source ensembl_havana \
  --seqid-map seqid_map.csv \
  --seqid-key "chr_abbreviation" \
  --seqid-value "number"
```


## surfa query

Query a pre-built uORF database for specific annotations. Specific query options are available via `surfa query --help.`

**Arguments**

| Argument | Description |
|----------|-------------|
| --db <FILE> | Path to uorfs.db (created using surfa build command). |
| --transcript <STR> | Target transcript (including version number). |
| --output | Output file name (incl JSON extension). Defaults to `query.json` |


**Example Query**

```
surfa query \
    --db /data/uorfs.db \
    --transcript ENST00000504921.7 \
    --output /data/uorfs.json
```