---
icon: lucide/folder-kanban
---

# Database Schema

The output database (`uorfs.db`) contains four tables.

## `transcripts`

One row per transcript parsed from the input GTF.

| Column | Type | Description |
|---|---|---|
| `transcript_id` | TEXT | Ensembl transcript ID (e.g. ENST00000000001) |
| `gene_id` | TEXT | Parent gene ID |
| `chromosome` | TEXT | Chromosome name |
| ... | | |

## `utr`

5' UTR regions associated with each transcript.

| Column | Type | Description |
|---|---|---|
| ... | | |

## `uorfs`

Upstream open reading frames identified within each UTR.

...

## `cds`

First CDS feature per transcript, used as the reference for uORF classification.

...