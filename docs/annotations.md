---
icon: lucide/folder-kanban
---

# Annotation Guide

This package makes some high-level decisions on how to handle the annotation of UTR regions.

- A noncoding exon is an exon which contains part of the 5'UTR sequence.
- We first compile a FASTA seq by aggregating all noncoding exons together.
- We build a sequence of the first N noncoding exons of a transcript, which usually contain some CDS at the end.
- Samples without a mapped CDS start site are removed from the database.
- Focus on 5'UTR exclusively
- N-terminal extensions?
- How uORF IDs are generated?
- How are database entries handled? Are there only 

# Abbreviations
- CDS
- UTR

