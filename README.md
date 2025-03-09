# d3-uORF-Viewer

## Tutorials

- I am following the d3 getting started guide [here](https://d3js.org/getting-started).
- It is recommending I use observable to get started, so I have begun [a project](https://observablehq.com/projects/@brad-b) using their free membership tier.

## Dev Log
1. I initialized the project using default Observable params.
2. I need to write a Python backend to perform the visualization.
    - At this time, I'll plan to start from GTF + FASTA using a dev Docker
    - Given an input gene, identify all uORF regions and their end points
    - Compute stats for each region (length, codon frame state, Kozak context)
    - Add optional filter for canonical versus non canonical start codons
3. Build Docker with required dependencies.
    - Currently using conda env named `d3-project`
    - Most appropriate Docker at this time is `bedtools.`

## Docker Usage
- Currently running a container `d3-dev` which has the dev directory mounted.

## Design Decisions
1. We will not provide initial support for MANE data, but will add later.
2. Users can specify source, one of 'havana', 'ensembl', 'ensembl_havana', 'insdc' but Havana is recommended because it uses manual review.
3. The definitions for uORF start/end sites are messy in Ensembl - there are multiple start/end definitions for a single 5'UTR in both Havana and Ensembl sources (part of this is caused by multiple exon 5'UTR regions). It may be best to use [MANE data](https://www.ncbi.nlm.nih.gov/refseq/MANE/) because it is intended as a single representative transcript for use in browsers.
    - MANE is only for human transcripts
    - We can add support later for nonhuman transcripts, using a compromise like "take the widest region when multiple definitions are encountered."
    - MANE FTP is [here](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/).
    - Can also use GENCODE Primary transcripts, which (like MANE) have stricter transcript definitions.
4. Currently we are only doing refseq but it should be straightforward to add support for Refseq, just import a different GTF file.