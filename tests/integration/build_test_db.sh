#!/bin/bash
surfa build --gtf "/app/tests/data/mini.gtf.gz" \
    --fasta  '/app/tests/data/minifasta.fa' \
    --output-dir "/app/tests/" \
    --ensembl-source "ensembl_havana"