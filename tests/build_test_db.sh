#!/bin/bash
python /app/make_uorf_db.py --gtf "/app/tests/mini.gtf.gz" \
    --fasta  '/app/tests/minifasta.fa' \
    --output-dir "/app/tests/" \
    --ensembl-source "ensembl_havana"