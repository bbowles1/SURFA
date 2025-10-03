#!/bin/bash
python /app/python/gtf_to_json.py --gtf "/app/tests/mini.gtf.gz" \
    --fasta  '/app/tests/minifasta.fa' \
    --output-dir "/app/tests/" \
    --ensembl-source "ensembl_havana"