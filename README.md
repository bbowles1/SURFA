# d3-uORF-Viewer

Goal: Create a sqlite database to support rendering of upstream open reading frames (uORFs).

# Environment
A UV environment is configured. All packages are pip-installed.

# Docker

## Local Testing

Mini-db build:
```
make_uorf_db.py --gtf "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/mini.gtf.gz" \
    --fasta  '/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/minifasta.fa' \
    --output-dir "/Users/bbowles/Documents/Code/tmp" \
    --ensembl-source "ensembl_havana"
```

Full db build:
```
make_uorf_db.py --gtf "/Users/bbowles/Documents/Code/refdata/ensembl/Homo_sapiens.GRCh38.115.gtf.gz" \
    --fasta  '/Users/bbowles/Documents/Code/refdata/FASTA/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa' \
    --output-dir "/Users/bbowles/Documents/Code/tmp" \
    --ensembl-source "ensembl_havana" 
```

JSON export:
```
make_json.py \
    --db /Users/bbowles/Documents/Code/GitHub/Upstream-Display/data/uorfs.db \
    --transcript ENST00000504921.7 \
    --output /Users/bbowles/Documents/Code/GitHub/Upstream-Display/data/uorfs.json
```