import sqlite3
import numpy as np
from collections import defaultdict
import pysam
import os
import logging

logger = logging.getLogger(__name__)

__all__ = [""]

vcf = pysam.VariantFile("variants.vcf.gz")  # requires .tbi index
con = sqlite3.connect("regions.db")

for chrom, start, end, region_id in con.execute(
    "SELECT chrom, start, end, id FROM regions"
):
    # pysam fetch is 0-based half-open, same convention as BED
    for rec in vcf.fetch(chrom, start, end):
        print(region_id, rec.chrom, rec.pos, rec.ref, rec.alts)


def ensure_tabix_vcf(vcf_path: str) -> str:
    """Given a path to a VCF (plain text, regular .gz, or already bgzipped),
    return a path to a bgzipped, tabix-indexed VCF, creating one if needed.

    :param vcf_path: Path to VCF file
    :type vcf_path: str
    :return: _description_
    :rtype: str
    """

    # already bgzipped + indexed
    if vcf_path.endswith(".gz") and os.path.exists(vcf_path + ".tbi"):
        return vcf_path

    # determine output path
    if vcf_path.endswith(".vcf.gz"):
        # It's gzipped but not necessarily *bgzipped* (regular gzip won't work with tabix)
        # Safest: decompress and recompress with bgzip to guarantee block-gzip format
        plain_path = vcf_path[:-3]
        if not os.path.exists(plain_path):
            import gzip, shutil

            with gzip.open(vcf_path, "rb") as f_in, open(plain_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        src = plain_path
        out_path = vcf_path if vcf_path.endswith(".bgz") else plain_path + ".gz"
    else:
        src = vcf_path
        out_path = vcf_path + ".gz"

    # bgzip-compress
    pysam.tabix_compress(src, out_path, force=True)

    # tabix-index (preset='vcf' handles VCF's chrom/pos columns automatically)
    pysam.tabix_index(out_path, preset="vcf", force=True)

    return out_path  # creates out_path + ".tbi" alongside it


vcf_path = ensure_tabix_vcf("variants.vcf")
vcf = pysam.VariantFile(vcf_path)


# main function
def main():
    vcf_path = "variants.vcf"
    db_path = "regions.db"

    # ensure VCF is tabix-indexed
    vcf_path_ensured = ensure_tabix_vcf(vcf_path)
    vcf = pysam.VariantFile(vcf_path_ensured)

    con = sqlite3.connect(db_path)


for chrom, start, end, region_id in con.execute(
    "SELECT chrom, start, end, id FROM regions"
):
    # pysam fetch is 0-based half-open, same convention as BED
    for rec in vcf.fetch(chrom, start, end):
        print(region_id, rec.chrom, rec.pos, rec.ref, rec.alts)
