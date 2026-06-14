import logging
import pandas as pd
import numpy as np
import os
import subprocess

__all__ = [
    "complement_function",
    "get_transcript_FASTA",
    "produce_seqid_dict",
    "fasta_from_stdout",
    "get_seq",
    "make_bed",
    "gtf_to_sequence",
    "generate_bed_failure_report",
    "get_contigs",
    "parse_contigs_from_fasta",
]

logger = logging.getLogger(__name__)


def complement_function(input_FASTA: str) -> str:
    """This function translates negative strand nucleotides into their complements,
    but does not reverse the reading frame - must do this manually

    :param input_FASTA: FASTA nucleotide sequence
    :type input_FASTA: str
    :return: Output nucleotide string converted to reverse compliment seq
    :rtype: str
    """
    nucleotide_dict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    input_FASTA = [nucleotide_dict[k.upper()] for k in input_FASTA]

    new_codon = "".join(input_FASTA)

    return new_codon  # output new codons


def get_transcript_FASTA(ensg_df):
    """Convert FASTA sequence to reverse complement for negative strands.
    Many to one result where individual exons are concatenated together into transcript sequences.

    :param ensg_df: Ensembl GTF converted to pandas dataframe
    :type ensg_df: pandas.DataFrame
    :raises Exception: Strand Parsing Exception
    :return: Strand corrected ensg-df, where FASTA seq has been converted to reverse comp for negative strands
    :rtype: pandas.DataFrame
    """

    logger.info("Retrieving transcript FASTA.")

    if len(ensg_df.groupby("transcript").strand.nunique().unique()) > 1:
        # there was an error inferring strand identity
        raise Exception("Multiple strand values associated with input transcripts!")

    # concatenate FASTA sequences
    ensg_df["transcript_FASTA"] = ensg_df.groupby("transcript")["FASTA"].transform(
        lambda x: "".join(x)
    )
    ensg_df = (
        ensg_df[
            [
                "transcript",
                "strand",
                "chrom",
                "gene_id",
                "transcript_FASTA",
                "gene_name",
            ]
        ]
        .drop_duplicates()
        .copy()
    )

    # subset based on strand
    pos_df = ensg_df.loc[ensg_df.strand == "+"].copy()
    neg_df = ensg_df.loc[ensg_df.strand == "-"].copy()

    # apply negative complement for negative strands
    neg_df["transcript_FASTA"] = neg_df.transcript_FASTA.str[::-1].apply(
        complement_function
    )

    # append sequences together
    ensg_df = pd.concat([pos_df, neg_df]).sort_index()

    return ensg_df


def produce_seqid_dict(seqid_map, key, value):
    """Generate a seqid dict from input .tsv of Refseq/Genbank identifiers for human chromosomes.
    Use with seqid_dict in get_seq.

    :param seqid_map: Dataframe containing Refseq/Genbank identifier mapping
    :type seqid_map: pandas.Dataframe
    :param key: Column to use as key in the seqid mapping.
    :type key: str
    :param value: Column to use as value in the seqid mapping. Should match FASTA seqid.
    :type value: str
    """

    logger.info("Generating a map of contig identifiers.")
    return seqid_map.set_index(key)[value].to_dict()


def fasta_from_stdout(fasta):
    """Unpack bedtools stdout to FASTA

    :param fasta: subprocess stdout stream
    :type fasta: str
    :yield: (index, sequence) pairs unpacked from stdout
    :rtype: tuple
    """

    current_index = None
    current_sequence = ""

    for line in fasta.split("\n"):
        if line.startswith(">"):
            # Yield the previous complete sequence if it exists
            if current_index is not None:
                yield (current_index, current_sequence)

            # Start a new sequence
            current_index = line.split("::")[0][1:]
            current_sequence = ""  # Reset for new sequence
        elif line.strip():  # Skip empty lines
            current_sequence += line  # This builds up the current sequence only

    # Yield the final sequence
    if current_index is not None:
        yield (current_index, current_sequence)


def get_seq(BED_df, FASTA_path, working_dir, seqid_dict=None):
    """Get FASTA sequence from a BED file using bedtools getfasta

    :param BED_df: Pandas dataframe with named chrom, chromStart, chromEnd columns
    :type BED_df: pandas.Dataframe
    :param FASTA_path: Path to FASTA file
    :type FASTA_path: str
    :param working_dir: Working dir for temporary files.
    :type working_dir: str
    :param seqid_dict: Dictionary to remap BED_df seq IDs to match input FASTA, defaults to None
    :type seqid_dict: dict, optional
    :return: Input BED file with 1) unmapped chrom values dropped and 2) a new FASTA sequence column
    :rtype: Pandas Dataframe
    """

    logger.info("Retrieving FASTA sequence using input BED file.")

    if seqid_dict:
        # raise warning if there are empty chrom values
        empty_chrom_fields = BED_df.chrom.isna().sum()
        if empty_chrom_fields > 0:
            logger.warning(
                f"GTF file used for FASTA retrieval includes {empty_chrom_fields} empty seqids."
            )

        # remap identifiers
        original_chroms = BED_df.chrom.copy()
        BED_df.loc[:, "chrom"] = BED_df.chrom.map(seqid_dict)
        dropped_chroms = original_chroms.loc[BED_df.chrom.isna()].drop_duplicates()
        if not dropped_chroms.empty:
            logger.warning(
                "Some contigs could not be mapped to their FASTA equivalents using the input files."
            )
            for i in dropped_chroms:
                logger.warning(
                    f"Contig {i} could not be mapped between GTF and BED file and will not be included in the output database."
                )
            logger.warning(
                "Please use the seqid-map / seqid-key / seqid-value args to map between the GTF seqid and your FASTA contigs."
            )

        # drop unmapped cols
        BED_df = BED_df.loc[BED_df.chrom.notna()]

    # add identifier to bedfile
    BED_df["index"] = BED_df.index.astype(str).copy()

    # write temporary bedfile
    tmp_bed = os.path.join(working_dir, "tmp.bed")
    logger.info(f"Writing temporary BED file to {tmp_bed}.")
    BED_df.to_csv(tmp_bed, sep="\t", index=False, header=None)

    # format bedtools call
    cmd = ["bedtools", "getfasta", "-name", "-fi", FASTA_path, "-bed", tmp_bed]
    logger.info(f"Running subprocess call {''.join(cmd)}.")

    # call bedtools, capture stdout
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    # get seq from stdout
    FASTA_list = pd.DataFrame(
        list(fasta_from_stdout(result.stdout)), columns=["index", "FASTA"]
    )
    FASTA_list.loc[:, "index"] = FASTA_list["index"].astype(str)

    # merge FASTA into original BED df
    logger.info("Merging retrieved FASTA sequence with input dataframe.")
    BED_df = BED_df.merge(FASTA_list, on="index", how="outer")

    # remove tmp file
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
        logger.info(f"Removing temporary BED file {tmp_bed}.")

    return BED_df


def make_bed(ensg_df):
    logger.info("Converting input dataframe to BED format.")

    # make BED-compatable dataframe
    BED_df = (
        ensg_df[["seqname", "start", "end"]]
        .rename(columns={"seqname": "chrom", "start": "chromStart", "end": "chromEnd"})
        .copy()
    )
    BED_df.loc[:, "chromStart"] = BED_df.chromStart.astype(int) - 1
    BED_df.loc[:, "chromStart"] = BED_df.chromStart.astype(int).astype(str)
    BED_df.loc[:, "chromEnd"] = BED_df.chromEnd.astype(int).astype(str)

    return BED_df


def gtf_to_sequence(
    input_df, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value
):
    if input_df.empty:
        logger.error(
            "Could not retrieve FASTA sequence. Dataframe used for gtf_to_sequence call is empty!"
        )
        raise Exception("Dataframe used for gtf_to_sequence call is empty!")

    logger.info("Retrieving FASTA sequence for GTF input.")

    # import BED file
    BED_df = make_bed(input_df)

    if seqid_path:
        # import seqid map
        logger.info(f"Importing contig mappings from {seqid_path}.")
        seqid_map = pd.read_csv(seqid_path)

        # convert to dict
        logger.info(
            "Remapping contig identifiers in the input BED file using seqid key/values."
        )
        seqid_dict = produce_seqid_dict(seqid_map, seqid_key, seqid_value)
        BED_df = get_seq(BED_df, FASTA_path, output_dir, seqid_dict=seqid_dict)

    else:
        # map FASTA return to input BED df
        BED_df = get_seq(BED_df, FASTA_path, output_dir)

    # join FASTA back to input data
    logger.info("Joining FASTA sequence back to input GTF.")
    BED_df.loc[:, "index"] = BED_df["index"].astype(int)
    BED_df = BED_df.set_index("index").FASTA

    # log the percentage of empty FASTA sequences
    percent_empty = round((BED_df.isna().sum() / len(BED_df) * 100), 2)
    logger.debug(
        f"{BED_df.isna().sum()} / {len(BED_df)} ({percent_empty}%) BED entries returned an empty FASTA."
    )

    # check for empty FASTA
    if BED_df.isna().all():
        logger.error(
            "No sequences were retrieved from the FASTA input! Is your GTF correctly formatted?"
        )
        generate_bed_failure_report(input_df.seqname, FASTA_path, seqid_dict)

        raise Exception(
            "No sequences were retrieved from the FASTA input! Is your GTF correctly formatted?"
        )

    # create return df
    output_df = input_df.merge(BED_df, left_index=True, right_index=True)

    return output_df


def generate_bed_failure_report(chrom_array, FASTA_path, seqid_dict):
    logger.error("Generating debugging report for FASTA chromosomes.")

    # write chrom_array
    print_chroms = ", ".join([str(i) for i in np.unique(chrom_array.astype(str))])
    logger.error(
        f"Identified the following chromosome names in the input GTF: {print_chroms}."
    )

    # get contigs from FASTA file
    fasta_contigs = parse_contigs_from_fasta(FASTA_path)
    print_contigs = ", ".join([str(i) for i in fasta_contigs])
    logger.error(
        f"Identified the following chromosome names in file {FASTA_path}: {print_contigs}."
    )

    # get keys
    if seqid_dict:
        print_keys = ", ".join([str(i) for i in seqid_dict.keys()])
        logger.error(
            f"The following keys were used to remap the keys in the input GTF file: {print_keys}."
        )

        # get values
        print_values = ", ".join([str(i) for i in seqid_dict.values()])
        logger.error(
            f"This database build attempted to map the GTF keys to the following values: {print_values}."
        )


def get_contigs(FASTA_path):
    """Returns a generator - processes one line at a time"""
    with open(FASTA_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                yield line.rstrip("\n")


def parse_contigs_from_fasta(FASTA_path):
    # Usage:
    contigs = []
    for contig in get_contigs(FASTA_path):
        contigs.append(contig)

    fasta_contigs = [i.split(" ")[0].strip(">") for i in contigs]

    return fasta_contigs
