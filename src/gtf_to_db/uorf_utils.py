import logging
import pandas as pd
import numpy as np
import os
import math
import sqlite3
import re
from gtf_to_db.fasta_utils import gtf_to_sequence, get_transcript_FASTA
from gtf_to_db.db_utils import write_to_db


__all__ = [
    "chunker",
    "score_kozak",
    "get_codons",
    "fasta_codon_search",
    "filter_codons",
    "match_codons",
    "return_FASTA",
    "return_FASTA_optimized",
    "get_uorfs",
    "import_reference",
    "get_exon_field_num",
    "get_transcript_field_num",
    "get_transcript_version_field_num",
    "check_identity",
    "unpack_transcript",
    "gtf_to_uorf_db",
]


def chunker(seq, size):
    # chunk data into n sizes
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def score_kozak(codon):
    """
    Score binding strength of each start codon. For now, returns ordinal values
    but is a place holder for more advanced functionality in the future.
    """

    codon_dict = {
        "AUG": 0,
        "CUG": 1,
        "GUG": 2,
        "UUG": 3,
        "ACG": 4,
        "AAG": 5,
        "AGG": 6,
        "AUC": 7,
        "AUU": 8,
        "AUA": 9,
    }

    return 1


def get_codons(seq, frame):
    """Split nucleotide sequence into codon blocks for a given frame

    :param seq: Nucleotide sequence to split
    :type seq: str
    :param frame: Frame (0,1,2) to split on.
    :type frame: int
    :return: List of split codons
    :rtype: list
    """

    # split sequence into codons for a given frame
    rna_seq = seq[frame:]
    return [rna_seq[i : i + 3] for i in range(0, len(rna_seq) - 2, 3)]


def fasta_codon_search(RNA, frame):
    """Split nucleotide sequence into codon blocks for a given frame.
    Faster version of get_codons.

    :param RNA: RNA nucleotide sequence
    :type RNA: str
    :param frame: Frame (0,1,2) to split on.
    :type frame: int
    :return: list of RNA nucleotides, split by frame.
    :rtype: list
    """

    # Begin RNA reading frame at Nth position, then iterate over in chunks of 3
    return list(map("".join, zip(*[iter(RNA[frame:])] * 3)))


def filter_codons(codons, targets):
    """Filter codons to a target list.

    :param codons: List of nucleotide codons.
    :type codons: list
    :param targets: List of target codons (ie AUG, CUG).
    :type targets: list
    :return: Filtered set of codons.
    :rtype: list
    """

    target_set = set(targets)  # O(1) lookups
    return [(codon, idx) for idx, codon in enumerate(codons) if codon in target_set]


def match_codons(start_codons, stop_codons):
    """Find next downstream stop codon for each start codon using numpy.

    :param start_codons: List of tuples with format (start_codon, position)
    :type start_codons: List
    :param stop_codons: List of tuples with format (stop_codon, position)
    :type stop_codons: list
    :return: List of codons pairs
    :rtype: list
    """

    if not stop_codons:
        return [(start, None) for start in start_codons]
    if not start_codons:
        return []

    # extract indices from input codons
    start_indices = np.array([idx for _, idx in start_codons])
    stop_indices = np.array([idx for _, idx in stop_codons])

    # use searchsorted to find next stop codon for each start
    next_stop_positions = np.searchsorted(stop_indices, start_indices, side="right")

    result = []
    for i, start in enumerate(start_codons):
        stop_number = next_stop_positions[i]
        if stop_number < len(stop_codons):
            result.append((start, stop_codons[stop_number]))
        else:
            result.append((start, None))  # no downstream stop codon

    return result


def return_FASTA(FASTA, codon_tuple):
    """Subset an input FASTA sequence to the supplied codon boundaries.

    :param FASTA: FASTA sequence of the UTR
    :type FASTA: str
    :param codon_tuple: tuple with format ((start_codon, position), (stop_codon, position))
    :type codon_tuple: tuple
    :return: uORF FASTA sequence
    :rtype: str
    """

    # check if stop is None
    if codon_tuple[1] is None:
        return FASTA[codon_tuple[0][1] :]
    else:
        return FASTA[codon_tuple[0][1] : codon_tuple[1][1] + 3]


def return_FASTA_optimized(FASTA, codon_tuple):
    """Optimized approach to return FASTA sequence from start and stop tuples

    :param FASTA: FASTA sequence of the UTR
    :type FASTA: str
    :param codon_tuple: tuple with format ((start_codon, position), (stop_codon, position))
    :type codon_tuple: tuple
    :return: uORF FASTA sequence
    :rtype: str
    """

    start = codon_tuple[0][1]
    end = codon_tuple[1]

    # Direct slice assignment is faster than conditional
    return FASTA[start:] if end is None else FASTA[start : end[1] + 3]


def get_uorfs(input_df):
    """
    Analyze a table of transcript information and return a table of uORFs

    :param FASTA_df: Table containing `transcript_FASTA` olumn
    :type FASTA_df: pandas Data.Frame
    :return: Pandas dataframe containing uORF information
    :rtype: pandas Data.Frame
    """

    logger = logging.getLogger(__name__)
    logger.info("Retrieving uORFs from input dataframe of FASTA sequence.")

    # convert input to upper
    FASTA_df = input_df.copy()
    FASTA_df["transcript_FASTA"] = FASTA_df.transcript_FASTA.str.upper()

    # set start and stop codons
    stop_codons = {"UAA", "UAG", "UGA"}
    start_codons = {
        "AUG",  # Canonical
        "CUG",  # near-cognate
        "GUG",
        "UUG",
        "ACG",
        "AAG",
        "AGG",
        "AUC",
        "AUU",
        "AUA",  # rare near-cognate
    }

    all_codons = stop_codons.union(start_codons)

    # split data into three frame states

    n_rows = FASTA_df.shape[0]
    FASTA_df = FASTA_df.loc[FASTA_df.index.repeat(3)].reset_index(
        drop=True
    )  # efficient way to copy dataframe 3x
    FASTA_df["frame"] = np.tile([0, 1, 2], n_rows)

    #########
    # STEPS #
    #########

    # find uORFs in frame:
    # 1. Convert to RNA seq
    # 1. get codons based on frame
    # 2. enumerate over to get start and indices

    # convert FASTA to RNA
    FASTA_df.loc[:, "transcript_FASTA"] = FASTA_df.transcript_FASTA.str.replace(
        "T", "U"
    )

    # apply fast codon search
    # faster than vectorized approaches for a chunk size of 1k
    FASTA_df["codons"] = FASTA_df.apply(
        lambda row: fasta_codon_search(row.transcript_FASTA, row.frame), axis=1
    )

    # retrieve start/stop codons from the list
    FASTA_df["start_codons"] = FASTA_df.codons.apply(
        lambda x: filter_codons(x, start_codons)
    )
    FASTA_df["stop_codons"] = FASTA_df.codons.apply(
        lambda x: filter_codons(x, stop_codons)
    )

    # determine which codon pairs are in frame with each other
    FASTA_df["codon_pairs"] = FASTA_df.apply(
        lambda row: match_codons(row.start_codons, row.stop_codons), axis=1
    )
    FASTA_df.drop(columns=["start_codons", "stop_codons"], inplace=True)

    # drop rows with no codon pairs
    FASTA_df = FASTA_df.loc[FASTA_df["codon_pairs"].apply(bool)]

    # explode rows --> 1 uORF per row
    uorf_df = FASTA_df.explode("codon_pairs")

    # unpack codon tuples -> need end, start POS
    uorf_df["start_codon"] = uorf_df.codon_pairs.str[0].str[0]
    uorf_df["start_codon_pos"] = uorf_df.codon_pairs.str[0].str[1]
    uorf_df["stop_codon"] = uorf_df.codon_pairs.str[1].str[0]
    uorf_df["stop_codon_pos"] = uorf_df.codon_pairs.str[1].str[1]

    # recode missing stop codons - None entries indicate uORF reads into main CDS.
    logger.debug("Recoding missing stop codons as NO_UTR_STOP.")
    uorf_df.loc[uorf_df.stop_codon.isna(), "stop_codon"] = "NO_UTR_STOP"

    # get uORF sequence from transcript FASTA
    # entries where stop codon is missing are read from the uORF to the CDS start
    sequences = []
    adjusted_bp_start = []
    adjusted_bp_stop = []
    for fasta, start, end in zip(
        uorf_df["transcript_FASTA"],
        uorf_df["start_codon_pos"],
        uorf_df["stop_codon_pos"],
    ):
        bp_start = start * 3
        if pd.isna(end):
            bp_end = np.nan
            sequences.append(fasta[bp_start:])
        else:
            bp_end = (end + 1) * 3
            sequences.append(fasta[bp_start : int(bp_end)])
        adjusted_bp_start.append(bp_start)
        adjusted_bp_stop.append(bp_end)

    # append values back to uorf df
    uorf_df["uORF_FASTA"] = sequences
    uorf_df["start_bp_pos"] = adjusted_bp_start
    uorf_df["stop_bp_pos"] = adjusted_bp_stop

    # get total uorf length
    uorf_df["uORF_length"] = uorf_df.stop_bp_pos - uorf_df.start_bp_pos

    # generate anotations
    # 1. is entry confined to uORF or does it extend into CDS? -> "uorf_location"
    # 2. is entry in-frame with upstream CDS or out of frame? -> "CDS_frame_status"
    uorf_df["uorf_location"] = np.where(
        uorf_df.stop_bp_pos.isna(), "CDS_overlap", "UTR_only"
    )

    uorf_df["CDS_frame_status"] = np.where(
        (uorf_df.utr_len - uorf_df.frame) % 3 == 0, "in-frame", "out-frame"
    )

    # drop transcript-level info: strand, FASTA, total_length (these are retrievable from the transcript table)
    # drop intermmediate info no longer needed: codon_pairs
    output_cols = [
        "transcript",
        "uORF_length",
        "start_bp_pos",
        "stop_bp_pos",
        "uORF_FASTA",
        "start_codon",
        "stop_codon",
        "frame",
        "uorf_location",
        "CDS_frame_status",
    ]
    uorf_df = uorf_df[output_cols]

    # rename columns to be more explicit
    # start/stop pos are relative to UTR start and not genomic coords
    uorf_df.rename(
        columns={"start_bp_pos": "rel_start_pos", "stop_bp_pos": "rel_stop_pos"},
        inplace=True,
    )

    return uorf_df


def import_reference(path):
    """Parse a reference FASTA file. Intended to import test data.

    :param path: Path to reference FASTA sequence.
    :type path: str
    :return: A concatenated FASTA sequence.
    :rtype: str
    """

    # convert FASTA entry to sequence
    sequence = []
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            if (">" in line) or ("#" in line):
                pass
            else:
                sequence.append(line.strip(" \n"))
    return "".join(sequence)


def get_exon_field_num(string):
    """Find the index of the exon number in a GTF file's attribute column.

    :param string: attribute string from a GTF file.
    :type string: str
    :raises Exception: Could not detect the index for the given string.
    :return: index of "exon_number" field.
    :rtype: int
    """

    # find index of exon number
    pattern = r"exon_number"
    index = next(
        (i for i, item in enumerate(string.split(";")) if re.search(pattern, item)),
        None,
    )
    if index:
        return index
    else:
        logging.error(f"No index field detected for {pattern}!")
        raise Exception(f"No index field detected for {pattern}!")


def get_transcript_field_num(string):
    """Find the index of the transcript ID in a GTF file's attribute column.

    :param string: attribute string from a GTF file.
    :type string: str
    :raises Exception: Could not detect the index for the given string.
    :return: index of "transcript_id" field.
    :rtype: int
    """

    # find index of exon number
    pattern = r"transcript_id"
    index = next(
        (i for i, item in enumerate(string.split(";")) if re.search(pattern, item)),
        None,
    )
    if index:
        return index
    else:
        logging.error(f"No index field detected for {pattern}!")
        raise Exception(f"No index field detected for {pattern}!")


def get_transcript_version_field_num(string):
    """Find the index of the transcript version in a GTF file's attribute column.

    :param string: attribute string from a GTF file.
    :type string: str
    :raises Exception: Could not detect the index for the given string.
    :return: index of "transcript_version" field.
    :rtype: int
    """

    # find index of exon number
    pattern = r"transcript_version"
    index = next(
        (i for i, item in enumerate(string.split(";")) if re.search(pattern, item)),
        None,
    )
    if index:
        return index
    else:
        logging.error(f"No index field detected for {pattern}!")
        raise Exception(f"No index field detected for {pattern}!")


def check_identity(region_start, strand, CDS_start):
    """Check whether an input GTF UTR region is 5' or 3'

    :param region_start: UTR region start position
    :type region_start: int
    :param strand: Strand (+ or -)
    :type strand: str
    :param CDS_start: Coding DNA sequence start position
    :type CDS_start: int
    :return: 5 or 3
    :rtype: int
    """

    if strand == "+":
        if region_start > CDS_start:
            return 3
        elif region_start < CDS_start:
            return 5
    if strand == "-":
        if region_start > CDS_start:
            return 5
        elif region_start < CDS_start:
            return 3


def unpack_transcript(input_df, df_name):
    """Unpack the transcript from the attribute column to a separate field.
    Includes handling for a variety of GTF formats.

    :param input_df: Pandas Dataframe.
    :type input_df: pd.DataFrame
    :param df_name: Name of dataframe for logging.
    :type df_name: str
    :return: Dataframe with transcript identity unpacked to specific column.
    :rtype: pd.DataDrame
    """

    logger = logging.getLogger(__name__)
    logger.debug(
        "Checking input dataframe for separate transcript and version definitions."
    )

    return_df = input_df.copy()

    # get sample attributes field
    sample_string = input_df.attribute.iloc[0]
    logger.debug(
        f"Sample string in input {df_name} GTF attributes field: {sample_string}"
    )

    # unpack transcript
    transcript_index = get_transcript_field_num(sample_string)
    logger.debug(f"Transcript detected in field {transcript_index}.")
    return_df["transcript"] = (
        return_df.attribute.str.split(";")
        .str[transcript_index]
        .str.split(" ")
        .str[2]
        .str.strip('"')
    )

    # determine if separate handling is required to append ENST + version number
    separate_transcript_version = "transcript_version" in input_df.attribute.iloc[0]
    logger.debug(
        f"Separate transcript version handling: {separate_transcript_version}."
    )

    if separate_transcript_version:
        version_index = get_transcript_version_field_num(sample_string)
        logger.debug(f"Transcript version number detected in field {version_index}.")
        return_df["transcript_version"] = (
            return_df.attribute.str.split(";")
            .str[version_index]
            .str.split(" ")
            .str[2]
            .str.strip('"')
        )

        logger.debug("Appending transcript number and transcript version.")
        return_df.loc[:, "transcript"] = (
            return_df.transcript + "." + return_df.transcript_version
        )

    else:
        logger.debug(
            "Transcript version is present in transcript definition. No appending required."
        )

    return return_df


def gtf_to_uorf_db(
    gtf_path,
    FASTA_path,
    output_dir,
    source,
    seqid_path=None,
    seqid_key=None,
    seqid_value=None,
):
    """Main function of the SURF-A package. Identify uORF sequences from an input GTF file, provide basic annotations,
    and supply a sqlite database of information.

    :param gtf_path: Path to Ensembl GTF file.
    :type gtf_path: str
    :param FASTA_path: Path to FASTA sequence.
    :type FASTA_path: str
    :param output_dir: Path to output directory for files and logs.
    :type output_dir: str
    :param source: Source to use within Ensembl GTF (ie Havanna)
    :type source: str
    :param seqid_path: Path to seqid dictionary which maps contigs in the FASTA to their identity in the GTF, defaults to None
    :type seqid_path: str, optional
    :param seqid_key: Key to use for mapping between contigs and GTF seqids, should match those in seqid_path file, defaults to None
    :type seqid_key: str, optional
    :param seqid_value: Value to use for mapping between contigs and GTF seqids, should match those in seqid_path file, defaults to None
    :type seqid_value: str, optional
    :raises Exception: Exception if an input path does not exist.
    """

    # setup logging
    logger = logging.getLogger(__name__)
    logger.info("Beginning database build.")

    # check that all paths exist
    for path in [gtf_path, output_dir, FASTA_path]:
        if not os.path.exists(path):
            logger.error(f"Input path {path} does not exist!")
            raise Exception(f"Path {path} does not exist!")
    logger.debug("All input paths exist.")

    ###############
    # DATA IMPORT #
    ###############

    # Annotate gene_start: the start of coding at the canonical transcript ATG. This is the CDS start position in the Ensembl .gff3.
    logger.debug(f"Importing input GTF {gtf_path}.")
    ensg_df = pd.read_csv(
        gtf_path,
        header=None,
        sep="\t",
        comment="#",
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    logger.debug(f"Input GTF imported. Dataframe has columns {ensg_df.columns}.")

    # subset to source
    logger.debug(f"Possible sources for input GTF: {ensg_df.source.unique()}.")
    ensg_df = ensg_df.loc[ensg_df.source == source]
    logger.info(f"Subsetting to dataframe source {source}.")

    ####################
    # EXTRACT CDS INFO #
    ####################

    # get length of first CDS bound
    cds = ensg_df.loc[ensg_df.feature == "CDS"].copy()

    # unpack transcript
    cds = unpack_transcript(cds, "CDS")

    # determine indices
    sample_string = cds.attribute.iloc[0]
    exon_index = get_exon_field_num(sample_string)

    # unpack fields
    cds["exon"] = (
        cds.attribute.str.split(";")
        .str[exon_index]
        .str.split(" ")
        .str[2]
        .str.strip('"')
        .astype(int)
    )
    cds["length"] = cds.end + 1 - cds.start
    cds.sort_values(by=["transcript", "exon"], ascending=True, inplace=True)

    # subset to first exon of CDS (which is split between 5' UTR and beginning of CDS)
    first_cds = cds.groupby("transcript").first().reset_index().drop_duplicates()
    first_cds = first_cds[["seqname", "transcript", "exon", "length", "start", "end"]]
    first_cds["type"] = "CDS"

    #################
    # GENE FEATURES #
    #################

    exons = ensg_df.loc[ensg_df.feature == "exon"].copy()

    # unpack transcripts
    exons = unpack_transcript(exons, "exons")

    # determine indices
    sample_string = exons.attribute.iloc[0]
    exon_index = get_exon_field_num(sample_string)

    # create exons table
    exons["length"] = exons.end + 1 - exons.start
    exons["exon"] = exons.attribute.str.split(";").str[exon_index].str.split(" ").str[2]
    exons["rel_end"] = exons.groupby("transcript")["length"].cumsum()
    exons["rel_start"] = exons["rel_end"] - exons["length"]

    ###################
    # 5' UTR FEATURES #
    ###################
    # Annotate 5'UTR features, define as all 5'UTR exons, excluding CDS sequences

    # bool: determine if we have a new or old 5'UTR annotation format
    new_utr_format = "five_prime_utr" in ensg_df.feature.values

    if new_utr_format:
        # ensembl dataframe is newer format which delimits 5'UTR
        utr_df = ensg_df.loc[ensg_df.feature == "five_prime_utr"].copy()

        # unpack transcripts
        utr_df = unpack_transcript(utr_df, "UTR")

        # determine indices
        sample_string = utr_df.attribute.iloc[0]

        # manually enumerate exon number
        utr_df.loc[utr_df.strand == "+", "exon"] = (
            utr_df.loc[utr_df.strand == "+"]
            .sort_values(by="start", ascending=True)
            .groupby("transcript")
            .cumcount()
            + 1
        )
        utr_df.loc[utr_df.strand == "-", "exon"] = (
            utr_df.loc[utr_df.strand == "-"]
            .sort_values(by="start", ascending=False)
            .groupby("transcript")
            .cumcount()
            + 1
        )

    else:
        # subset to UTRs
        utr_df = ensg_df.loc[ensg_df.feature == "UTR"].copy()

        # unpack transcripts
        utr_df = unpack_transcript(utr_df, "UTR")

        # determine indices
        sample_string = utr_df.attribute.iloc[0]
        exon_index = get_exon_field_num(sample_string)

        utr_df["exon"] = (
            utr_df.attribute.str.split(";").str[exon_index].str.split(" ").str[2]
        )

    utr_df.sort_values(by=["transcript", "exon"], ascending=True, inplace=True)

    # merge in CDS start
    utr_df = pd.merge(
        utr_df,
        first_cds[["transcript", "start"]].rename(columns={"start": "cds_start"}),
        on=["transcript"],
        how="left",
    )

    # drop rows where CDS start could not be mapped
    no_cds_rows = utr_df.loc[utr_df.cds_start.isna()].shape[0]
    if no_cds_rows > 0:
        logger.warning(
            f"Removed {no_cds_rows} rows where CDS start is not defined in the input GTF."
        )
    utr_df = utr_df.loc[utr_df.cds_start.notna()]

    # determine UTR status
    if new_utr_format:
        utr_df["utr_type"] = utr_df.apply(
            lambda row: check_identity(row.start, row.strand, row.cds_start), axis=1
        )
        utr_df = utr_df.loc[utr_df.utr_type == 5]

    #############
    # CDS FRAME #
    #############
    # now that we have UTR info, calculate CDS frame state

    utr_df["length"] = utr_df.end + 1 - utr_df.start

    frame_state = (utr_df.groupby("transcript")["length"].sum() % 3).reset_index()
    frame_state.rename(columns={"length": "frame_state"}, inplace=True)
    utr_df = pd.merge(utr_df, frame_state)

    ###################
    # UTR ANNOTATIONS #
    ###################

    # calculate relative start and stop of each UTR exon
    utr_df["rel_stop"] = utr_df.groupby("transcript").length.cumsum()
    utr_df["rel_start"] = utr_df["rel_stop"] - utr_df.length

    #########
    # FASTA #
    #########

    logger.info("Matching GTF to FASTA sequences...")

    # new method for retrieving FASTA seq
    utr_df = gtf_to_sequence(
        utr_df, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value
    )
    first_cds = gtf_to_sequence(
        first_cds, FASTA_path, output_dir, seqid_path, seqid_key, seqid_value
    )

    logger.info("FASTA sequences retrieved.")

    # drop columns where FASTA sequence could not be mapped
    drop_rows = utr_df.loc[utr_df.FASTA.isna()].index
    if not drop_rows.empty:
        logger.warnings(
            f"{len(drop_rows)} FASTA sequences could not be mapped to the input GTF."
        )
    utr_df = utr_df.loc[utr_df.FASTA.notna()]

    # convert exon FASTA to transcript FASTA
    # 1. concatenate FASTA sequences for all exons of the same transcript
    # 2. Apply reverse complement for negative strands
    utr_df["gene_id"] = utr_df.attribute.str.split(" ").str[1].str.strip('"; ,')
    utr_df["gene_name"] = (
        utr_df.attribute.str.split(";").str[3].str.split(" ").str[2].str.strip('"; ,')
    )
    transcript_df = get_transcript_FASTA(utr_df.rename(columns={"seqname": "chrom"}))
    # determine transcript length
    transcript_df["length"] = transcript_df.transcript_FASTA.str.len()

    ##############
    # UTR LENGTH #
    ##############
    # calc UTR length and append to transcript table
    utr_len = utr_df.groupby("transcript")["length"].sum().reset_index()
    transcript_df = pd.merge(
        transcript_df, utr_len.rename(columns={"length": "utr_len"})
    )

    ##############
    # CDS LENGTH #
    ##############

    # calculate rel start/stop annotations for CDS display
    first_cds = pd.merge(
        first_cds, transcript_df[["transcript", "utr_len"]].drop_duplicates()
    )
    first_cds.rename(columns={"utr_len": "rel_start"}, inplace=True)
    first_cds["rel_stop"] = first_cds.rel_start + first_cds.length

    #############
    # TRANSFORM #
    #############
    # transform data into SQLite
    # 1. group by transcript
    # 2. generate transcript exon table
    # 3. generate uORF sequence table
    # 4. Write to sqlite
    # 5. Make function to unpack sqlite -> JSON

    logger.info("Searching transcript sequences for upstream open reading frames...")

    # iterate over transcripts in chunks
    # extrat uORF sequenes
    uorf_list = []
    size = 1000
    count = 0
    n_chunks = math.ceil(transcript_df.shape[0] / size)
    for chunk in chunker(transcript_df, size):
        progress = count / n_chunks * 100
        print(f"\rSearching for uORFs: {progress:.1f}%", end="", flush=True)

        tmp_table = get_uorfs(chunk)
        uorf_list.append(tmp_table)

        count += 1
    print("\rSearching for uORFs: 100%\n", end="", flush=True)
    uorf_table = pd.concat(uorf_list)

    logger.info("Upstream open reading frames identified.")

    # sort, add a unique uORF identifier
    uorf_table = uorf_table.sort_values(
        by=["transcript", "start_codon", "rel_start_pos", "rel_stop_pos", "uORF_length"]
    ).rename(columns={"uORF_length": "uorf_length"})

    # add a unique uORF identifier: ENST + codon + number
    uorf_table["uorf_count"] = (
        uorf_table.groupby(["transcript", "start_codon"]).cumcount() + 1
    )
    uorf_table["uorf_id"] = (
        uorf_table.transcript
        + "."
        + uorf_table.start_codon
        + "."
        + uorf_table.uorf_count.astype(str)
    )
    uorf_table.drop(columns=["uorf_count"], inplace=True)

    # count number of uORFs per transcript
    orfs_per_enst = uorf_table.groupby("transcript").size().reset_index()
    orfs_per_enst.columns = ["transcript", "uorf_count"]
    transcript_df = pd.merge(transcript_df, orfs_per_enst)

    ###############
    # NA HANDLING #
    ###############

    # set handling for end position, length
    length_dict = transcript_df.set_index("transcript")["length"].to_dict()

    uorfs_with_unmapped_start = (uorf_table.stop_codon == "NO_UTR_STOP").sum()
    logger.debug(
        f"{uorfs_with_unmapped_start} uORFs with non-UTR stop codons. Handling unmapped length, end anotations."
    )

    uorf_table.loc[uorf_table.stop_codon == "NO_UTR_STOP", "rel_stop_pos"] = (
        uorf_table.loc[uorf_table.stop_codon == "NO_UTR_STOP"].transcript.map(
            length_dict
        )
    )
    uorf_table.loc[uorf_table.stop_codon == "NO_UTR_STOP", "uorf_length"] = (
        uorf_table.loc[uorf_table.stop_codon == "NO_UTR_STOP"].transcript.map(
            length_dict
        )
    )
    uorf_table.loc[:, "rel_stop_pos"] = uorf_table.rel_stop_pos.astype(int)
    uorf_table.loc[:, "uorf_length"] = uorf_table.uorf_length.astype(int)

    ##########
    # SQLite #
    ##########
    # export the following tables:
    # 1. transcript table
    # 2. UTR exon table
    # 3. uORF table

    # rename columns in utr_df
    utr_df.rename(columns={"seqname": "chrom"}, inplace=True)

    db_path = os.path.join(output_dir, "uorfs.db"
    write_to_db(
        dataframes={
            "transcripts": transcript_df,
            "utr": utr_df,
            "uorfs": uorf_table,
            "cds": first_cds,
        },
        db_path=db_path),
    )
    logger.info(f"Database saved to {db_path}.")
    print(f"Database saved to {db_path}")
