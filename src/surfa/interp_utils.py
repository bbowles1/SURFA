#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 16:41:59 2021

Sequencing tools for uORF analysis

1.3.2022 update: clarified wording for certain effects:
novel_start > novel_near_cognate_start
start_downregulation > start_Kozak_downregulation
start_upregulation > start_Kozak_upregulation
novel_CTG > novel_CTG_start
novel_ATG > novel_ATG_start

@author: m181414
"""

import pandas as pd
import numpy as np
import os


stop_codons = ["TAA", "TAG", "TGA"]
start_codons = [
    "ATG",
    "TTG",
    "GTG",
    "CTG",
    "AGG",
    "ACG",
    "ATA",
    "ATT",
    "ATC",
]  # all possible start codons, AAG is not being considered due to very low ribosome loading in literature

low_start_codons = [
    "TTG",
    "GTG",
    "AGG",
    "ACG",
    "ATA",
    "ATT",
    "ATC",
]  # low activity start codons per literature
high_start_codons = ["ATG", "CTG"]  # highly active start codons per literature


class Error(Exception):
    """Base class for other exceptions"""

    pass


class Processing_Error(Error):
    """Raised when an issue in dataframe processing has occurred"""

    pass


def FASTA_function(POS, FASTA_range, FASTA, strand):  # calculates new_FASTA for SNVs
    for interval in FASTA_range:
        if POS in interval:
            range_index = FASTA_range.index(interval)  # n of interval

            relative_POS = POS - interval.start  # relative POS within interval

            if range_index > 0:
                FASTA_range = FASTA_range[
                    :(range_index)
                ]  # cuts FASTA_range to all ranges before the range_index
                FASTA_POS = [
                    (i.stop - i.start) for i in FASTA_range
                ]  # calculates # of nucleotides in each range
                FASTA_POS = sum(
                    FASTA_POS
                )  # calculates all nucleotides in preceeding ranges

                relative_POS = (
                    relative_POS + FASTA_POS
                )  # Add preceeding nucleotides to FASTA_POS - this is the str index of the middle FASTA nucleotide

                if (
                    (strand == "+")
                    and (range_index == (len(FASTA_range) - 1))
                    and (relative_POS > (len(FASTA) - 4))
                ):  # variant is close to gene UTR start site
                    return [
                        ("Kozak_effect" + ":" + str(relative_POS - len(FASTA))),
                        relative_POS,
                    ]

                else:
                    return [
                        (FASTA[(relative_POS - 2) : (relative_POS + 3)]),
                        relative_POS,
                    ]  # should always equal the REF nucleotide

            else:  # range index = 0
                if (strand == "-") and (relative_POS < 3):
                    return [
                        ("Kozak_effect" + ":" + str(relative_POS - len(FASTA))),
                        relative_POS,
                    ]

                elif (
                    (strand == "+")
                    and (range_index == (len(FASTA_range) - 1))
                    and (relative_POS > (len(FASTA) - 4))
                ):
                    return [
                        ("Kozak_effect" + ":" + str(relative_POS - len(FASTA))),
                        relative_POS,
                    ]

                else:
                    return [
                        (FASTA[(relative_POS - 2) : (relative_POS + 3)]),
                        relative_POS,
                    ]  # should always equal the REF nucleotide


def FASTA_deletion(POS, FASTA_range, FASTA, REF):  # calculates new_FASTA for deletions
    for interval in FASTA_range:
        if POS in interval:
            range_index = FASTA_range.index(interval)  # n of interval

            relative_POS = (POS) - interval.start  # relative POS within interval

            if range_index > 0:
                FASTA_range = FASTA_range[
                    :(range_index)
                ]  # cuts FASTA_range to all ranges before the range_index
                FASTA_POS = [
                    (i.stop - i.start) for i in FASTA_range
                ]  # calculates # of nucleotides in each range
                FASTA_POS = sum(
                    FASTA_POS
                )  # calculates all nucleotides in preceeding ranges

                relative_POS = (
                    relative_POS + FASTA_POS
                )  # Add preceeding nucleotides to FASTA_POS - this is the str index of the middle FASTA nucleotide

                return_str = (
                    (FASTA[: (relative_POS + 1)])
                    + (FASTA[(relative_POS + (len(REF))) :])
                )

                return return_str

            else:
                return_str = (
                    (FASTA[: (relative_POS + 1)])
                    + (FASTA[(relative_POS + (len(REF))) :])
                )

                return return_str


def FASTA_insertion(
    POS, FASTA_range, FASTA, REF, ALT
):  # injects ALT calls into the UTR_sequence
    for interval in FASTA_range:
        if POS in interval:
            range_index = FASTA_range.index(interval)  # n of interval

            relative_POS = (POS) - interval.start  # relative POS within interval

            if range_index > 0:
                FASTA_range = FASTA_range[
                    :(range_index)
                ]  # cuts FASTA_range to all ranges before the range_index
                FASTA_POS = [
                    (i.stop - i.start) for i in FASTA_range
                ]  # calculates # of nucleotides in each range
                FASTA_POS = sum(
                    FASTA_POS
                )  # calculates all nucleotides in preceeding ranges

                relative_POS = (
                    relative_POS + FASTA_POS
                )  # Add preceeding nucleotides to FASTA_POS - this is the str index of the middle FASTA nucleotide

                return_str = (
                    (FASTA[:(relative_POS)]) + ALT + (FASTA[(relative_POS + 1) :])
                )

                return return_str

            else:
                return_str = (
                    (FASTA[:(relative_POS)]) + ALT + (FASTA[(relative_POS + 1) :])
                )

                return return_str


def frame_function(
    START, FASTA_range, FASTA
):  # returns the relative_POS of the uORF START within the UTR
    # START = uORF_start, function is strand-agnostic - must occur after START/STOP sites have been switched back for negative strands

    for interval in FASTA_range:
        if START in interval:
            range_index = FASTA_range.index(interval)  # n of interval

            relative_POS = START - interval.start  # relative POS within interval

            if range_index > 0:
                FASTA_range = FASTA_range[
                    :(range_index)
                ]  # cuts FASTA_range to all ranges before the range_index
                FASTA_POS = [
                    (i.stop - i.start) for i in FASTA_range
                ]  # calculates # of nucleotides in each range
                FASTA_POS = sum(
                    FASTA_POS
                )  # calculates sum of all nucleotides in preceeding ranges

                relative_POS = (
                    relative_POS + FASTA_POS
                )  # Add preceeding nucleotides to FASTA_POS - this is the str index of the middle FASTA nucleotide

                return relative_POS

            else:
                return relative_POS


def splice_function(
    POS, FASTA_range, ALT
):  # identify variants where a deletion may disrupt UTR splicing
    splice_flag = False

    for interval in FASTA_range:
        if POS in interval:
            if (POS + (len(ALT) - 1)) not in interval:
                splice_flag = True

        elif (POS + (len(ALT) - 1)) in interval:
            if POS not in interval:
                splice_flag = True

    return splice_flag


def indel_interpreter(
    start_POS, FASTA
):  # Indels need full ORF context to interpret. This function compresses uORF information to a single codon_change annotation.
    UTR_sequence = FASTA[
        start_POS:
    ]  # truncate the FASTA sequence (old_FASTA or new_FASTA)

    codons = [UTR_sequence[i : i + 3] for i in range(0, len(UTR_sequence), 3)]

    stop_flag = False

    for i in codons:
        if i in ["TAA", "TAG", "TGA"]:
            stop_flag = True

            stop_POS = codons.index(i) + 1

            return codons[0] + "*" + str(int(stop_POS))

            break

    if stop_flag == False:
        try:
            return codons[0] + "*"
        except:
            return "None"  # typically occurs when deletion spans the UTR boundaries


def codon_shift_function(codon_change):
    if "None" not in codon_change:
        left_int = codon_change.split(">")[0].split("*")[1]
        right_int = codon_change.split(">")[1].split("*")[1]

        if (left_int != "") & (right_int != ""):
            return abs((int(left_int) - int(right_int)))

        else:
            return 0

    else:
        return 0


def novel_start_collector(
    old_FASTA, new_FASTA, REF, ALT, strand, relative_POS, UTR_sequence
):
    # this function identifies novel start_codon creation within the UTR
    # potential bug --> a novel start is detected even when an indel does not introduce a new start. IF statement to fiter these out is not very thorough
    # could arguably just be over-interpreting, since my analysis is limited to the lesion site for indels

    new_orf = ""

    if (len(REF) == 1) and (len(ALT) == 1):  # variant is a SNP
        old_codon_list = [old_FASTA[i : i + 3] for i in range(len(old_FASTA) - 2)]
        new_codon_list = [new_FASTA[i : i + 3] for i in range(len(new_FASTA) - 2)]
        adj_POS = relative_POS - 2

    elif len(REF) != 1:  # variant is deletion
        # subset to lesion site

        # adj_relative_POS = (relative_POS - len(REF) + 1)

        if strand == "+":
            ref_lesion_site = UTR_sequence[(relative_POS - 2) : (relative_POS + 3)]
            alt_lesion_site = new_FASTA[(relative_POS - 2) : (relative_POS + 3)]
            adj_POS = relative_POS - 2

        elif strand == "-":
            ref_lesion_site = UTR_sequence[(relative_POS - 2) : (relative_POS + 3)]
            alt_lesion_site = new_FASTA[
                ((relative_POS - (len(REF))) - 1) : ((relative_POS - len(REF)) + 4)
            ]
            adj_POS = (relative_POS - (len(REF))) - 1

        old_codon_list = [
            ref_lesion_site[i : i + 3] for i in range(len(ref_lesion_site) - 2)
        ]
        new_codon_list = [
            alt_lesion_site[i : i + 3] for i in range(len(alt_lesion_site) - 2)
        ]

    elif len(ALT) != 1:  # variant is an insertion
        if strand == "+":
            ref_lesion_site = UTR_sequence[(relative_POS - 2) : (relative_POS + 3)]
            alt_lesion_site = new_FASTA[
                (relative_POS - 2) : (relative_POS + len(ALT) + 2)
            ]
            adj_POS = relative_POS - 2

        elif strand == "-":
            ref_lesion_site = UTR_sequence[(relative_POS - 2) : (relative_POS + 3)]
            alt_lesion_site = new_FASTA[
                (relative_POS - 2) : (relative_POS + (len(ALT)) + 2)
            ]
            adj_POS = relative_POS - 2

        old_codon_list = [
            ref_lesion_site[i : i + 3] for i in range(len(ref_lesion_site) - 2)
        ]
        new_codon_list = [
            alt_lesion_site[i : i + 3] for i in range(len(alt_lesion_site) - 2)
        ]

    # ---------------------------------------------------------------------- #

    for i in range(0, (len(new_codon_list))):
        if i <= (len(old_codon_list) - 1):  # avoids list index error
            if (new_codon_list[i] != old_codon_list[i]) and new_codon_list[
                i
            ] in start_codons:  # if region contains novel start codon
                try:
                    Kozak = (
                        UTR_sequence[adj_POS + i - 3] + UTR_sequence[adj_POS + i + 3]
                    )
                    # grabs -3 and then +4 Kozak sequence, in that order, see wikipedia article for rationale https://en.wikipedia.org/wiki/Kozak_consensus_sequence
                except:
                    Kozak = "error"

                if (len(REF) == 1) and (len(ALT) == 1):  # variant is a SNP
                    UTR = (
                        new_FASTA[i : i + 3] + UTR_sequence[(adj_POS + i + 3) :]
                    )  # pull UTR sequence to search for stop codons

                else:
                    UTR = new_FASTA[(adj_POS + i) :]  # truncate UTR sequence

                codons = [
                    UTR[i : i + 3] for i in range(0, len(UTR), 3)
                ]  # split UTR sequene into in-frame codons

                stop_flag = False

                for x in codons:
                    if x in ["TAA", "TAG", "TGA"]:
                        stop_flag = True

                        stop_POS = codons.index(x) + 1

                        new_orf = (
                            new_orf
                            + ";"
                            + Kozak
                            + "|"
                            + codons[0]
                            + "*"
                            + str(int(stop_POS))
                        )
                        break

                if stop_flag == False:
                    # place in/out annotation here
                    if len(codons[-1]) != 3:
                        frame = "out"

                    else:
                        frame = "in"

                    try:
                        new_orf = (
                            new_orf + ";" + Kozak + "|" + codons[0] + "*" + "/" + frame
                        )
                    except:
                        pass  # typically occurs when deletion spans the UTR boundaries

        elif new_codon_list[i] in start_codons:
            try:
                Kozak = UTR_sequence[adj_POS + i - 3] + UTR_sequence[adj_POS + i + 3]
                # grabs -3 and then +4 Kozak sequence, in that order, see wikipedia article for rationale https://en.wikipedia.org/wiki/Kozak_consensus_sequence
            except:
                Kozak = "error"

            if (len(REF) == 1) and (len(ALT) == 1):  # variant is a SNP
                UTR = (
                    new_FASTA[i : i + 3] + UTR_sequence[(adj_POS + i + 3) :]
                )  # pull UTR sequence to search for stop codons

            else:
                UTR = new_FASTA[(adj_POS + i) :]  # truncate UTR sequence

            codons = [
                UTR[i : i + 3] for i in range(0, len(UTR), 3)
            ]  # split UTR sequene into in-frame codons

            stop_flag = False

            for x in codons:
                if x in ["TAA", "TAG", "TGA"]:
                    stop_flag = True

                    stop_POS = codons.index(x) + 1

                    new_orf = (
                        new_orf
                        + ";"
                        + Kozak
                        + "|"
                        + codons[0]
                        + "*"
                        + str(int(stop_POS))
                    )
                    break

            if stop_flag == False:
                # place in/out annotation here
                if len(codons[-1]) != 3:
                    frame = "out"

                else:
                    frame = "in"

                try:
                    new_orf = (
                        new_orf + ";" + Kozak + "|" + codons[0] + "*" + "/" + frame
                    )

                except:
                    pass  # typically occurs when deletion spans the UTR boundaries

    return new_orf


def complement_function(
    input_FASTA,
):  # This function translates negative strand nucleotides into their complements, but does not reverse the reading frame - must do this manually
    nucleotide_dict = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

    input_FASTA = [nucleotide_dict[k.upper()] for k in input_FASTA]

    new_codon = "".join(input_FASTA)

    return new_codon  # output new codons


def is_whole(n):
    return n % 1 == 0


def site_sensor(
    new_FASTA, site_list
):  # returns True if new_FASTA contains an element in the passed list
    site = False  # must pass new_FASTA as a list structure
    if type(new_FASTA) == list:
        for i in new_FASTA:
            if i in site_list:
                site = True
        return site
    else:
        if new_FASTA in site_list:
            site = True
        return site


# uORF interpreter
# recommended usage: uorfs(vcf, False, False, 0)
def interpret(
    vcf_path,
    Exclude_non_ATG=False,
    Process_compound_variants=False,
    codon_shift_threshold=0,
):
    # define paths for uORF reference data, which uorf_setup.py created in the /build directory
    ORFA_DIR = "/Users/bbowles/Documents/Code/uORF_Code/ORF-Annotation/"
    BUILD_DIR = os.path.join(ORFA_DIR, "build")
    cdspath = os.path.join(BUILD_DIR, "CDS_dictionary.tsv")
    rangepath = os.path.join(BUILD_DIR, "range_df.pkl")

    if type(vcf_path) == pd.core.frame.DataFrame:
        if False in [i in vcf_path.columns for i in ["CHROM", "POS", "REF", "ALT"]]:
            raise pd.errors.ParserError(
                "Input dataframe is missing VCF variant info colums (CHROM, POS, REF, ALT)!"
            )

        VCF_file = vcf_path

    elif type(vcf_path) == str:
        # import vcf - all columns but POS as string
        VCF_file = pd.read_csv(
            vcf_path, comment="#", sep="\t", dtype=str, index_col=False
        )

    if True:  # VCF import, dataframe setup
        # format input dataframe
        VCF_file.loc[:, "POS"] = VCF_file.POS.astype(int)

        # strip whitespace
        VCF_file = VCF_file.apply(
            lambda x: x.astype(str).str.strip() if x.dtype == "object" else x
        )

        # standardize 'None' instances across the dataframe
        VCF_file.replace("None", np.NaN, inplace=True)
        VCF_file.replace(".", np.NaN, inplace=True)

        # rename columns to handle bior quirks (if applicable)
        name_dict = {}
        for i in VCF_file.columns:
            name_dict.update({i: str(i).split(".")[-1]})
        VCF_file.rename(columns=name_dict, inplace=True)

        # Drop uORFs with -inf score
        VCF_file = VCF_file.loc[VCF_file.uORF_score.astype(float) != np.NINF]

        # Define new columns to place annotations into
        VCF_file["effect"] = ""  # effect of new codon, ie "start site loss"
        VCF_file["distance_to_start"] = np.NaN  # distance to canonical gene start site
        VCF_file["start_site"] = np.NaN
        VCF_file["stop_site"] = np.NaN
        VCF_file["frame"] = np.NaN
        VCF_file["compound_error"] = False

        # Convert columns to proper datatypes
        VCF_file.gene_start = VCF_file.gene_start.astype(
            float
        )  # VCF_file gene_start needs to be float to handle datatype issues in UTR iterator
        VCF_file = VCF_file.loc[VCF_file.START.notna()]
        VCF_file.START = VCF_file.START.astype(int)
        VCF_file = VCF_file.loc[VCF_file.STOP.notna()]
        VCF_file.STOP = VCF_file.STOP.astype(int)
        VCF_file = VCF_file.loc[VCF_file.POS.notna()]
        VCF_file.POS = VCF_file.POS.astype(int)

        # Handle multiple calls in a single ALT position
        VCF_file["multi_error"] = (
            False  # This column is a Bool column to track whether the entry was originally multi-allelic, before being exploded
        )
        VCF_file.loc[VCF_file.ALT.str.contains(","), "multi_error"] = True
        VCF_file["ALT"] = VCF_file["ALT"].str.split(",")
        VCF_file = VCF_file.explode(
            "ALT"
        )  # Handles instances where the 'ALT' caller returned multiple results - explodes value and creates additional rows in df

        # To run with BEDtools getFASTA, I needed to switch start and stop codon locations for negative strands (so that START > STOP). This code block reverts that operation.
        pos_df = VCF_file.loc[VCF_file.strand == "+"]  # split dataframe based on strand
        neg_df = VCF_file.loc[VCF_file.strand == "-"]
        neg_df.rename(
            columns={"START": "STOP", "STOP": "START"}, inplace=True
        )  # rename neg_df column headers

        VCF_file = pd.concat([pos_df, neg_df], sort=True)  # re-merge dataframe columns
        VCF_file.sort_index(inplace=True)  # sorts the dataframe on index

        if Exclude_non_ATG == True:
            VCF_file = VCF_file.loc[
                VCF_file.start_codon == "ATG"
            ]  # removes all non-ATG uORFs for hyper-aggressive filtering (can be selected in 'user defined inputs' section)

        # drop incorrectly formatted gene names
        VCF_file = VCF_file.loc[~VCF_file.uORF_ID.str.contains("ENSG")]

        VCF_file.reset_index(
            drop=True, inplace=True
        )  # reset index to avoid issues with the UTR iterator - avoids duplicate index values

    if True:  # load dictionaries
        # load in titer_df which contains 100nt of CDS sequence for the TITER ML tool
        titer_df = pd.read_csv(cdspath, sep="\t", index_col=False)
        titer_df["titer_str"] = titer_df.FASTA.str[0:103]
        titer_df = titer_df[["transcript", "titer_str"]]

        # convert titer_df to dictionary
        titer_dict = {}
        for row in titer_df.itertuples():
            titer_dict.update({row.transcript: row.titer_str})

        # load in CDS dictionary, used to calculate total length of novel uORFs
        CDS_dict = (
            pd.read_csv(cdspath, sep="\t", index_col=False)
            .set_index("transcript")
            .FASTA.to_dict()
        )

        # load in range df, which contains additional, python-object annotations
        range_df = pd.read_pickle(rangepath)

        UTR_range_dict = {}
        UTR_FASTA_dict = {}
        for row in range_df.itertuples():
            UTR_range_dict.update({row.transcript: row.FASTA_range})
            UTR_FASTA_dict.update({row.transcript: row.FASTA})

        # set variables for error tracking
        pos_indel_start_mismatch = 0
        neg_indel_start_mismatch = 0
        pos_indel_ref_mismatch = 0
        neg_indel_ref_mismatch = 0
        pos_SNP_ref_mismatch = 0
        neg_SNP_ref_mismatch = 0

    if True:  # map FASTA sequences
        # Determine UTR FASTA sequence

        # map FASTA annotations
        print("Unpacking UTR FASTA annotation.")
        VCF_file["UTR_sequence"] = (
            VCF_file.uORF_ID.str.split(".").str[0].map(UTR_FASTA_dict)
        )
        VCF_file["UTR_range"] = (
            VCF_file.uORF_ID.str.split(".").str[0].map(UTR_range_dict)
        )

        # handle errors where UTR_sequences or UTR ranges have not successfully mapped - occurs because some Ensembl genes have "ENSTR" prefix instead of "ENST"e
        if (
            not VCF_file.loc[VCF_file.UTR_sequence.isna()].empty
            or not VCF_file.loc[VCF_file.UTR_range.isna()].empty
        ):
            VCF_file.loc[VCF_file.UTR_sequence.isna(), "UTR_sequence"] = (
                "ENST" + VCF_file.uORF_ID.str.split(".").str[0].str.replace("ENSTR", "")
            ).map(UTR_FASTA_dict)
            VCF_file.loc[VCF_file.UTR_range.isna(), "UTR_range"] = (
                "ENST" + VCF_file.uORF_ID.str.split(".").str[0].str.replace("ENSTR", "")
            ).map(UTR_range_dict)

        # Move ALU element variants to a separate dataframe
        ALU_df = VCF_file.loc[
            (VCF_file.REF.str.contains(">")) | (VCF_file.ALT.str.contains(">"))
        ]
        VCF_file = VCF_file.loc[
            ~((VCF_file.REF.str.contains(">")) | (VCF_file.ALT.str.contains(">")))
        ]

        if not ALU_df.empty:
            ALU_df["effect"] = ";transposon movement;"

        # Identify variants that may affect splicing within the UTR
        splice_bool = (
            (VCF_file.REF.str.len() != 1)
            & (VCF_file.ALT.str.len() == 1)
            & (
                VCF_file.apply(
                    lambda row: splice_function(row.START, row.UTR_range, row.ALT),
                    axis=1,
                )
            )
        )
        splice_df = VCF_file.loc[splice_bool.astype(bool)]

        splice_df["effect"] = "UTR_splicing_deletion"
        splice_df["codon_change"] = np.NaN
        VCF_file = VCF_file.loc[np.logical_not(splice_bool)]

        # Drop variants outside of UTR range - currently removes a small subset of indels that occur within -1 of UTR boundaries and could be considered splice site variants
        VCF_file = VCF_file.loc[
            VCF_file.apply(
                lambda row: True
                if True in [row.POS in i for i in row.UTR_range]
                else False,
                axis=1,
            )
        ]
        # VCF_file.loc[np.logical_not(VCF_file.apply(lambda row : True if True in [row.POS in i for i in row.FASTA_range] else False, axis = 1))]
        # ^ debugging only, displays entries cut from df

        # Determine relative_POS of uORF within UTR_sequence - should give location of first nucleotide in START codon

        VCF_file["start_POS"] = np.NaN
        VCF_file.loc[VCF_file.type != "gene_UTR", "start_POS"] = VCF_file.apply(
            lambda row: frame_function(row.START, row.UTR_range, row.UTR_sequence),
            axis=1,
        )

        # drop or separate compound_variants
        if Process_compound_variants == True:
            compound_df = VCF_file.loc[
                (VCF_file.REF.str.len() != 1) & (VCF_file.ALT.str.len() != 1)
            ]  # remove compound variants from VCF file - require separate manipulation
            VCF_file = VCF_file.loc[
                (VCF_file.REF.str.len() == 1) | (VCF_file.ALT.str.len() == 1)
            ]

        else:
            VCF_file = VCF_file.loc[
                (VCF_file.REF.str.len() == 1) | (VCF_file.ALT.str.len() == 1)
            ]

        # calculate new UTR sequence for indels
        indel_df = VCF_file.loc[
            (VCF_file.REF.str.len() != 1) | (VCF_file.ALT.str.len() != 1)
        ]
        VCF_file = VCF_file.loc[
            ((VCF_file.REF.str.len() == 1) & (VCF_file.ALT.str.len() == 1))
        ]

        if not indel_df.empty:  # hot fix for very small VCF files that do not return indels - these are low coverage files that seem unlikely to contain high coverage UTR info
            # determine UTR sequence for deletions - it is matching rel_POS calls to REF calls successfully
            del_df = indel_df.loc[indel_df.REF.str.len() != 1]

            del_df["old_FASTA"] = del_df.UTR_sequence
            del_df["relative_POS"] = (
                (
                    del_df.apply(
                        lambda row: FASTA_function(
                            (row.POS + 1), row.UTR_range, row.UTR_sequence, row.strand
                        ),
                        axis=1,
                    ).str[1]
                )
                - 1
            )

            cistron_df = del_df.loc[del_df.relative_POS.isna()]

            if not (
                (
                    cistron_df.apply(
                        lambda row: True in [row.POS in i for i in row.UTR_range],
                        axis=1,
                    )
                ).empty
            ):
                cistron_df = cistron_df.loc[
                    cistron_df.apply(
                        lambda row: True in [row.POS in i for i in row.UTR_range],
                        axis=1,
                    )
                ]
                cistron_df.loc[:, "effect"] = (
                    ";deletion_near_UTR_bounds;+1 site within intron;"  # deletions identified here are usually at the +1 site within an intron
                )

            del_df = del_df.loc[del_df.relative_POS.notna()]

            del_df["new_FASTA"] = del_df.apply(
                lambda row: FASTA_deletion(
                    row.POS, row.UTR_range, row.UTR_sequence, row.REF
                ),
                axis=1,
            )
            del_df = del_df.loc[del_df.new_FASTA.str.len() >= 5]
            del_df = del_df.loc[~del_df.old_FASTA.str.contains("k")]
            del_df.dropna(subset=["old_FASTA"], inplace=True)
            del_df.dropna(subset=["new_FASTA"], inplace=True)
            # currently drops deletions which had issues with full FASTA retrieval - this occurs when deletion is close to UTR boundaries

            # Determine UTR sequence for insertions - it is matching rel_POS calls to REF calls successfully
            ins_df = indel_df.loc[indel_df.ALT.str.len() != 1]
            ins_df["old_FASTA"] = ins_df.UTR_sequence
            ins_df["relative_POS"] = (
                (
                    ins_df.apply(
                        lambda row: FASTA_function(
                            (row.POS + 1), row.UTR_range, row.UTR_sequence, row.strand
                        ),
                        axis=1,
                    ).str[1]
                )
                - 1
            )

            cistron_df = pd.concat(
                [(ins_df.loc[ins_df.relative_POS.isna()]), cistron_df]
            )

            if not (
                (
                    cistron_df.apply(
                        lambda row: True in [row.POS in i for i in row.UTR_range],
                        axis=1,
                    )
                ).empty
            ):
                cistron_df = cistron_df.loc[
                    cistron_df.apply(
                        lambda row: True in [row.POS in i for i in row.UTR_range],
                        axis=1,
                    )
                ]
                cistron_df.loc[:, "effect"] = (
                    ";deletion_near_UTR_bounds;+1 site within intron;"  # deletions identified here are usually at the +1 site within an intron
                )

            ins_df = ins_df.loc[ins_df.relative_POS.notna()]

            ins_df["new_FASTA"] = ins_df.apply(
                lambda row: FASTA_insertion(
                    row.POS, row.UTR_range, row.UTR_sequence, row.REF, row.ALT
                ),
                axis=1,
            )
            ins_df.dropna(subset=["old_FASTA"], inplace=True)
            ins_df.dropna(subset=["new_FASTA"], inplace=True)

            # concat indel data back together
            indel_df = pd.concat([ins_df, del_df], sort=True)

            # isolate indels that affect Kozak context - uORF_tools.FASTA_function returns a Kozak annotation when nucleotides are close to boundaries
            indel_Kozak = indel_df.loc[indel_df.old_FASTA.str.contains("Kozak")]
            indel_df = indel_df.loc[~indel_df.old_FASTA.str.contains("Kozak")]

        if not VCF_file.empty:
            # calculate FASTA context for SNPs
            VCF_file["FASTA"] = VCF_file.apply(
                lambda row: FASTA_function(
                    row.POS, row.UTR_range, row.UTR_sequence, row.strand
                ),
                axis=1,
            )
            VCF_file["relative_POS"] = VCF_file.FASTA.str[1]
            VCF_file.loc[:, "FASTA"] = VCF_file.FASTA.str[0]

            # collect Kozak context variants from dataframe
            Kozak_df = VCF_file.loc[VCF_file.FASTA.str.len() == 15]
            VCF_file = VCF_file.loc[(VCF_file.FASTA.str.len() == 5)]
            # cuts out: 5'cap variants

            # Determine effect of codon sequence changes --> only considering -1 and -3 positions, literature considers them more important
            Kozak_df["codon_change"] = np.NaN
            Kozak_df["effect"] = Kozak_df.FASTA
            Kozak_df.loc[(Kozak_df.strand == "+"), "codon_change"] = (
                Kozak_df.REF + ">" + Kozak_df.ALT
            )
            Kozak_df.loc[(Kozak_df.strand == "-"), "codon_change"] = (
                Kozak_df.REF.apply(complement_function)
                + ">"
                + Kozak_df.ALT.apply(complement_function)
            )

            # Determine effect for -3 position
            Kozak_df.loc[
                (Kozak_df.FASTA.str[-2:] == "-3")
                & (Kozak_df.codon_change.str[0].isin(["C", "T"]))
                & (Kozak_df.codon_change.str[2].isin(["A", "G"])),
                "effect",
            ] = Kozak_df.effect + ";upregulation"
            Kozak_df.loc[
                (Kozak_df.FASTA.str[-2:] == "-3")
                & (Kozak_df.codon_change.str[0].isin(["A", "G"]))
                & (Kozak_df.codon_change.str[2].isin(["C", "T"])),
                "effect",
            ] = Kozak_df.effect + ";downregulation"

            # Determine effect for -1 position
            Kozak_df.loc[
                (Kozak_df.FASTA.str[-2:] == "-1")
                & (Kozak_df.codon_change.str[0].isin(["A", "T"]))
                & (Kozak_df.codon_change.str[2].isin(["C", "G"])),
                "effect",
            ] = Kozak_df.effect + ";upregulation"
            Kozak_df.loc[
                (Kozak_df.FASTA.str[-2:] == "-1")
                & (Kozak_df.codon_change.str[0].isin(["C", "G"]))
                & (Kozak_df.codon_change.str[2].isin(["A", "T"])),
                "effect",
            ] = Kozak_df.effect + ";downregulation"

            # Keep only variants where Kozak signal is expected to change
            Kozak_df = Kozak_df.loc[Kozak_df.effect.notna()]
            Kozak_df = Kozak_df.loc[Kozak_df.effect.str.contains(";")]

        if not indel_df.empty:
            Kozak_df = pd.concat([Kozak_df, indel_Kozak], sort=True)

        if not VCF_file.empty:
            # Determine old vs new FASTA
            print("Determining old vs new FASTA codons.")

            VCF_file["old_FASTA"] = (
                VCF_file.FASTA
            )  # for many entries, new and old FASTA files will be equal
            VCF_file["new_FASTA"] = np.NaN

            # split dataframe, each chunk needs to be analyzed separately
            POS_df = VCF_file.loc[
                (VCF_file.strand == "+") & (VCF_file.REF.str.len() == 1)
            ]  # positive strand SNPs/insertions
            NEG_df = VCF_file.loc[
                (VCF_file.strand == "-") & (VCF_file.REF.str.len() == 1)
            ]  # negative strand SNPs/insertions

            POS_df.loc[:, "new_FASTA"] = (
                POS_df.FASTA.str[0:2] + POS_df.ALT.astype(str) + POS_df.FASTA.str[3:5]
            )  # positive strand SNPs/insertions
            NEG_df.loc[:, "new_FASTA"] = (
                NEG_df.FASTA.str[0:2] + NEG_df.ALT.astype(str) + NEG_df.FASTA.str[3:5]
            ).apply(lambda x: x[::-1])
            NEG_df.loc[:, "old_FASTA"] = NEG_df.FASTA.apply(lambda x: x[::-1])
            NEG_df.loc[:, "old_FASTA"] = NEG_df.old_FASTA.apply(list).apply(
                complement_function
            )
            NEG_df.loc[:, "new_FASTA"] = NEG_df.new_FASTA.apply(list).apply(
                complement_function
            )  # finished reversing all entries in NEG_df

        if Process_compound_variants == True and not compound_df.empty:
            # isolate any compound variants present - typically result from multi-variant calls at a single site, possibly as a result of poor read mapping
            compound_df.loc[:, ["compound_error"]] = True
            compound_df["old_FASTA"] = compound_df.FASTA
            compound_df["new_FASTA"] = (
                compound_df.FASTA.str[0:2] + compound_df.ALT
            )  # attempt to force variant parsimony
            compound_pos = compound_df.loc[compound_df.strand == "+"]
            compound_neg = compound_df.loc[compound_df.strand == "-"]
            compound_neg.loc[:, "old_FASTA"] = compound_neg.old_FASTA.apply(
                lambda x: x[::-1]
            )
            compound_neg.loc[:, "old_FASTA"] = compound_neg.old_FASTA.apply(list).apply(
                complement_function
            )
            compound_neg.loc[:, "new_FASTA"] = compound_neg.new_FASTA.apply(
                lambda x: x[::-1]
            )
            compound_neg.loc[:, "new_FASTA"] = compound_neg.new_FASTA.apply(list).apply(
                complement_function
            )  # finished reversing all entries in compound_neg

            # drop compound_df entries that could not be corrected
            compound_df = compound_df.loc[
                ~((compound_df.REF.str.len() != 1) & (compound_df.ALT.str.len() != 1))
            ]

        if not indel_df.empty:
            # calculate old and new codons for deletions
            POS_indel = indel_df.loc[indel_df.strand == "+"]  # positive strand indels
            neg_indel = indel_df.loc[indel_df.strand == "-"]  # negative strand indels

            # no manipulation needed for POS strand deletions
            neg_indel.loc[:, "old_FASTA"] = neg_indel.old_FASTA.apply(lambda x: x[::-1])
            neg_indel.loc[:, "new_FASTA"] = neg_indel.new_FASTA.apply(lambda x: x[::-1])
            neg_indel.loc[:, "old_FASTA"] = neg_indel.old_FASTA.apply(list).apply(
                complement_function
            )
            neg_indel.loc[:, "new_FASTA"] = neg_indel.new_FASTA.apply(list).apply(
                complement_function
            )  # finished reversing all entries in indel_df

            # Need to apply reverse complement to UTR sequence - allows the script to determine if UTR mutations cause uORFs to read into downstream CDS, or if they encounter an in-frame stop
            # this adjustment is not currently being performed for SNPs
            neg_indel.loc[:, "UTR_sequence"] = neg_indel.UTR_sequence.apply(list).apply(
                complement_function
            )
            neg_indel.loc[:, "UTR_sequence"] = neg_indel.UTR_sequence.apply(
                lambda x: x[::-1]
            )

        # Need to cat all dataframes back together - NEG_del, POS_del, POS_df, NEG_DF, then reindex

        if Process_compound_variants == True:
            VCF_file = pd.concat(
                [
                    i
                    for i in [
                        neg_indel,
                        POS_indel,
                        POS_df,
                        NEG_df,
                        compound_neg,
                        compound_pos,
                    ]
                    if not i.empty
                ],
                sort=True,
            )

        elif indel_df.empty:
            VCF_file = pd.concat(
                [i for i in [POS_df, NEG_df] if not i.empty], sort=True
            )

        else:
            VCF_file = pd.concat(
                [i for i in [neg_indel, POS_indel, POS_df, NEG_df] if not i.empty],
                sort=True,
            )

        # split dataframe to uORF vs UTR entries
        UTR_file = VCF_file.loc[
            VCF_file.type == "gene_UTR"
        ]  # splits dataframe into 2 files - interpretation of uORF sites depends on frame, while...
        VCF_file = VCF_file.loc[
            VCF_file.type != "gene_UTR"
        ]  # ...interpretation of gene_UTR sites only depends on the surrounding FASTA context

        if VCF_file.start_POS.isna().sum() != 0:
            print(
                "Dropping",
                VCF_file.start_POS.isna().sum(),
                "rows due to missing uORF start codon annotation.",
            )
            print(
                "This is likely caused by different transcript definitions between the .gff file used to build uORF catalog and the McGillivray et al. uORF catalog."
            )
        VCF_file = VCF_file.loc[
            VCF_file.start_POS.notna()
        ]  # drop entries with nan start_POS - may be caused by difference between Refseq and Ensembl UTR definitions
        VCF_file.loc[:, "start_POS"] = VCF_file.start_POS.astype(int)
        VCF_file.loc[:, "relative_POS"] = VCF_file.relative_POS.astype(int)

        print(len(VCF_file.gene.unique()), "unique genes for file.")

    if True:  # begin interpreting intra-uORF variants
        # calculate variant frame

        indel_df = VCF_file.loc[
            (VCF_file.REF.str.len() != 1) | (VCF_file.ALT.str.len() != 1)
        ]
        VCF_file = VCF_file.loc[
            (VCF_file.REF.str.len() == 1) & (VCF_file.ALT.str.len() == 1)
        ]

        # calculate frame for SNPs
        if not VCF_file.empty:
            pos_SNP = VCF_file.loc[VCF_file.strand == "+"]
            neg_SNP = VCF_file.loc[VCF_file.strand == "-"]

            pos_SNP.loc[
                ((pos_SNP.relative_POS - pos_SNP.start_POS) / 3).apply(is_whole),
                "frame",
            ] = 1  # codon_frame = 1, nucleotide at POS is the first one in the codon
            pos_SNP.loc[
                (((pos_SNP.relative_POS - pos_SNP.start_POS) - 1) / 3).apply(is_whole),
                "frame",
            ] = 2  # codon_frame = 1, nucleotide at POS is the first one in the codon
            pos_SNP.loc[pos_SNP.frame.isna(), "frame"] = (
                3  # codon_frame = 3, nucleotide at POS is the third one in the codon
            )

            neg_SNP.loc[
                ((neg_SNP.start_POS - neg_SNP.relative_POS) / 3).apply(is_whole),
                "frame",
            ] = 1  # codon_frame = 1, nucleotide at POS is the first one in the codon
            neg_SNP.loc[
                (((neg_SNP.start_POS - neg_SNP.relative_POS) - 1) / 3).apply(is_whole),
                "frame",
            ] = 2  # codon_frame = 1, nucleotide at POS is the first one in the codon
            neg_SNP.loc[neg_SNP.frame.isna(), "frame"] = (
                3  # codon_frame = 3, nucleotide at POS is the third one in the codon
            )

        # calculate frame for indels (relative_POS + 1 used instead of relative_POS)
        if not indel_df.empty:
            pos_indel = indel_df.loc[indel_df.strand == "+"]
            neg_indel = indel_df.loc[indel_df.strand == "-"]

            # adjust the start position for the pos_indel dataframe to reflect how VCF files handle indel POS definitions
            pos_indel.loc[:, "start_POS"] = pos_indel.start_POS + 1

            # Adjust the neg_indel POS to reflect the reverse compliment
            neg_indel.loc[:, "relative_POS"] = (
                neg_indel.UTR_sequence.str.len() - neg_indel.relative_POS - 1
            )
            neg_indel.loc[:, "start_POS"] = (
                neg_indel.UTR_sequence.str.len() - neg_indel.start_POS - 2
            )  # extra 1-nucleotide adjustment b/c of how VCF files represent indels

            # Confirms that relative_POS and start_POS adjustments return the same nucleotide.
            pos_indel_start_mismatch = pos_indel.apply(
                lambda row: row.UTR_sequence[int(row.start_POS)]
                != (row.start_codon[0]),
                axis=1,
            ).sum()
            if pos_indel_start_mismatch != 0:
                print("New start_POS in neg_indel does not match REF.")
                # raise Processing_Error

            neg_indel_start_mismatch = neg_indel.apply(
                lambda row: row.UTR_sequence[int(row.start_POS)]
                != (row.start_codon[0]),
                axis=1,
            ).sum()
            if neg_indel_start_mismatch != 0:
                print("New start_POS in neg_indel does not match REF.")
                # raise Processing_Error

            # Confirms that VCF REF and nucleotides and uORF catalog nucleotides are the same
            # currently coded hot fix; drop affected regions
            neg_indel_ref_mismatch = neg_indel.apply(
                lambda row: complement_function(row.UTR_sequence[int(row.relative_POS)])
                != (row.REF[0]),
                axis=1,
            ).sum()
            if neg_indel_ref_mismatch != 0:
                print("New relative_POS in neg_indel does not match REF.")
                # for debugging:
                # neg_indel.loc[neg_indel.apply(lambda row : row.UTR_sequence[int(row.start_POS)] != (row.start_codon[0]), axis = 1)]
                neg_indel.loc[
                    neg_indel.apply(
                        lambda row: row.UTR_sequence[int(row.start_POS)]
                        == (row.start_codon[0]),
                        axis=1,
                    )
                ]
            pos_indel_ref_mismatch = pos_indel.apply(
                lambda row: row.UTR_sequence[row.relative_POS] != (row.REF[0]), axis=1
            ).sum()
            if pos_indel_ref_mismatch != 0:
                print("New relative_POS in pos_indel does not match REF.")

            # recat indel dataframes
            indel_df = pd.concat([pos_indel, neg_indel], sort=True)

        # Determine SNP site (START or STOP)
        if not VCF_file.empty:
            # Perform START/STOP indexing for SNP variants
            pos_SNP.loc[
                (
                    (pos_SNP["START"].values <= pos_SNP["POS"].values)
                    & (pos_SNP["POS"].values <= (pos_SNP["START"].values + 2))
                ),
                "start_site",
            ] = True
            pos_SNP.loc[pos_SNP.start_site.isna(), "start_site"] = False
            pos_SNP.loc[
                (
                    ((pos_SNP["STOP"].values - 2) <= pos_SNP["POS"].values)
                    & (pos_SNP["POS"].values <= (pos_SNP["STOP"].values))
                ),
                "stop_site",
            ] = True
            pos_SNP.loc[pos_SNP.stop_site.isna(), "stop_site"] = False
            neg_SNP.loc[
                (
                    ((neg_SNP["START"].values - 2) <= neg_SNP["POS"].values)
                    & (neg_SNP["POS"].values <= (neg_SNP["START"].values))
                ),
                "start_site",
            ] = True
            neg_SNP.loc[neg_SNP.start_site.isna(), "start_site"] = False
            neg_SNP.loc[
                (
                    (neg_SNP["STOP"].values <= neg_SNP["POS"].values)
                    & (neg_SNP["POS"].values <= (neg_SNP["STOP"].values + 2))
                ),
                "stop_site",
            ] = True
            neg_SNP.loc[neg_SNP.stop_site.isna(), "stop_site"] = False

            # concat dataframes back together - downstream analysis is largely agnostic of strand
            SNP_df = pd.concat([pos_SNP, neg_SNP], sort=True)

            # Drop dataframe rows where START and STOP codons overlap - likely represent errors in the original catalog
            SNP_df.drop(
                (
                    SNP_df.loc[
                        (SNP_df.start_site.astype(bool))
                        & (SNP_df.stop_site.astype(bool))
                    ]
                ).index,
                inplace=True,
            )

            # Begin calculating SNP effects

            # subset SNP codons based on uORF reading frame
            SNP_df.loc[SNP_df.frame == 1, "old_FASTA"] = SNP_df.loc[
                SNP_df.frame == 1
            ].old_FASTA.str[2:5]
            SNP_df.loc[SNP_df.frame == 1, "new_FASTA"] = SNP_df.loc[
                SNP_df.frame == 1
            ].new_FASTA.str[2:5]

            SNP_df.loc[SNP_df.frame == 2, "old_FASTA"] = SNP_df.loc[
                SNP_df.frame == 2
            ].old_FASTA.str[1:4]
            SNP_df.loc[SNP_df.frame == 2, "new_FASTA"] = SNP_df.loc[
                SNP_df.frame == 2
            ].new_FASTA.str[1:4]

            SNP_df.loc[SNP_df.frame == 3, "old_FASTA"] = SNP_df.loc[
                SNP_df.frame == 3
            ].old_FASTA.str[0:3]
            SNP_df.loc[SNP_df.frame == 3, "new_FASTA"] = SNP_df.loc[
                SNP_df.frame == 3
            ].new_FASTA.str[0:3]

            # Patch errors in site calling
            SNP_start_mask = (
                (SNP_df.start_site.astype(bool))
                & (np.logical_not(SNP_df.old_FASTA.isin(start_codons)))
            )  # Positions where the reported start site does not match an appropriate start codon
            SNP_df = SNP_df.loc[~SNP_start_mask]  # drops above sites from dataframe

            # correct STOP site entries - some STOP sites were truncated if they overlapped a downstream gene
            SNP_df.loc[
                (SNP_df.stop_site.astype(bool))
                & (np.logical_not(SNP_df.old_FASTA.isin(stop_codons))),
                "stop_site",
            ] = False

            # Interpret variants

            # is variant synonymous?
            SNP_df.loc[(SNP_df.old_FASTA == SNP_df.new_FASTA), "effect"] = (
                SNP_df.effect + ";synonymous"
            )
            SNP_df.loc[
                (
                    (SNP_df.old_FASTA.isin(stop_codons))
                    & (SNP_df.new_FASTA.isin(stop_codons))
                ),
                "effect",
            ] = SNP_df.effect + ";synonymous"

            # Does variant cause loss of start or stop site?
            SNP_df.loc[
                (
                    (SNP_df.stop_site.astype(bool))
                    & (np.logical_not(SNP_df.new_FASTA.isin(stop_codons)))
                ),
                "effect",
            ] = SNP_df.effect + ";stop_loss"

            SNP_df.loc[
                (
                    (SNP_df.start_site.astype(bool))
                    & (np.logical_not(SNP_df.new_FASTA.isin(start_codons)))
                ),
                "effect",
            ] = SNP_df.effect + ";start_loss"

            # Does the variant cause a premature stop?
            SNP_df.loc[
                (
                    (np.logical_not(SNP_df.start_site.astype(bool)))
                    & (np.logical_not(SNP_df.stop_site.astype(bool)))
                    & (SNP_df.new_FASTA.isin(stop_codons))
                ),
                "effect",
            ] = SNP_df.effect + ";premature_stop"

            # Does the variant alter the start context but still preserve similar ribosome binding ability? --> start_missense
            SNP_df.loc[
                (
                    (SNP_df.start_site.astype(bool))
                    & (SNP_df.old_FASTA.isin(low_start_codons))
                    & (SNP_df.new_FASTA.isin(low_start_codons))
                ),
                "effect",
            ] = SNP_df.effect + ";start_missense"
            SNP_df.loc[
                (
                    (SNP_df.start_site.astype(bool))
                    & (SNP_df.old_FASTA.isin(high_start_codons))
                    & (SNP_df.new_FASTA.isin(high_start_codons))
                ),
                "effect",
            ] = SNP_df.effect + ";start_missense"

            # Does the variant upregulate start site ribosome binding ability?
            SNP_df.loc[
                (
                    (SNP_df.start_site.astype(bool))
                    & (SNP_df.old_FASTA.isin(low_start_codons))
                    & (SNP_df.new_FASTA.isin(high_start_codons))
                ),
                "effect",
            ] = SNP_df.effect + ";start_Kozak_upregulation"

            # Does the variant downregulate start site ribosome binding ability?
            SNP_df.loc[
                (
                    (SNP_df.start_site.astype(bool))
                    & (SNP_df.old_FASTA.isin(high_start_codons))
                    & (SNP_df.new_FASTA.isin(low_start_codons))
                ),
                "effect",
            ] = SNP_df.effect + ";start_Kozak_downregulation"

            # remaining variants are predicted low interest - likely missense, although some may alter Kozak sequence or have non-uORF effects
            SNP_df.loc[
                ((SNP_df.effect.isna()) | np.logical_not((SNP_df.effect.astype(bool)))),
                "effect",
            ] = SNP_df.effect + ";missense"

            # Process negative strand SNP_df relative_POS and UTR_sequence: important for later TITER str processing
            # apply reverse complement to SNP_df negative strands, adjust relative_POS
            SNP_df.loc[(SNP_df.strand == "-"), "UTR_sequence"] = (
                SNP_df.UTR_sequence.apply(list).apply(complement_function)
            )
            SNP_df.loc[(SNP_df.strand == "-"), "UTR_sequence"] = (
                SNP_df.UTR_sequence.apply(lambda x: x[::-1])
            )
            SNP_df.loc[(SNP_df.strand == "-"), "relative_POS"] = (
                SNP_df.UTR_sequence.str.len() - SNP_df.relative_POS - 1
            )

            # check VCF reference/catalog reference for positive SNPs
            pos_SNP_ref_mismatch = (
                SNP_df.loc[SNP_df.strand == "+"]
                .apply(
                    lambda row: row.UTR_sequence[row.relative_POS] != row.REF, axis=1
                )
                .sum()
            )
            if pos_SNP_ref_mismatch != 0:
                print(
                    "New relative_POS in negative strand SNP_df sites does not match REF."
                )

            neg_SNP_ref_mismatch = (
                SNP_df.loc[SNP_df.strand == "-"]
                .apply(
                    lambda row: row.UTR_sequence[row.relative_POS]
                    != complement_function(row.REF),
                    axis=1,
                )
                .sum()
            )
            if neg_SNP_ref_mismatch != 0:
                print(
                    "New relative_POS in negative strand SNP_df sites does not match REF."
                )

        # Begin calculating indel effects
        if not indel_df.empty:
            # create a list of new codons introduced by each insertion - used to analyze in-frame insertions
            # Uses POS + 1 string indexing (see below note)
            # relative_POS is the 1st nucleotide in the indel (the unchanged one), but the old_FASTA is now centered on the first nucleotide of the lesion
            # frame is also calculated based on relative_POS + 1 --> these changes are made in an effort to standardize indel reporting across this script

            # treat indels just like SNPs when calling variants
            # trying something new: generate a "start_codon*[distance_to_stop]" annotation for each uORF. An old vs new annotation to compare effect, tack on frame at end
            # can then determine if start codon has changed by grabbing first triplet in sequence
            # just need to make sure I am truly adjusting start_POS and relative_POS for indels

            pos_ins = indel_df.loc[
                (indel_df.strand == "+") & (indel_df.REF.str.len() == 1)
            ]
            neg_ins = indel_df.loc[
                (indel_df.strand == "-") & (indel_df.REF.str.len() == 1)
            ]

            pos_del = indel_df.loc[
                (indel_df.strand == "+") & (indel_df.REF.str.len() != 1)
            ]
            neg_del = indel_df.loc[
                (indel_df.strand == "-") & (indel_df.REF.str.len() != 1)
            ]

            # Perform START/STOP indexing for indel variants (POS is +1 from what would be used in the SNP_df calculation)
            # uses absolute genomic position because surrounding UTR context is not important for this calculation
            # analysis is different for positive and negative strands
            pos_ins.loc[
                (
                    (pos_ins["START"].values <= (pos_ins["POS"].values + 1))
                    & ((pos_ins["POS"].values + 1) <= (pos_ins["START"].values + 2))
                ),
                "start_site",
            ] = True
            pos_ins.loc[pos_ins.start_site.isna(), "start_site"] = False
            pos_ins.loc[
                (
                    ((pos_ins["STOP"].values - 2) <= (pos_ins["POS"].values + 1))
                    & ((pos_ins["POS"].values + 1) <= (pos_ins["STOP"].values))
                ),
                "stop_site",
            ] = True
            pos_ins.loc[pos_ins.stop_site.isna(), "stop_site"] = False
            neg_ins.loc[
                (
                    ((neg_ins["START"].values - 2) <= (neg_ins["POS"].values + 1))
                    & ((neg_ins["POS"].values + 1) <= (neg_ins["START"].values))
                ),
                "start_site",
            ] = True
            neg_ins.loc[neg_ins.start_site.isna(), "start_site"] = False
            neg_ins.loc[
                (
                    (neg_ins["STOP"].values <= (neg_ins["POS"].values + 1))
                    & ((neg_ins["POS"].values + 1) <= (neg_ins["STOP"].values + 2))
                ),
                "stop_site",
            ] = True
            neg_ins.loc[neg_ins.stop_site.isna(), "stop_site"] = False

            # Drop dataframe rows where START and STOP codons overlap - likely represent errors in the original catalog
            pos_ins.drop(
                (
                    pos_ins.loc[
                        (pos_ins.start_site.astype(bool))
                        & (pos_ins.stop_site.astype(bool))
                    ]
                ).index,
                inplace=True,
            )
            neg_ins.drop(
                (
                    neg_ins.loc[
                        (neg_ins.start_site.astype(bool))
                        & (neg_ins.stop_site.astype(bool))
                    ]
                ).index,
                inplace=True,
            )

            # site logic for deletions - determines if there is intersect between deleted nucleotides and START/STOP site for uORF
            # when looking for genomic region intersect, need to expand bounds by +1 on the right side - same issue I've been having with range
            pos_del.loc[:, "start_site"] = pos_del.apply(
                lambda row: bool(
                    (
                        set(range((row.START), ((row.START) + 3))).intersection(
                            set(range((row.POS + 1), (row.POS + len(row.REF))))
                        )
                    )
                ),
                axis=1,
            )
            pos_del.loc[:, "stop_site"] = pos_del.apply(
                lambda row: bool(
                    (
                        set(range((row.STOP - 2), (row.STOP + 1))).intersection(
                            set(range((row.POS + 1), (row.POS + len(row.REF))))
                        )
                    )
                ),
                axis=1,
            )
            neg_del.loc[:, "start_site"] = neg_del.apply(
                lambda row: bool(
                    (
                        set(range((row.START - 2), (row.START + 1))).intersection(
                            set(range((row.POS + 1), (row.POS + len(row.REF))))
                        )
                    )
                ),
                axis=1,
            )
            neg_del.loc[:, "stop_site"] = neg_del.apply(
                lambda row: bool(
                    (set(range(row.STOP, (row.STOP + 3)))).intersection(
                        set(range((row.POS + 1), (row.POS + len(row.REF))))
                    )
                ),
                axis=1,
            )

            # reconcat dataframe
            indel_df = pd.concat([pos_ins, neg_ins, pos_del, neg_del], sort=True)

            # debug check - make sure that interpreter always has CDS uORFs reading past the gene start site
            if (
                indel_df.loc[indel_df.type != "CDSpartial"]
                .apply(
                    lambda row: indel_interpreter(row.start_POS, row.UTR_sequence),
                    axis=1,
                )
                .str[-1]
                == "*"
            ).sum() != 0:
                print("UTRonly type uORFs have a frame error!")
                indel_df = indel_df.loc[
                    (
                        (indel_df.type == "UTRonly")
                        & (
                            indel_df.apply(
                                lambda row: indel_interpreter(
                                    row.start_POS, row.UTR_sequence
                                ),
                                axis=1,
                            ).str[-1]
                            != "*"
                        )
                    )
                    | (
                        (indel_df.type == "CDSpartial")
                        & (
                            indel_df.apply(
                                lambda row: indel_interpreter(
                                    row.start_POS, row.UTR_sequence
                                ),
                                axis=1,
                            ).str[-1]
                            == "*"
                        )
                    )
                ]

            indel_df["codon_change"] = np.NaN
            indel_df.loc[:, "codon_change"] = (
                indel_df.apply(
                    lambda row: indel_interpreter(row.start_POS, row.old_FASTA), axis=1
                )
                + ">"
                + indel_df.apply(
                    lambda row: indel_interpreter(row.start_POS, row.new_FASTA), axis=1
                )
            )

            # Drop synonymous changes from indel_df - variant changes maintained uORF status
            indel_df = indel_df.loc[
                indel_df.codon_change.str.split(">").str[0]
                != indel_df.codon_change.str.split(">").str[1]
            ]

            # interpret start codon changes
            indel_df.loc[:, "effect"] = ""

            # locate any indels that completely remove a start codon
            indel_df.loc[
                (indel_df.start_site)
                & np.logical_not(
                    (
                        indel_df.codon_change.str.split(">")
                        .str[1]
                        .str.split("*")
                        .str[0]
                        .isin(start_codons)
                    )
                ),
                "effect",
            ] = indel_df.effect + ";start_loss"
            # flag these variants if deletions occurs near UTR bounds
            indel_df.loc[
                (indel_df.codon_change.str.contains("None"))
                & (indel_df.start_site)
                & np.logical_not(
                    (
                        indel_df.codon_change.str.split(">")
                        .str[1]
                        .str.split("*")
                        .str[0]
                        .isin(start_codons)
                    )
                ),
                "effect",
            ] = indel_df.effect + ";deletion_near_UTR_bounds"
            indel_df.loc[
                (indel_df.REF.str.len() != 1)
                & (
                    (indel_df.relative_POS) + (indel_df.REF.str.len())
                    >= indel_df.UTR_sequence.str.len()
                )
                & (
                    np.logical_not(
                        indel_df.effect.str.contains("deletion_near_UTR_bounds")
                    )
                ),
                "effect",
            ] = (
                indel_df.loc[
                    (indel_df.REF.str.len() != 1)
                    & (
                        (indel_df.relative_POS) + (indel_df.REF.str.len())
                        >= indel_df.UTR_sequence.str.len()
                    )
                    & (
                        np.logical_not(
                            indel_df.effect.str.contains("deletion_near_UTR_bounds")
                        )
                    )
                ].effect
                + ";deletion_near_UTR_bounds"
            )

            # locate any indels that downregulate a start codon
            indel_df.loc[
                (indel_df.start_site)
                & (
                    indel_df.codon_change.str.split(">")
                    .str[1]
                    .str.split("*")
                    .str[0]
                    .isin(high_start_codons)
                )
                & np.logical_not(
                    (
                        indel_df.codon_change.str.split(">")
                        .str[1]
                        .str.split("*")
                        .str[0]
                        .isin(low_start_codons)
                    )
                ),
                "effect",
            ] = indel_df.effect + ";start_Kozak_downregulation"

            # locate any indels that upregulate a start codon
            indel_df.loc[
                (indel_df.start_site)
                & (
                    indel_df.codon_change.str.split(">")
                    .str[1]
                    .str.split("*")
                    .str[0]
                    .isin(low_start_codons)
                )
                & np.logical_not(
                    (
                        indel_df.codon_change.str.split(">")
                        .str[1]
                        .str.split("*")
                        .str[0]
                        .isin(high_start_codons)
                    )
                ),
                "effect",
            ] = indel_df.effect + ";start_Kozak_upregulation"

            # interpret stop codon changes

            # locate indels that change CDSpartial to UTRonly
            indel_df.loc[
                (
                    (
                        indel_df.codon_change.str.split(">")
                        .str[0]
                        .str.split("*")
                        .str[1]
                        == ""
                    )
                    & (
                        indel_df.codon_change.str.split(">")
                        .str[1]
                        .str.split("*")
                        .str[1]
                        != ""
                    )
                ),
                "effect",
            ] = indel_df.effect + ";stop_introduction;CDSpartial>UTRonly_conversion"

            # locate indels that change UTRonly to CDSpartial
            indel_df.loc[
                (
                    (
                        indel_df.codon_change.str.split(">")
                        .str[0]
                        .str.split("*")
                        .str[1]
                        != ""
                    )
                    & (
                        indel_df.codon_change.str.split(">")
                        .str[1]
                        .str.split("*")
                        .str[1]
                        == ""
                    )
                ),
                "effect",
            ] = indel_df.effect + ";stop_loss;UTRonly>CDSpartial_conversion"

            # locate deletions that change uORF stop codon location (shifts below codon_shift_threshold are discarded)
            indel_df.loc[
                (
                    indel_df.apply(
                        lambda row: codon_shift_function(row.codon_change), axis=1
                    )
                    > codon_shift_threshold
                ),
                "effect",
            ] = indel_df.effect + ";stop_codon_shift"
            # these are difficult to interpret and I will most likely need to use a ML tool such as TITER

            # debug check - locate rows where indel interpreter returns None for an unclear reason. "None" returns that are due to a start site deletion will not throw an exception here.
            if (
                ~indel_df.loc[
                    indel_df.apply(
                        lambda row: indel_interpreter(row.start_POS, row.new_FASTA),
                        axis=1,
                    )
                    == "None"
                ].effect.str.contains("deletion_near_UTR_bounds|start_loss")
            ).sum() != 0:
                print(
                    'indel_interpreter is delivering "None" for rows that do not have start_loss.'
                )
                raise Processing_Error

            # drop entries where no annotation was returned - typically because stop codon was shifted below the codon_shift_threshold
            indel_df = indel_df.loc[indel_df.effect.astype(bool)]

        # completeinterpretation, merge files back together
        if indel_df.empty:
            VCF_file = SNP_df
        elif VCF_file.empty:
            VCF_file = indel_df
        else:
            VCF_file = pd.concat([SNP_df, indel_df], sort=True)

    if True:  # interpret UTR variants
        # Goal: revise uORF approach
        # create a codon_change annotation for each new start codon introduced, feed into a delimited string, then explode into a new dataframe
        # filter new reading frames appropriately - most will need to be fed into TITER or another tool, but I can use select a high interest subset for manual review

        if not UTR_file.empty:
            # handle exceptions caused by empty UTR_file - likely happens due to pre-processing the WES file to isolate exonic regions only

            # drop entries with NaN relative_POS - returned when an indel occurs at the 5' end of the UTR
            UTR_file = UTR_file.loc[UTR_file.relative_POS.notna()]
            UTR_file.loc[:, "relative_POS"] = UTR_file.relative_POS.astype(int)

            UTR_file["codon_change"] = (
                ""  # set codon_change to blank instead of NaN - makes bool filtering easier
            )

            # separate dataframe into SNPs and indels
            SNP_df = UTR_file.loc[
                (UTR_file.REF.str.len() == 1) & (UTR_file.ALT.str.len() == 1)
            ]
            indel_df = UTR_file.loc[
                (UTR_file.REF.str.len() != 1) | (UTR_file.ALT.str.len() != 1)
            ]

            if not SNP_df.empty:
                # apply reverse complement to SNP_df negative strands, adjust relative_POS
                SNP_df.loc[(SNP_df.strand == "-"), "UTR_sequence"] = (
                    SNP_df.UTR_sequence.apply(list).apply(complement_function)
                )
                SNP_df.loc[(SNP_df.strand == "-"), "UTR_sequence"] = (
                    SNP_df.UTR_sequence.apply(lambda x: x[::-1])
                )
                SNP_df.loc[(SNP_df.strand == "-"), "relative_POS"] = (
                    SNP_df.UTR_sequence.str.len() - SNP_df.relative_POS - 1
                )

                # confirm that relative_POS adjustments were successful
                if not SNP_df.loc[SNP_df.strand == "-"].empty:
                    if (
                        SNP_df.loc[SNP_df.strand == "-"].apply(
                            lambda row: complement_function(
                                row.UTR_sequence[row.relative_POS]
                            ),
                            axis=1,
                        )
                        != SNP_df.loc[SNP_df.strand == "-"].REF
                    ).sum() != 0:
                        print("New relative_POS in SNP_df does not match REF.")
                        raise Processing_Error

                # Determine codon_changes - do variants introduce new ORFs within the UTR?
                SNP_df.loc[:, "codon_change"] = SNP_df.apply(
                    lambda row: novel_start_collector(
                        row.old_FASTA,
                        row.new_FASTA,
                        row.REF,
                        row.ALT,
                        row.strand,
                        row.relative_POS,
                        row.UTR_sequence,
                    ),
                    axis=1,
                ).str.strip(";")

                # calculate uORF length for variants that introduce a new start (both UTRonly and CDSpartial regions)
                SNP_df["uORF_length"] = np.NaN
                SNP_df["gene_FASTA"] = np.NaN

                # map gene_FASTA to a unique column
                SNP_df.loc[SNP_df.codon_change.notna(), "gene_FASTA"] = (
                    SNP_df.uORF_ID.map(CDS_dict)
                )

                # apply novel_start_collector (NSC) function
                SNP_df.loc[SNP_df.codon_change.notna(), "uORF_length"] = (
                    SNP_df.loc[SNP_df.codon_change.notna()]
                    .apply(
                        lambda row: novel_start_collector(
                            (row.old_FASTA + row.gene_FASTA),
                            (row.new_FASTA + row.gene_FASTA),
                            row.REF,
                            row.ALT,
                            row.strand,
                            row.relative_POS,
                            (row.UTR_sequence + row.gene_FASTA),
                        ),
                        axis=1,
                    )
                    .str.strip(";")
                )

                # gather codon number from NSC output
                SNP_df.loc[SNP_df.uORF_length.notna(), "uORF_length"] = (
                    SNP_df.loc[SNP_df.uORF_length.notna()]
                    .uORF_length.str.split("*")
                    .str[1]
                )

                # drop gene_FASTA column, no longer needed
                SNP_df.drop(columns=["gene_FASTA"], inplace=True)

                if (
                    SNP_df.loc[SNP_df.codon_change.astype(bool)]
                    .apply(
                        lambda row: row.codon_change.split("|")[1][0:3]
                        not in start_codons,
                        axis=1,
                    )
                    .sum()
                    != 0
                ):
                    print("Error calling UTR_df SNP variants!")
                    raise Processing_Error

                # drop entries where codon_change is empty
                SNP_df = SNP_df.loc[SNP_df.codon_change.astype(bool)]

            if not indel_df.empty:
                # correct relative_POS for negative strand indels
                indel_df.loc[(indel_df.strand == "-"), "relative_POS"] = (
                    indel_df.UTR_sequence.str.len() - indel_df.relative_POS - 1
                )

                # confirm that relative_POS still matches first nucleotide in REF
                if (
                    indel_df.loc[indel_df.strand == "-"]
                    .apply(
                        lambda row: row.UTR_sequence[int(row.relative_POS)]
                        != complement_function(row.REF[0]),
                        axis=1,
                    )
                    .sum()
                ) != 0:
                    print("New relative_POS in neg_indel does not match REF.")
                    raise Processing_Error

                # Determine codon_changes - do variants introduce new ORFs within the UTR?
                indel_df.loc[:, "codon_change"] = indel_df.apply(
                    lambda row: novel_start_collector(
                        row.old_FASTA,
                        row.new_FASTA,
                        row.REF,
                        row.ALT,
                        row.strand,
                        row.relative_POS,
                        row.UTR_sequence,
                    ),
                    axis=1,
                ).str.strip(";")

                # calculate uORF length for variants that introduce a new start (both UTRonly and CDSpartial regions)
                indel_df["uORF_length"] = np.NaN
                indel_df["gene_FASTA"] = np.NaN

                # map gene_FASTA to a unique column
                indel_df.loc[indel_df.codon_change.notna(), "gene_FASTA"] = (
                    indel_df.uORF_ID.map(CDS_dict)
                )

                # apply novel_start_collector (NSC) function
                indel_df.loc[indel_df.codon_change.notna(), "uORF_length"] = (
                    indel_df.loc[indel_df.codon_change.notna()]
                    .apply(
                        lambda row: novel_start_collector(
                            (row.old_FASTA + row.gene_FASTA),
                            (row.new_FASTA + row.gene_FASTA),
                            row.REF,
                            row.ALT,
                            row.strand,
                            row.relative_POS,
                            (row.UTR_sequence + row.gene_FASTA),
                        ),
                        axis=1,
                    )
                    .str.strip(";")
                )

                # gather codon number from NSC output
                indel_df.loc[indel_df.uORF_length.notna(), "uORF_length"] = (
                    indel_df.loc[indel_df.uORF_length.notna()]
                    .uORF_length.str.split("*")
                    .str[1]
                )

                # drop gene_FASTA column, no longer needed
                indel_df.drop(columns=["gene_FASTA"], inplace=True)

                if (
                    indel_df.loc[indel_df.codon_change.astype(bool)]
                    .apply(
                        lambda row: row.codon_change.split("|")[1][0:3]
                        not in start_codons,
                        axis=1,
                    )
                    .sum()
                    != 0
                ):
                    print("Error calling UTR_df indel variants!")
                    raise Processing_Error

                # drop entries where codon_change is empty
                indel_df = indel_df.loc[indel_df.codon_change.astype(bool)]

            # reconcat dataframes
            if SNP_df.empty:
                UTR_file = indel_df
            elif indel_df.empty:
                UTR_file = SNP_df
            else:
                UTR_file = pd.concat([SNP_df, indel_df], sort=True)

            # explode codon_change entries containing multiple calls
            if not UTR_file.empty:
                UTR_file.loc[:, "codon_change"] = UTR_file.codon_change.str.split(
                    ";"
                )  # splits effect entries to a semi-colon delimited list
                UTR_file = UTR_file.explode(
                    "codon_change"
                )  # explodes multiple effect entries to separate columns
                UTR_file = UTR_file.loc[
                    ~UTR_file.effect.isna()
                ]  # drops all UTR_file entries for which there was no predicted effect (missense variants within the 5'UTR)

            # manually classify high impact UTR_effects

            # read Kozak annotation to separate dataframe
            UTR_file["Kozak"] = UTR_file.codon_change.str.split("|").str[0]
            UTR_file.loc[:, "codon_change"] = UTR_file.codon_change.str.split("|").str[
                1
            ]

            # flag CTG and ATG start codons
            UTR_file.loc[
                UTR_file.codon_change.str.split("*").str[0].str.contains("ATG"),
                "effect",
            ] = (
                UTR_file.loc[
                    UTR_file.codon_change.str.split("*").str[0].str.contains("ATG")
                ].effect
                + ";novel_ATG_start"
            )
            UTR_file.loc[
                UTR_file.codon_change.str.split("*").str[0].str.contains("CTG"),
                "effect",
            ] = (
                UTR_file.loc[
                    UTR_file.codon_change.str.split("*").str[0].str.contains("CTG")
                ].effect
                + ";novel_CTG_start"
            )
            UTR_file.loc[
                UTR_file.codon_change.str.split("*").str[0].isin(low_start_codons),
                "effect",
            ] = (
                UTR_file.loc[
                    UTR_file.codon_change.str.split("*").str[0].isin(low_start_codons)
                ].effect
                + ";novel_near_cognate_start"
            )

            # Determine CDSpartial or UTR only
            UTR_file.loc[
                UTR_file.codon_change.str.split("/")
                .str[0]
                .str.split("*")
                .str[1]
                .astype(bool),
                "effect",
            ] = (
                UTR_file.loc[
                    UTR_file.codon_change.str.split("/")
                    .str[0]
                    .str.split("*")
                    .str[1]
                    .astype(bool),
                    "effect",
                ]
                + ";UTRonly"
            )
            UTR_file.loc[
                (
                    np.logical_not(
                        UTR_file.codon_change.str.split("/")
                        .str[0]
                        .str.split("*")
                        .str[1]
                        .astype(bool)
                    )
                )
                & (UTR_file.codon_change.str.contains("/out")),
                "effect",
            ] = (
                UTR_file.loc[
                    (
                        np.logical_not(
                            UTR_file.codon_change.str.split("/")
                            .str[0]
                            .str.split("*")
                            .str[1]
                            .astype(bool)
                        )
                    )
                    & (UTR_file.codon_change.str.contains("/out")),
                    "effect",
                ]
                + ";CDSpartial;out-frame"
            )
            UTR_file.loc[
                (
                    np.logical_not(
                        UTR_file.codon_change.str.split("/")
                        .str[0]
                        .str.split("*")
                        .str[1]
                        .astype(bool)
                    )
                )
                & (UTR_file.codon_change.str.contains("/in")),
                "effect",
            ] = (
                UTR_file.loc[
                    (
                        np.logical_not(
                            UTR_file.codon_change.str.split("/")
                            .str[0]
                            .str.split("*")
                            .str[1]
                            .astype(bool)
                        )
                    )
                    & (UTR_file.codon_change.str.contains("/in")),
                    "effect",
                ]
                + ";CDSpartial;in-frame"
            )

            # Determine Kozak context (1st element in str is -3 position, second element is the +4 position)
            UTR_file.loc[
                (UTR_file.Kozak.str[0].isin(["A", "G"]))
                & (UTR_file.Kozak.str[1] == "G"),
                "Kozak",
            ] = "strong"
            UTR_file.loc[
                (UTR_file.Kozak.str[0].isin(["A", "G"]))
                | (UTR_file.Kozak.str[1] == "G"),
                "Kozak",
            ] = "adequate"
            UTR_file.loc[
                ~UTR_file.Kozak.isin(["strong", "adequate", "error"]), "Kozak"
            ] = "weak"

            UTR_file.loc[(UTR_file.Kozak == "strong"), "effect"] = (
                UTR_file.loc[(UTR_file.Kozak == "strong")].effect + ";strong_Kozak"
            )
            UTR_file.loc[UTR_file.Kozak == "adequate", "effect"] = (
                UTR_file.loc[(UTR_file.Kozak == "adequate")].effect + ";adequate_Kozak"
            )
            UTR_file.loc[UTR_file.Kozak == "weak", "effect"] = (
                UTR_file.loc[(UTR_file.Kozak == "weak")].effect + ";weak_Kozak"
            )

        # merge all dataframes
        concat_list = []

        for name in ["VCF_file", "UTR_file", "splice_df", "Kozak_df"]:
            if name in vars():
                concat_list.append(vars()[name])

        VCF_file = pd.concat(concat_list, sort=True)  # remerge UTR and uORF dataframes
        VCF_file.sort_index(inplace=True)  # sorts the dataframe
        VCF_file.reset_index(drop=True, inplace=True)

    if True:  # post processing
        # TITER processing

        # some elements in titer_str are less than 100 nt - need to investigate further

        # additionally, TITER_str code here does not fill in some genes that have <100nt distance to start :(

        # Updates: neg strand relative_POS and UTR_sequence should be adjusted for all SNPs

        # calculate distance_to_start annotation
        VCF_file.loc[:, "distance_to_start"] = VCF_file.apply(
            lambda row: (len(row.UTR_sequence[int(row.relative_POS) :]) - 1), axis=1
        )

        # map 100 nt of CDS sequence from titer_df
        VCF_file["titer_str"] = np.nan
        VCF_file["titer_ctrl"] = np.nan
        VCF_file.loc[:, "titer_str"] = (
            VCF_file.uORF_ID.str.split(".").str[0].map(titer_dict)
        )
        VCF_file.loc[:, "titer_ctrl"] = (
            VCF_file.uORF_ID.str.split(".").str[0].map(titer_dict)
        )
        # all rows now contain identical 103nt sequence in titer_str and titer_ctrl

        # assemble 203nt TITER_str

        # compile new_UTR sequence for SNP changes within 100nt of canonical gene start
        VCF_file.loc[
            (
                (VCF_file.titer_str.notna())
                & (VCF_file.REF.str.len() == 1)
                & (VCF_file.ALT.str.len() == 1)
                & (VCF_file.strand == "+")
            ),
            "new_FASTA",
        ] = VCF_file.loc[
            (
                (VCF_file.titer_str.notna())
                & (VCF_file.REF.str.len() == 1)
                & (VCF_file.ALT.str.len() == 1)
                & (VCF_file.strand == "+")
            )
        ].apply(
            lambda row: row.UTR_sequence[: int(row.relative_POS)]
            + row.ALT
            + row.UTR_sequence[(int(row.relative_POS) + 1) :],
            axis=1,
        )
        VCF_file.loc[
            (
                (VCF_file.titer_str.notna())
                & (VCF_file.REF.str.len() == 1)
                & (VCF_file.ALT.str.len() == 1)
                & (VCF_file.strand == "-")
            ),
            "new_FASTA",
        ] = VCF_file.loc[
            (
                (VCF_file.titer_str.notna())
                & (VCF_file.REF.str.len() == 1)
                & (VCF_file.ALT.str.len() == 1)
                & (VCF_file.strand == "-")
            )
        ].apply(
            lambda row: row.UTR_sequence[: int(row.relative_POS)]
            + complement_function(row.ALT)
            + row.UTR_sequence[(int(row.relative_POS) + 1) :],
            axis=1,
        )

        # append pre-ATG and post-ATG strings - UTR sequence for negative strands is raw sequence (not reverse complemented)
        VCF_file.loc[(VCF_file.titer_str.notna()), "titer_str"] = (
            VCF_file.loc[(VCF_file.titer_str.notna())].new_FASTA.str[-100:]
            + VCF_file.loc[(VCF_file.titer_str.notna())].titer_str
        )
        VCF_file.loc[(VCF_file.titer_ctrl.notna()), "titer_ctrl"] = (
            VCF_file.loc[(VCF_file.titer_ctrl.notna())].UTR_sequence.str[-100:]
            + VCF_file.loc[(VCF_file.titer_ctrl.notna())].titer_ctrl
        )

        # reset titer_str values of the wrong length
        VCF_file.loc[
            (
                (VCF_file.titer_str.str.len() != 203)
                | (VCF_file.titer_ctrl.str.len() != 203)
            ),
            "titer_str",
        ] = np.nan
        VCF_file.loc[
            (
                (VCF_file.titer_str.str.len() != 203)
                | (VCF_file.titer_ctrl.str.len() != 203)
            ),
            "titer_ctrl",
        ] = np.nan

        # missense variants are not labelled with an explicit effect but may be pertinent to TITER_analysis - label these as "missense"
        # this only occurs for nucleotides within 100nt of the CDS AUG - more distant missense variants remain unlabelled
        VCF_file.loc[
            (
                (VCF_file.titer_str.notna())
                & (VCF_file.titer_str != VCF_file.titer_ctrl)
                & np.logical_not(VCF_file.effect.astype(bool))
            ),
            "effect",
        ] = "missense"

        # Some indels occur in regions with repetitive DNA - label these
        # they could be sequencing artifacts or represent possible repeat expansions?
        VCF_file.loc[
            (
                (VCF_file.titer_str == VCF_file.titer_ctrl)
                & (VCF_file.distance_to_start < 100)
                & (VCF_file.REF.str.len() != 1)
                & (VCF_file.ALT.str.len() != 1)
            ),
            "effect",
        ] = (
            VCF_file.loc[
                (
                    (VCF_file.titer_str == VCF_file.titer_ctrl)
                    & (VCF_file.distance_to_start < 100)
                    & (VCF_file.REF.str.len() != 1)
                    & (VCF_file.ALT.str.len() != 1)
                )
            ].effect
            + ";repetitive_DNA"
        )

        # Assign allele genotype if available

        VCF_file["genotype"] = np.nan

        genotype_flag = False

        # identify columns containing genotype information
        for col in VCF_file.columns.to_list():
            if (VCF_file[col].astype(str).str.contains(r"^1/1:").sum()) > 0:
                genotype_flag = True
                VCF_file["genotype"] = VCF_file[col].str.split(":").str[0]

        # error check for genotype information
        if (
            genotype_flag == True
        ):  # iterator identified a column with genotype information
            if VCF_file.genotype.sum() == 0:
                print("genotype contains nan!")
                raise Processing_Error

        # if no genotype information found, read in a nan genotype column
        else:
            VCF_file["genotype"] = np.nan

        # Add in cistron_df and ALU_df - relative_POS is nan, which causes issues for upstream analysis
        if "cistron_df" in locals():
            VCF_file = pd.concat(
                [VCF_file, cistron_df], sort=True
            )  # remerge UTR and uORF dataframes
        if "ALU_df" in locals():
            VCF_file = pd.concat([VCF_file, ALU_df], sort=True)

        VCF_file.sort_index(inplace=True)  # sorts the dataframe
        VCF_file.reset_index(drop=True, inplace=True)

        # strip excess semicolons
        VCF_file.loc[:, "effect"] = VCF_file.effect.str.strip(";")

        print("\n")
        print("##### ORFA Post-Processing Report #####")
        if (
            sum(
                [
                    pos_indel_start_mismatch,
                    neg_indel_start_mismatch,
                    pos_indel_ref_mismatch,
                    neg_indel_ref_mismatch,
                    pos_SNP_ref_mismatch,
                    neg_SNP_ref_mismatch,
                ]
            )
            != 0
        ):
            print("--------------")
            print("|")
            print(
                "Computed positive indel start sites from uORF catalog that did not match VCF REF nucleotide:",
                pos_indel_start_mismatch,
            )
            print(
                "Computed negative indel start sites from uORF catalog that did not match VCF REF nucleotide:",
                neg_indel_start_mismatch,
            )
            print("|")
            print(
                "Computed positive indel REF sites from uORF catalog that did not match VCF REF nucleotide:",
                pos_indel_ref_mismatch,
            )
            print(
                "Computed negative ndel REF sites from uORF catalog that did not match VCF REF nucleotide:",
                neg_indel_ref_mismatch,
            )
            print("|")
            print(
                "Computed positive SNP REF sites that did not match VCF REF nucleotide:",
                pos_SNP_ref_mismatch,
            )
            print(
                "Computed negative SNP REF sites that did not match VCF REF nucleotide:",
                neg_SNP_ref_mismatch,
            )
            print("|")
            print(
                "##### Mismatch errors are likely a result of different reference genomes used for VCF variant calling and uORF Catalog Construction ######"
            )
            print("|")
            print("--------------")
            print("Mismatch summary:")
            print(
                pd.DataFrame(
                    data={
                        "variable": [
                            "pos_indel_start_mismatch",
                            "neg_indel_start_mismatch",
                            "pos_indel_ref_mismatch",
                            "neg_indel_ref_mismatch",
                            "pos_SNP_ref_mismatch",
                            "neg_SNP_ref_mismatch",
                        ],
                        "value": [
                            pos_indel_start_mismatch,
                            neg_indel_start_mismatch,
                            pos_indel_ref_mismatch,
                            neg_indel_ref_mismatch,
                            pos_SNP_ref_mismatch,
                            neg_SNP_ref_mismatch,
                        ],
                    }
                ).to_string(index=False)
            )
        else:
            print("\n   No catalog build errors detected.\n")
        print("#######################################")

        return VCF_file
