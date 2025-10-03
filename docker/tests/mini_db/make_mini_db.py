#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 22:20:50 2025

Goal: make mini input files for testing database builds.

NRG1 = chr8
MEF2C = chr5

For the purposes of this script, exon = 5'UTR portion of the exon, and CDS = 
    coding portion of the exon. Together, "exon" + "CDS" forms what is the traditional
    exon definition.

@author: bbowles
"""

import pandas as pd
import math
import subprocess

def return_position(start, fasta):
    return start + len(fasta)

def import_gtf(gtf_path):
    
    # import gtf
    ensg_df = pd.read_csv(gtf_path, 
                        header=None, 
                        sep='\t', comment='#', names=['seqname', 'source', 'feature', 
                                                        'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    return ensg_df
    


def chunk_str(fasta_str, size):
    # used to write FASTA lines of set len
    if len(fasta_str)>size:
        ceiling = math.ceil(len(fasta_str)/size)
        chunk_str = [fasta_str[0+(i*size):(i+1)*size]+"\n" if i!= ceiling else fasta_str[0+(i*size):] for i in range(ceiling) ]
        return "".join(chunk_str)
    else:
        return fasta_str

    

def complement_function(input_FASTA):  # This function translates negative strand nucleotides into their complements, but does not reverse the reading frame - must do this manually
    """FASTA string

    :param input_FASTA: FASTA nucleotide sequence
    :type input_FASTA: str
    :return: Output nucleotide string converted to reverse compliment seq
    :rtype: str
    """
    nucleotide_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
        
    input_FASTA = [nucleotide_dict[k.upper()] for k in input_FASTA]
       
    new_codon = ''.join(input_FASTA)
    
    return new_codon        # output new codons


####################
# MAKE MEF2C FASTA #
####################

# pad sequence by 10bp on either side
padding = "N"*10

# MEF2C - ENST00000504921.7
mef2c_exon_1 = """GCAGTCACAGACACTTGAGCACACGCGTACACCCAGACATCTTCGGGCTGCTATTGGATT
GACTTTGAAGGTTCTGTGTGGGTCGCCGTGGCTGCATGTTTGAATCAGGTGGAGAAGCAC
TTCAACGCTGGACGAAGTAAAGATTATTGTTGTTATTTTTTTTTTCTCTCTCTCTCTCTC
TTAAGAAAGGAAAATATCCCAAGGACTAATCTGATCGGGTCTTCCTTCAT""".replace('\n','').strip(' ').upper()

mef2c_exon_2 = """CAGGAACGAATGCAGGAATTTGGGAACTGAGCTGTGCAAGTGCTGAAGAAGGAGATTTGT
TTGGAGGAAACAGGAAAGAGAAAGAAAAGGAAGGAAAAAATACATAATTTCAGGGACGAG
AGAGAGAAGAAAAACGGGGACT""".replace('\n','').strip(' ').upper()

# NOT THE ACTUAL SEQ, NEED TO UPDATE
mef2c_cds = """ATGGGGAGAAAAAAGATTCAGATTACGAGGATTATGGA
TGAACGTAACAGACAG""".replace('\n','').strip(' ').upper()

# reverse sequences (convert mrna to dna)
mef2c_fasta = complement_function(
    mef2c_cds)[::-1] + complement_function(mef2c_exon_2)[::-1] + complement_function(mef2c_exon_1)[::-1]

# append padding
mef2c_fasta = padding + mef2c_fasta + padding

# get header
"""
>1 dna_sm:chromosome chromosome:GRCh37:1:1:249250621:1 REF
coord_system:version:name:start:end:strand
"""
mef2c_chrom = "5"
mef2c_header = f">{mef2c_chrom} dna_sm:chromosome chromosome:GRCh37:{mef2c_chrom}:1:{len(mef2c_fasta)+1}:1 REF"
mef2c_header = f">5 dna:chromosome chromosome:GRCh38:5:1:{len(mef2c_fasta)+1}:1 REF"
mef2c_header = ">5 dna:chromosome chromosome:GRCh38:5:1:181538259:1 REF"
nrg1_header = ">5 dna:chromosome chromosome:GRCh38:5:1:181538259:1 REF"


mef2c_df = pd.DataFrame([
    ['mef2c','ENST00000504921.7',1,'exon', mef2c_exon_1],
    ['mef2c','ENST00000504921.7',2,'exon', mef2c_exon_2+mef2c_cds],
    ['mef2c','ENST00000504921.7',1,'utr', mef2c_exon_1],
    ['mef2c','ENST00000504921.7',2,'utr', mef2c_exon_2],
     ['mef2c','ENST00000504921.7',2,'cds', mef2c_cds]
    ], columns=['gene','transcript','exon', 'type', 'FASTA'])
mef2c_df["length"] = mef2c_df.FASTA.str.len()


# convert sequence coordinates to test coords
fasta_len = len(mef2c_fasta)
target_exon = 1
mef2c_df.loc[(mef2c_df.type=="exon") & (mef2c_df.exon==target_exon), 'start'] = fasta_len - 9 - len(mef2c_exon_1)
mef2c_df.loc[(mef2c_df.type=="exon") & (mef2c_df.exon==target_exon), 'end'] = fasta_len - 10

target_exon = 2
mef2c_df.loc[(mef2c_df.type=="exon") & (mef2c_df.exon==target_exon), 'start'] = fasta_len - 9 - len(mef2c_exon_1) - len(mef2c_exon_2) - len(mef2c_cds)
mef2c_df.loc[(mef2c_df.type=="exon") & (mef2c_df.exon==target_exon), 'end'] = fasta_len - 10 - len(mef2c_exon_1)

target_exon = 1
mef2c_df.loc[(mef2c_df.type=="utr") & (mef2c_df.exon==target_exon), 'start'] = fasta_len - 9 - len(mef2c_exon_1)
mef2c_df.loc[(mef2c_df.type=="utr") & (mef2c_df.exon==target_exon), 'end'] = fasta_len - 10

target_exon = 2
mef2c_df.loc[(mef2c_df.type=="utr") & (mef2c_df.exon==target_exon), 'start'] = fasta_len - 9 - len(mef2c_exon_1) - len(mef2c_exon_2)
mef2c_df.loc[(mef2c_df.type=="utr") & (mef2c_df.exon==target_exon), 'end'] = fasta_len - 10 - len(mef2c_exon_1)


mef2c_df.loc[(mef2c_df.type=="cds") , 'start'] = 11
mef2c_df.loc[(mef2c_df.type=="cds") , 'end'] = fasta_len - 10 - len(mef2c_exon_1) - len(mef2c_exon_2)


print("MEF2C DATA:\n")
print(mef2c_df.to_string(index=False))

mef2c_df[['start','end']]




###################
# MAKE NRG1 FASTA #
###################

# NRG1 - ENST00000405005.8

nrg1_exon_1 = """AAACTTGTTGGAACTCCGGGCTCGCGCGGAGGCCAGGAGCTGAGCGGCGGCGGCTGCCGG
ACGATGGGAGCGTGAGCAGGACGGTGATAACCTCTCCCCGATCGGGTTGCGAGGGCGCCG
GGCAGAGGCCAGGACGCGAGCCGCCAGCGGTGGGACCCATCGACGACTTCCCGGGGCGAC
AGGAGCAGCCCCGAGAGCCAGGGCGAGCGCCCGTTCCAGGTGGCCGGACCGCCCGCCGCG
TCCGCGCCGCGCTCCCTGCAGGCAACGGGAGACGCCCCCGCGCAGCGCGAGCGCCTCAGC
GCGGCCGCTCGCTCTCCCCCTCGAGGGACAAACTTTTCCCAAACCCGATCCGAGCCCTTG
GACCAAACTCGCCTGCGCCGAGAGCCGTCCGCGTAGAGCGCTCCGTCTCCGGCGAG""".replace('\n','').strip(' ').upper()

nrg1_cds = """ATGT
CCGAGCGCAAAGAAGGCAGAGGCAAAGGGAAGGGCAAGAAGAAGGAGCGAGGCTCCGGCA
AGAAGCCGGAGTCCGCGGCGGGCAGCCAGAGCCCAG""".replace('\n','').strip(' ').upper()

nrg1_fasta = padding + nrg1_exon_1 + nrg1_cds + padding


# get header
nrg1_chrom = "8"
nrg1_header = f">{nrg1_chrom} dna_sm:chromosome chromosome:GRCh38:{nrg1_chrom}:1:{len(nrg1_fasta)+1}:1 REF"
nrg1_header = f">8 dna:chromosome chromosome:GRCh38:8:1:{len(nrg1_fasta)+1}:1 REF"
nrg1_header = ">8 dna:chromosome chromosome:GRCh38:8:1:145138636:1 REF"

nrg1_df = pd.DataFrame([
    ['nrg1','ENST00000405005.8',1,'exon',nrg1_exon_1+nrg1_cds],
    ['nrg1','ENST00000405005.8',1,'utr',nrg1_exon_1],
    ['nrg1','ENST00000405005.8',1,'cds',nrg1_cds]
    ], columns=['gene','transcript','exon', 'type', 'FASTA'])

# adjust coordinates
fasta_len = len(nrg1_fasta)
nrg1_df.loc[nrg1_df.type=="exon", 'start'] = 11
nrg1_df.loc[nrg1_df.type=="exon", 'end'] = fasta_len - 10

nrg1_df.loc[nrg1_df.type=="utr", 'start'] = 11
nrg1_df.loc[nrg1_df.type=="utr", 'end'] = fasta_len - len(nrg1_cds) - 10

nrg1_df.loc[nrg1_df.type=="cds", 'start'] = 11 + len(nrg1_exon_1)
nrg1_df.loc[nrg1_df.type=="cds", 'end'] = fasta_len - 10


print("NRG1 DATA:\n")
print(nrg1_df.to_string(index=False))

 
# write FASTA output
outfasta = "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/mini_db/minifasta.fa"
with open(outfasta, "w") as f:
    lines = [chunk_str(i,60) for i in [
        mef2c_header, '\n', mef2c_fasta, '\n', nrg1_header, '\n', nrg1_fasta]]
    f.writelines(lines)
        
# Index the FASTA file
subprocess.run([
    'samtools',
    'faidx',
    outfasta
])


############
# MAKE GTF #
############
target_transcripts = {"mef2c":"ENST00000504921.7",
                      "nrg1":"ENST00000405005.8"}

gtf_path = "/Users/bbowles/Documents/Code/refdata/MANE/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"
ensg_df = import_gtf(gtf_path)
ensg_cols = ensg_df.columns

# subset to havana build
ensg_df = ensg_df.loc[ensg_df.source=="ensembl_havana"]

# create subset dataframes
cds = ensg_df.loc[ensg_df.feature=="CDS"].copy()
exons = ensg_df.loc[ensg_df.feature=='exon'].copy()
utr_df = ensg_df.loc[ensg_df.feature=='UTR'].copy()

# annotate transcripts
cds["transcript"] = cds.attribute.str.split(";").str[1].str.split(' ').str[2].str.strip('"')
exons['transcript'] = exons.attribute.str.split(';').str[1].str.split(' ').str[2].str.strip('"')
utr_df['transcript'] = utr_df.attribute.str.split(';').str[1].str.split(' ').str[2].str.strip('"')
cds["type"] = "cds"
utr_df["type"] = "utr"
exons["type"] = "exon"


# select required CDS entries
cds = cds.loc[cds.transcript.isin(target_transcripts.values())]
cds.sort_values(by=["seqname","start","end"], ascending=True, inplace=True)

# take first cds for NRG1, last for MEF2C
cds = pd.DataFrame([cds.loc[cds.transcript == target_transcripts['mef2c']].iloc[-1,:],
          cds.loc[cds.transcript == target_transcripts['nrg1']].iloc[0,:]])
cds["exon"] = cds.attribute.str.split(';').str[6].str.split(' ').str[2].astype(int)

# select required UTR entries
utr_df = utr_df.loc[utr_df.transcript.isin(target_transcripts.values())]
utr_df.sort_values(by=["seqname","start","end"], ascending=True, inplace=True)

# get UTR exon count
utr_df["exon"] = utr_df.attribute.str.split(';').str[6].str.split(' ').str[2].astype(int)

# take first UTR entry for NRG1, first two for MEF2C
utr_df.sort_values(by=["seqname","transcript","exon"], ascending=True, inplace=True)

mef2c_utrs = utr_df.loc[(utr_df.transcript == target_transcripts['mef2c']) & 
                        (utr_df.exon.isin([1,2]))]
nrg1_utr = utr_df.loc[utr_df.transcript == target_transcripts['nrg1']].iloc[0,:]
utr_df = pd.concat([mef2c_utrs,
          pd.DataFrame(nrg1_utr).transpose()],
                  axis=0)

# select required exon positions
exons = exons.loc[exons.transcript.isin(target_transcripts.values())]
exons.sort_values(by=["seqname","start","end"], ascending=True, inplace=True)

# get exon count
exons["exon"] = exons.attribute.str.split(';').str[6].str.split(' ').str[2].astype(int)

# take first UTR entry for NRG1, first two for MEF2C
mef2c_exons = exons.loc[(exons.transcript == target_transcripts['mef2c']) & 
                        (exons.exon.isin([1,2]))]
nrg1_exon = exons.loc[exons.transcript == target_transcripts['nrg1']].iloc[0,:]
exons = pd.concat([mef2c_exons,
          pd.DataFrame(nrg1_exon).transpose()],
                  axis=0)

# concat entries
output_gtf = pd.concat([utr_df, exons, cds]).sort_values(by=["seqname","start","end"], ascending=True)

# adjust coordinates
adjust_df = pd.concat([nrg1_df, mef2c_df])[["gene","transcript","type","exon","start","end"]]
adjust_df.rename(columns={"start":"new_start","end":"new_end"}, inplace=True)
adjust_df.loc[:,'exon'] = adjust_df.exon.astype(int)
output_gtf.loc[:,'exon'] = output_gtf.exon.astype(int)



output_gtf = pd.merge(output_gtf, adjust_df, how="outer")
output_gtf = output_gtf.drop(columns=["start","end"]).rename(columns={"new_start":"start","new_end":"end"})

# drop empty rows
output_gtf.dropna(inplace=True)

# collect crosstabs as a check
crosstab = pd.crosstab(output_gtf.feature, output_gtf.transcript)

# format output gtf
output_gtf = output_gtf[ensg_cols]

output_gtf.sort_values(by=['seqname','start','end'], ascending = True)

print(output_gtf.shape)

outpath = "/Users/bbowles/Documents/Code/GitHub/d3-uORF-Viewer/tests/mini_db/mini.gtf.gz"
output_gtf.to_csv(outpath, 
                  index=False, sep='\t', compression="gzip", header=None)

print(import_gtf(outpath).shape)


