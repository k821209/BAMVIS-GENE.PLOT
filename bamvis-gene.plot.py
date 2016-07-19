from __future__ import print_function
import numpy as np
import sys
sys.path.append('/ref/analysis/pipelines/')
import kang
import pysam
import re
import pickle
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot') 
import sys
import pandas as pd
import subprocess
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.cm as cm
from tqdm import tqdm 


## data preparation
#file_gff = '/ref/analysis/References/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.fix.intron.gff3'
file_fa  = '/ref/analysis/References/Creinhardtii/Creinhardtii_281_v5.0.fa'
file_bam = '/ref/analysis/Cre/braker/braker.try5_mario/intron3000.merge.sorted.bam'

samfile = pysam.Samfile( file_bam, "rb" )
dicFa   = kang.Fasta2dic(file_fa)

df_gff    = pd.read_pickle('/ref/analysis/pipelines/pandas_df/Creinhardtii_281_v5.5.gene.gff3.pandas.df.pk')
df_gff_ix = df_gff.reset_index().set_index(['genename','longest',2])
## data preparation end


## input data
genename = sys.argv[1] #'Cre05.g238687.v5.5'
chromosome,left,right,strand = list(df_gff_ix.loc[(genename,'1','mRNA')][[0,3,4,6]].values[0])
print(chromosome,left,right,strand)
## input data end


## definitions 
def cigar_parse(cigar):
    match = re.findall(r'(\d+)(\w)', cigar)
    return match
## definitions end 



## gene element array construct
'''
genenamelist = set(df_gff_ix.reset_index()[df_gff_ix.reset_index()[2] == 'mRNA']['genename'])
dic_intron_array = {}
for genename in tqdm(genenamelist):
    df = df_gff_ix.loc[(genename,'1','CDS')].sort_values(by=3).reset_index()
    chromosome  = df[0][0]
    CDSpos      = []
    intron_list = []
    for i in df.index:
        CDSpos.append(df.loc[i][3])
        CDSpos.append(df.loc[i][4])
    for n,p in enumerate(CDSpos):
        if n%2 == 1:
            if n == (len(CDSpos) - 1):
                continue
            intron_list.append([CDSpos[n],CDSpos[n+1]])
    for intLeft, intRight in intron_list:
        try:
            dic_intron_array[chromosome][intLeft-1:intRight] = 1
        except KeyError:
            chrlen     = len(dicFa[chromosome])
            dic_intron_array[chromosome] = np.zeros(chrlen,dtype=np.int)
            dic_intron_array[chromosome][intLeft-1:intRight] = 1
dic_cds_array = {}
for line in tqdm(open(file_gff)):
    if line[0] == '#':
        continue
    cell       = line.strip().split('\t')
    chromosome = cell[0]
    strType    = cell[2]
    intLeft    = int(cell[3])
    intRight   = int(cell[4])
    if strType == 'CDS':
        try:
            dic_cds_array[chromosome][intLeft-1:intRight] = 1
        except KeyError:
            chrlen     = len(dicFa[chromosome])
            dic_cds_array[chromosome] = np.zeros(chrlen,dtype=np.int)
            dic_cds_array[chromosome][intLeft-1:intRight] = 1
'''
## gene element array construct end


## target region RNAseq support array construct

samfile = pysam.Samfile( file_bam, "rb" )
it = samfile.fetch(chromosome,left,right)
dic_tlen_cover_f   = {}
dic_intron_cover_f = {}
dic_match_cover_f  = {}
dic_tlen_cover_r   = {}
dic_intron_cover_r = {}
dic_match_cover_r  = {}
for line in tqdm(it):
    # Check qual
    
    if line.is_proper_pair == False:
        continue
    if line.is_duplicate   == True:
        continue
    if line.is_qcfail      == True:
        continue
    if line.is_secondary   == True:
        continue
    
    # Check qual end
    
    # Check strandness
    ## dUTP strandness : fr_firststrand 
    ## first read  : reverse
    ## second read : forward 
    ## -- > Forward support 
    ## else: Reverse support
    bGOODTOGO = 1
    if line.is_read1 == True and line.is_reverse == True:
        strandness = '+'
    elif line.is_read1 == True and line.is_reverse == False:
        strandness = '-'
    elif line.is_read2 == True and line.is_reverse == True:
        strandness = '-'
    elif line.is_read2 == True and line.is_reverse == False:
        strandness = '+'
    else:
        strandness = '0'
    # Check strandness end      
    
    chromosome   = line.reference_name
    dicFlag      = kang.flagparser(line.flag)
    cigar        = line.cigarstring
    cigarM       = cigar_parse(cigar)
    cigarstrings = [x[1] for x in cigarM]
    cigarvalues  = [x[0] for x in cigarM]
    startP       = line.pos
    fragsize     = line.tlen
    qname        = line.qname
    properpaired = line.is_proper_pair
    length       = startP
    match_list   = [strandness]
    intron_list  = [strandness]
    gap_length   = 0
    for n, cigarstring in enumerate(cigarstrings):
        if  cigarstring == 'M':
            match = [length]
            length += int(cigarvalues[n])
            match.append(length)
            match_list.append(match)
        elif cigarstring == 'I':
            length     += int(cigarvalues[n])
            gap_length += int(cigarvalues[n])
        elif cigarstring == 'N':
            intron = [length]
            if int(cigarvalues[n]) > 5000: # cut-off for maximum intron length
                bGOODTOGO = 0
            length     += int(cigarvalues[n])
            gap_length += int(cigarvalues[n])
            intron.append(length)
            intron_list.append(intron)
    if bGOODTOGO == 0:
        continue
        
    if match_list[0] == '+':
        for l,r in match_list[1:]:
            try:
                dic_match_cover_f[chromosome][l-1:r] += 1
            except KeyError:
                dic_match_cover_f[chromosome] = np.zeros(len(dicFa[chromosome]))
                dic_match_cover_f[chromosome][l-1:r] += 1
    else:
        for l,r in match_list[1:]:
            try:
                dic_match_cover_r[chromosome][l-1:r] += 1
            except KeyError:
                dic_match_cover_r[chromosome] = np.zeros(len(dicFa[chromosome]))
                dic_match_cover_r[chromosome][l-1:r] += 1
    
    if intron_list[0] == '+':
        for l,r in intron_list[1:]:
            try:
                dic_intron_cover_f[chromosome][l-1:r] += 1
            except KeyError:
                dic_intron_cover_f[chromosome] = np.zeros(len(dicFa[chromosome]))
                dic_intron_cover_f[chromosome][l-1:r] += 1
    else:
        for l,r in intron_list[1:]:
            try:
                dic_intron_cover_r[chromosome][l-1:r] += 1
            except KeyError:
                dic_intron_cover_r[chromosome] = np.zeros(len(dicFa[chromosome]))
                dic_intron_cover_r[chromosome][l-1:r] += 1

## target region RNAseq support array construct end


## draw def 

def draw_support(genename,chromosome,left,right,expected_strandness):
    if expected_strandness == '+':
        try:
            sense_intron_support     = dic_intron_cover_f[chromosome][left-1:right]
            antisense_intron_support = dic_intron_cover_r[chromosome][left-1:right]
            sense_cds_support        = dic_match_cover_f[chromosome][left-1:right]
            antisense_cds_support    = dic_match_cover_r[chromosome][left-1:right]
        except KeyError:
            return None
    else:
        try:
            sense_intron_support     = dic_intron_cover_r[chromosome][left-1:right]
            antisense_intron_support = dic_intron_cover_f[chromosome][left-1:right]
            sense_cds_support        = dic_match_cover_r[chromosome][left-1:right]
            antisense_cds_support    = dic_match_cover_f[chromosome][left-1:right]
        except KeyError:
            print('None')
            return None
    fig,ax = plt.subplots(1,figsize=(30,10))
    ax.plot(np.arange(len(antisense_intron_support)),-antisense_intron_support,linewidth=2,label='antisense_intron_support')
    ax.plot(np.arange(len(antisense_cds_support)),-antisense_cds_support,linewidth=2,label='antisense_cds_support')
    ax.plot(np.arange(len(sense_intron_support)),sense_intron_support,linewidth=2,label='sense_intron_support')
    ax.plot(np.arange(len(sense_cds_support)),sense_cds_support,linewidth=2,label='sense_cds_support')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.ylabel('Count per base')
    plt.xlabel('Position in gene')
    plt.tight_layout()
    plt.title("%s, %s:%d-%d (%s)"%(genename,chromosome,left,right,expected_strandness))
    plt.gcf().subplots_adjust(right=0.75)
    print('end')
    plt.savefig("%s.svg"%genename)

draw_support(genename,chromosome,left,right,strand)

