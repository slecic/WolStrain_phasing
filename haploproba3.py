import os
import pandas as pd
import numpy as np
import subprocess
import fileinput
import re

##################################### make counts file from vcf #################################################

def vcf_extract(vcf_file):
    ### extract chrom, pos, ref, alt, ref count and alt count from the vcf file
    cmd=f"bcftools view --max-alleles 3 --types snps {vcf_file} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%RO:%AO]\n'"

    ### write this into a tab delimited file
    with open('file3.txt', 'w') as outfile:
        outfile.write(os.popen(cmd).read()+"\n")

vcf_extract('all.var.wCer1.poly5-0.05.Q20.vcf')
#vcf_extract('var.diplo-0.1.NoWol.snps.q20.noMiss.vcf.recode.vcf')

### names of the first four columns
cols = ['CHROM','POS', 'REF', 'ALT']
### nucleotides
bases= ['A','C','G','T','N']

def vcf_colnames(vcf_file):
    ### extract sample names from the vcf file
    cmd3=f"bcftools query -l {vcf_file}"
    samples = os.popen(cmd3).read().split()
    sampname = cols + samples
    #print(sampname)

    ### add sample names as a header to the file.txt
    with open('file3.txt','r+') as infile:
        content = infile.read()
        infile.seek(0,0)
        infile.write("\t".join(sampname)+"\n"+content)


vcf_colnames('all.var.wCer1.poly5-0.05.Q20.vcf')

### melt the table using custom names for new long columns 'POP' containing population names and 'COUNTS' containing the values  

dt = pd.read_table('file3.txt')
dtm = pd.melt(dt, id_vars=cols, var_name='POP', value_name='COUNTS')
#print(dtm)
dtm[['COUNT1','COUNT2', 'COUNT3']] = dtm['COUNTS'].str.split('[:,]', expand=True)
dtm = dtm.drop(columns = ['COUNTS'])
print(dtm)
dtm.replace('\.+', np.nan, regex=True, inplace=True)
dtm=dtm.dropna(subset=['COUNT1','COUNT2', 'COUNT3'], how='all')
print(dtm)
dtm.fillna(value=0, inplace=True)
print(dtm)

### Add a column for each nucleotide (A, C, G, T) with counts 
def basecountsA(row):

    ### if REF or ALT carry nucleotide A add counts; otherwise add zero
    if row['REF']=='A':
        return (row['COUNT1'])
    elif row['ALT']=='A':
        return (row['COUNT2'])
    elif row['ALT']=='A,C':
        return (row['COUNT2'])
    elif row['ALT']=='A,G':
        return (row['COUNT2'])
    elif row['ALT']=='A,T':
        return (row['COUNT2'])
    elif row['ALT']=='C,A':
        return (row['COUNT3'])
    elif row['ALT']=='G,A':
        return (row['COUNT3'])
    elif row['ALT']=='T,A':
        return (row['COUNT3'])
    else:
        return (0)


def basecountsC(row):

    if row['REF']=='C':
        return (row['COUNT1'])
    elif row['ALT']=='C':
        return (row['COUNT2'])
    elif row['ALT']=='C,A':
        return (row['COUNT2'])
    elif row['ALT']=='C,G':
        return (row['COUNT2'])
    elif row['ALT']=='C,T':
        return (row['COUNT2'])
    elif row['ALT']=='A,C':
        return (row['COUNT3'])
    elif row['ALT']=='G,C':
        return (row['COUNT3'])
    elif row['ALT']=='T,C':
        return (row['COUNT3'])
    else:
        return (0)

def basecountsG(row):
    if row['REF']=='G':
        return (row['COUNT1'])
    elif row['ALT']=='G':
        return (row['COUNT2'])
    elif row['ALT']=='G,A':
        return (row['COUNT2'])
    elif row['ALT']=='G,C':
        return (row['COUNT2'])
    elif row['ALT']=='G,T':
        return (row['COUNT2'])
    elif row['ALT']=='A,G':
        return (row['COUNT3'])
    elif row['ALT']=='C,G':
        return (row['COUNT3'])
    elif row['ALT']=='T,G':
        return (row['COUNT3'])
    else:
        return (0)

def basecountsT(row):
    if row['REF']=='T':
        return (row['COUNT1'])
    elif row['ALT']=='T':
        return (row['COUNT2'])
    elif row['ALT']=='T,A':
        return (row['COUNT2'])
    elif row['ALT']=='T,C':
        return (row['COUNT2'])
    elif row['ALT']=='T,G':
        return (row['COUNT2'])
    elif row['ALT']=='A,T':
        return (row['COUNT3'])
    elif row['ALT']=='C,T':
        return (row['COUNT3'])
    elif row['ALT']=='G,T':
        return (row['COUNT3'])
    else:
        return (0)


def addcounts():

    dtm['countsA'] = dtm.apply (lambda row: basecountsA(row), axis=1)
    dtm['countsC'] = dtm.apply (lambda row: basecountsC(row), axis=1)
    dtm['countsG'] = dtm.apply (lambda row: basecountsG(row), axis=1)
    dtm['countsT'] = dtm.apply (lambda row: basecountsT(row), axis=1)

    return (dtm)

dtm_counts = addcounts()
print(dtm_counts)

### remove MNVs (if called with freebayes; MNPs are passed as type=snp if only one base is changed) where counts are 0 for all 4 nucleotides
indexNames = dtm_counts[(dtm_counts['countsA'] == 0) & (dtm_counts['countsC'] == 0) & (dtm_counts['countsG'] == 0) & (dtm_counts['countsT'] == 0)].index
dtm_counts.drop(indexNames, inplace=True)

### save to a text file
dtm_counts.to_csv(r'counts3.txt', header=None, index=None, sep=' ', mode='w')

### collapse count columns into one colon(:) delimited column
dtm_counts['countstr'] = dtm_counts[dtm_counts.columns[8:]].apply(lambda x: ':'.join(x.dropna().astype(str)),axis=1)
print(dtm_counts)
dtm_counts.drop(['ALT','COUNT1', 'COUNT2', 'COUNT3', 'countsA', 'countsC', 'countsG', 'countsT'], axis = 1, inplace=True)
print(dtm_counts)

### unmelt the data frame
sync = dtm_counts.set_index(['CHROM', 'POS', 'REF', 'POP'])['countstr'].unstack().reset_index()
sync.columns=sync.columns.tolist()
sync.dropna(inplace=True)
sync.to_csv("syncounts3.txt", sep="\t", index=None)
print(sync)

### change column order to match the original data frame
oldcols = list(dt.drop(['ALT'], axis = 1).columns)
print(oldcols)
sync = sync.reindex(columns=oldcols)
sync.to_csv("syncounts13.txt", sep="\t", index=None, header=True)
sync.to_csv("syncounts23.txt", sep="\t", index=None, header=None)
#print(sync)
#sync = dtm_counts.pivot_table(index=['CHROM', 'POS', 'REF'], columns='POP')
#sync.columns = sync.columns.droplevel().rename(None)
#sync.reset_index().fillna("null").to_csv("syncounts.txt", sep="\t", index=None)
#print(sync)
### save to a text file                                                                
#sync.to_csv(r'syncounts.txt', header=None, index=None, sep=' ', mode='w+')
