
import pandas as pd
import numpy as np
import scipy 

# Create one file with all the information we need

info            = 'data/pnas_patient_info.csv'
reads           = 'data/readcounts_96_nodup.tsv'
genes           = 'data/preselectedList'
reads           = 'data/pnas_readcounts_96_nodup.txt'
tpm             = 'data/pnas_tpm_96_nodup.txt'
reads_norm      = 'data/pnas_normal_readcounts.txt'
tpm_norm        = 'data/pnas_normal_tpm.txt'

info            = pd.read_csv(info)
info            = info.reset_index()
info['gene_id'] = info['index']
info            = info.set_index('gene_id')

genes           = pd.read_csv(genes, header=None)
genes           = genes[0].tolist()

# Cancer Patients
reads           = pd.read_csv(reads, sep='\t', header=None)
reads           = reads.T
reads.columns   = reads.iloc[0]
reads           = reads[1:]
cols            = [x for x in reads.columns if 'ENSG' in x]
reads           = reads[cols]

tpm             = pd.read_csv(tpm, sep='\t', header=None)
tpm             = tpm.T
tpm.columns     = tpm.iloc[0]
tpm             = tpm[1:]
cols            = [x for x in tpm.columns if 'ENSG']
tpm             = tpm[cols]

# Normal Patients
reads_norm         = pd.read_csv(reads_norm, sep='\t')
reads_norm         = reads_norm.T
reads_norm.columns = reads_norm.iloc[0]
reads_norm         = reads_norm[1:]
cols               = [x for x in reads_norm.columns if 'ENSG' in x]
reads_norm         = reads_norm[cols]

tpm_norm           = pd.read_csv(tpm_norm, sep='\t')
tpm_norm           = tpm_norm.T
tpm_norm.columns   = tpm_norm.iloc[0]
tpm_norm           = tpm_norm[1:]
cols               = [x for x in tpm_norm.columns if 'ENSG']
tpm_norm           = tpm_norm[cols]


reads['status'] = 1
tpm[  'status'] = 1

reads_norm['status'] = 0
tpm_norm[  'status'] = 0

all_reads = pd.concat([reads, reads_norm])
all_tpm   = pd.concat([tpm,   tpm_norm])

all_reads.to_csv('data/all_reads.csv', index=None)
all_tpm.to_csv(  'data/all_tpm.csv',   index=None)

# For Recurrence
reads['index'] = reads.index-1
reads_info     = pd.merge(info, reads, on='index')
reads_info.to_csv('data/reads_info.csv', index=None)


