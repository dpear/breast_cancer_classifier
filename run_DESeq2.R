install.packages('devtools') #assuming it is not already installed
library(devtools)
install_github('andreacirilloac/updateR')
library(updateR)
updateR(admin_password = 'dsperry')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
 ##### IDK how i installed DESeq2 it was a pain lol


setwd('Google Drive/My Drive/UCSD/q-spring22/Vineet/breast_cancer_classifier/')

counts = 'all_reads.csv'

counts = read.table(counts)

dds = DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy = TRUE)

counts = 'data/pnas_readcounts_96_nodup.txt'
counts = read.csv(counts, sep='\t')

counts_norm = 'data/pnas_normal_readcounts.txt'
counts_norm = read.csv(counts_norm, sep='\t')

all_counts = merge(counts, counts_norm)

samples = colnames(all_counts)
samples = samples[2:length(samples)]
samples
dex     = c(rep(1,96), rep(0,32))
meta    = data.frame(samples=samples, dex=dex)


