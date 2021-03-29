#you may need to install these packages first

library(DESeq2)
library(tximport)
library(tidyverse)
library(vsn)
library(ggplot2)

#quality assessment - assess ALL samples
#define samples
samples <- read_tsv('groups.tsv', col_names = c('sample', 'group'))
#import rsem files
files <- file.path("RSEM_files", paste0(samples$sample, ".rsem.genes.results"))
names(files) <- samples$sample
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

#filter by TPM
gene.list <- rownames(txi.rsem$abundance[rowMaxs(txi.rsem$abundance, value=FALSE)>5,])
txi.rsem$abundance <- txi.rsem$abundance[rowMaxs(txi.rsem$abundance, value=FALSE)>5,]
txi.rsem$counts <- subset(txi.rsem$counts, rownames(txi.rsem$counts) %in% gene.list)
txi.rsem$length <- subset(txi.rsem$length, rownames(txi.rsem$length) %in% gene.list)

#generate sample table
sampleTable <- data.frame(condition = factor(samples$group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

#fix zero lengths
txi.rsem$length[txi.rsem$length == 0] <- 1

#import to deseq2
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

#prefilter (remove very low expressed genes)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#analyze
dds <- DESeq(dds)

## PCA plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition") + theme_bw() 
ggsave("/results/allsamples_PCA.png", width = unit(7, 'in'), height  = unit(8, 'in'))
