# BDO_rnaseq

RNA sequencing data processing and analysis of bile-derived organoids

#Processing on Yale HPC cluster

The FASTQ files were aligned to hg38 human reference genome with GENCODE v26 annotations by STAR (using the pipeline described at https://github.com/leklab/RNAseq), and RSEM was used to quantify gene expression levels from the STAR-aligned bam files.

#Downstream analysis in R

DESeq2 was used to identify differentially expressed genes, defined as those with an adjusted p-value of <0.05. Lowly expressed genes with Transcripts Per Kilobase Million (TPM) <5 in all samples, as well as those with a total count <10 across all samples, were filtered prior to differential gene expression analysis. 

Goseq was used for gene ontology analysis to identify pathways enriched in, or depleted of, significant changes in gene expression. 


