# BDO_rnaseq

RNA sequencing data processing and analysis of bile-derived organoids

**Processing on Yale HPC cluster**

The FASTQ files were aligned to hg38 human reference genome with GENCODE v26 annotations by STAR using the pipeline described at https://github.com/leklab/RNAseq - `gencode.v26.GRCh38.annotation.gtf` is the GTF from GENCODE downloaded from https://www.gencodegenes.org/human/release_26.html. The following additional processing step is required for Goseq:

```sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' gencode.v26.GRCh38.ERCC.genes.gtf > gencode.v26.GRCh38.ERCC.genes.fix.gtf
```

RSEM was used to quantify gene expression levels from the STAR-aligned bam files (`rsem.sh`). The ```rsem.genes.results``` output files were used for downstream analyses. The following additional processing step is required for Goseq:

```ls -1 *results | while read line
  do
  	cat $line | sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' > $line.fix
  done
```

**Downstream analysis in R (using v3.6.1)**

[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) was used to identify differentially expressed genes, defined as those with an adjusted p-value of <0.05. Lowly expressed genes with Transcripts Per Kilobase Million (TPM) <5 in all samples, as well as those with a total count <10 across all samples, were filtered prior to the differential gene expression analysis (`BDO_deseq.R`).

[Goseq](https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf) was used for gene ontology analysis to identify pathways enriched in, or depleted of, significant changes in gene expression (`BDO_goenrichment.R`). 

**Requirements**

Using HPC:
`STAR`
`Picard`
`SAMtools`
`Python`
`RSEM`

In R:
`DESeq2`
`GoSeq`
`org.Hs.eg.db`
`tximport`
`tidyverse`
`ggplot2`
`pheatmap`
`EnhancedVolcano`
`rtracklayer`

DESeq2, tximport, EnhancedVolcano, and rtracklayer were installed through the bioconductor package manager.
