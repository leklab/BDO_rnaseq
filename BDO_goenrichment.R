#may need to install these packages

library(DESeq2)
library(tximport)
library(tidyverse)
library(goseq)
library(org.Hs.eg.db)

#format rsem files- remove version numbers, this is done from command line
  #ls -1 *results | while read line
  #do
  #cat $line | sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' > $line.fix
  #done
##

#samples file
samples <- read_tsv('groups.tsv', col_names = c('sample', 'group'))

#sample comparisons
samples_CvsA <- samples[samples$group=="CTL-org" | samples$group=="IL17-org",]

#needed for goseq
#sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' gencode.v26.GRCh38.ERCC.genes.gtf > gencode.v26.GRCh38.ERCC.genes.fix.gtf
txdb <- GenomicFeatures::makeTxDbFromGFF('gencode.v26.GRCh38.ERCC.genes.fix.gtf', format="gtf")
txsByGene=GenomicFeatures::transcriptsBy(txdb,"gene")
lengthData=median(width(txsByGene))

run_goseq <- function(samples,group1,group2,comparison){
  print("import rsem files")
  files <- file.path("RSEM_files", paste0(samples$sample, ".rsem.genes.results.fix")) #notice, now it's the "fix" files, meaning without transcript versions
  names(files) <- samples$sample
  txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  
  print("generate sample table")
  sampleTable <- data.frame(condition = factor(samples$group))
  rownames(sampleTable) <- colnames(txi.rsem$counts)
  print(sampleTable)
  
  print("fix zero lengths")
  txi.rsem$length[txi.rsem$length == 0] <- 1
  
  print("filter by TPM")
  gene.list <- rownames(txi.rsem$abundance[rowMaxs(txi.rsem$abundance, value=FALSE)>5,])
  txi.rsem$abundance <- txi.rsem$abundance[rowMaxs(txi.rsem$abundance, value=FALSE)>5,]
  txi.rsem$counts <- subset(txi.rsem$counts, rownames(txi.rsem$counts) %in% gene.list)
  txi.rsem$length <- subset(txi.rsem$length, rownames(txi.rsem$length) %in% gene.list)
  
  print("import to deseq2")
  dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
  
  print("prefilter (remove very low expressed genes)")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  print("set up comparisons")
  dds$condition <- factor(dds$condition, levels = c(group2,group1))
  
  print("analyze")
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.05, contrast = c("condition",group1,group2)) #alpha is the adjusted p value cutoff, default is 0.1 

  print("format for goseq")
  assayed.genes <- rownames(res)
  fdr.threshold <- 0.05
  sig.genes <- res[!is.na(res$padj) & (res$padj < fdr.threshold),]
  sig.genes.500 <- head(sig.genes[order(sig.genes$padj),],n=500)
  #de.genes <- rownames(sig.genes)[ which(sig.genes$padj < fdr.threshold) ]
  if (nrow(sig.genes) > 500) {
    de.genes <- rownames(sig.genes.500)
    } else { 
      de.genes <- rownames(sig.genes)[ which(sig.genes$padj < fdr.threshold) ]
    }
  print(length(de.genes))
  
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  
  lengthDataDE <- lengthData[intersect(names(lengthData),names(gene.vector))]
  gene.vectorDE <- gene.vector[intersect(names(gene.vector),names(lengthData))]
  
  print("run goseq, quantify the length bias")
  pwf=nullp(gene.vectorDE,"hg38","ensGene", bias.data = lengthDataDE)
  print(head(pwf))
  
  print("run goseq, default method")
  GO.wall=goseq(pwf,"hg38","ensGene")
  print(head(GO.wall, n=20))
  write.table(as.data.frame(GO.wall), file = sprintf("/results/%s/%s_GO.enrichment_results.txt",comparison,comparison), append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

run_goseq(samples_CvsA,"IL17-org","CTL-org","GroupCvsGroupA")


