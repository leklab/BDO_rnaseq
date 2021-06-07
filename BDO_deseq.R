### you may need to install these packages first
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(vsn)
library(ggplot2)

#format gtf
gtf <- as_tibble(rtracklayer::import('gencode.v26.GRCh38.ERCC.genes.gtf'))
gtf <- gtf %>% select(gene_id, gene_name) %>% distinct()

#define samples
samples <- read_tsv('groups.tsv', col_names = c('sample', 'group'))

#comparisons
samples_BvsA <- samples[samples$group=="CTL-org" | samples$group=="Org+Th17",]
samples_CvsA <- samples[samples$group=="CTL-org" | samples$group=="IL17-org",]
samples_BvsD <- samples[samples$group=="Org+Th17" | samples$group=="Th17",]

#samples is relevant sample sheet, and comparison is the name of the comparison to be used in file saving
#groups are which 2 groups to be compared, treated group first
run_deseq <- function(samples,group1,group2,comparison){
  dir.create(sprintf("results/%s",comparison)) 
  
  print("import rsem files")
  files <- file.path("RSEM_files", paste0(samples$sample, ".rsem.genes.results"))
  names(files) <- samples$sample
  txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  
  print("filter by TPM")
  gene.list <- rownames(txi.rsem$abundance[rowMaxs(txi.rsem$abundance, value=FALSE)>5,])
  txi.rsem$abundance <- txi.rsem$abundance[rowMaxs(txi.rsem$abundance, value=FALSE)>5,]
  txi.rsem$counts <- subset(txi.rsem$counts, rownames(txi.rsem$counts) %in% gene.list)
  txi.rsem$length <- subset(txi.rsem$length, rownames(txi.rsem$length) %in% gene.list)
  
  print("generate sample table")
  sampleTable <- data.frame(condition = factor(samples$group))
  rownames(sampleTable) <- colnames(txi.rsem$counts)
  print(sampleTable)
  
  print("fix zero lengths")
  txi.rsem$length[txi.rsem$length == 0] <- 1
  
  print("import to deseq2")
  dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
  
  print("prefilter (remove very low expressed genes)")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  print("set up comparisons")
  #need to be levels = c("untreated","treated")
  dds$condition <- factor(dds$condition, levels = c(group2,group1))
  print(dds$condition)
  
  print("analyze")
  dds <- DESeq(dds)
  #for res, should be contrast=c("condition","treated","untreated")
  res <- results(dds, alpha=0.05, contrast = c("condition",group1,group2)) #alpha is the adjusted p value cutoff, default is 0.1 
  print(res)
  print(summary(res))
  sink(file = sprintf("results/%s/%s_results_summary.txt",comparison,comparison))
  summary(res)
  sink(file = NULL)
  
  print("create heatmap")
  rld = rlogTransformation(dds)
  num.genes.plot <- ifelse((sum(res$padj < 0.05, na.rm=TRUE)<50),sum(res$padj < 0.05, na.rm=TRUE),50)
  print(num.genes.plot)  
  mat = assay(rld)[ head(order(res$padj),num.genes.plot), ] # select the top 50 genes with the lowest padj (or less if less significant)
  mat = mat - rowMeans(mat) # Subtract the row means from each value
  # Optional, but to make the plot nicer:
  df = as.data.frame(colData(rld)[,c("condition")]) # Create a dataframe with a column of the conditions
  colnames(df) = "condition" # Rename the column header
  rownames(df) = colnames(mat) # add rownames
  #change the rownames from ENSG to gene names
  rownames_df <- as.data.frame(rownames(mat))
  colnames(rownames_df) <- "gene_id"
  rownames_df <- merge(rownames_df,gtf,by='gene_id', all.x=TRUE)
  rownames_df <- rownames_df[order(rownames_df$gene_id),] #make sure these are both in same order
  mat <- mat[order(row.names(mat)), ] #note this means plot is in order of ENSG number
  rownames(mat) = rownames_df$gene_name #change from ENSG to gene id
  # and plot the actual heatmap
  pheatmap(mat, annotation_col=df, show_rownames = T, legend = F, annotation_names_col = F, filename = sprintf("/results/%s/%s_heatmap.png",comparison,comparison))
    
  print("save")
  res <- as_tibble(res) %>% mutate(gene_id = rownames(res)) %>% left_join(gtf) %>% data.frame() #add gene name
  resOrdered <- as.data.frame(res[order(res$pvalue),])
  write.table(resOrdered, file = sprintf("results/%s/%s_results.txt",comparison,comparison), append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  res.sig <- resOrdered[!is.na(resOrdered$padj) & (resOrdered$padj<0.05),]
  write.table(res.sig, file = sprintf("results/%s/%s_results_sig_only.txt",comparison,comparison), append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  print(nrow(res.sig))
  
  print("Volcano plot")
  rownames(res) <- res$gene_id
  
  EnhancedVolcano(res,
                  lab = res$gene_name,
                  x = 'log2FoldChange',
                  y = 'padj',#other option is padj
                  xlim = c(-10, 10), #edit as required
                  ylim = c(0, 17.5), #edit as required
                  title = sprintf("%s vs %s\n(plotting adjusted p-values)",group1,group2),
                  pCutoff = 0.05, #10e-12 #default is 1e-05
                  #FCcutoff = 2, #default is 1
                  transcriptPointSize = 1.5,
                  transcriptLabSize = 3.0,
                  col=c('black', 'black', 'black', 'red3'),
                  colAlpha = 1, 
                  subtitle = NULL, 
                  legendVisible = F)
  ggsave(sprintf("results/%s/%s_volcano.png",comparison,comparison), width = unit(7, 'in'), height  = unit(8, 'in'))
}

run_deseq(samples_BvsA,"Org+Th17","CTL-org","GroupB_allvsGroupA")
run_deseq(samples_CvsA,"IL17-org","CTL-org","GroupCvsGroupA")
run_deseq(samples_BvsD,"Org+Th17","Th17","GroupBvsGroupD_Tcells")
