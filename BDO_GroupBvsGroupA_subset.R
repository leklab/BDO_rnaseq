### you may need to install these packages first
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(vsn)
library(ggplot2)
library(goseq)
library(org.Hs.eg.db)
library(dplyr)

#group D - Th17 alone, group A - control organoids alone, group B - organoids + Th17 cells
#this is a specialized analysis to identify genes upregulated in group B vs group A, and in group B vs group D
GrpBvsGrpD <- read.table("/results/GroupBvsGroupD_Tcells/GroupBvsGroupD_Tcells_results.txt", header=TRUE)
GrpBvsGrpA <- read.table("/results/GroupB_allvsGroupA/GroupB_allvsGroupA_results.txt", header=TRUE)
#subset
GrpBvsGrpD_subset <- GrpBvsGrpD[GrpBvsGrpD$log2FoldChange>1,]
#subset groupB vs groupA by genes that are upregulated in B vs D, then subset to genes increased in groupB vs groupA
GrpBvsGrpA_subset <- GrpBvsGrpA[which(GrpBvsGrpA$gene_id %in% GrpBvsGrpD_subset$gene_id),]
GrpBvsGrpA_subset <- GrpBvsGrpA_subset[GrpBvsGrpA_subset$log2FoldChange>1 & order(GrpBvsGrpA_subset$padj),]
#write to file
write.table(GrpBvsGrpA_subset, file = "/results/GroupBvsGroupA_subset/GroupBvsGroupA_subset_results.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#volcano plot
EnhancedVolcano(GrpBvsGrpA_subset,
                lab = GrpBvsGrpA_subset$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(0, 15),
                ylim = c(0, 20),
                title = "GroupB vs GroupA upregulated genes\nthat are also upregulated in GroupB vs GroupD (plotting adjusted p-values)",
                pCutoff = 0.05, #10e-12 #default is 1e-05
                #FCcutoff = 2, #default is 1
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1, 
                subtitle = NULL, 
                legendVisible = F,
                drawConnectors = TRUE,
                pLabellingCutoff = 1e-7) #this is equal to -log10P of 7
ggsave("/results/GroupBvsGroupA_subset/GroupBvsGroupA_subset_volcano.png", width = unit(10, 'in'), height  = unit(8, 'in'))

#to get the heatmap, need the dds object
gtf <- as_tibble(rtracklayer::import('gencode.v26.GRCh38.ERCC.genes.gtf'))
gtf <- gtf %>% dplyr::select(gene_id, gene_name) %>% distinct()
samples <- read_tsv('groups.tsv', col_names = c('sample', 'group'))
samples_BvsA <- samples[samples$group=="CTL-org" | samples$group=="Org+Th17",]

run_deseq_modified <- function(samples,group1,group2,comparison){
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
  return(dds)
} 

dds <- run_deseq_modified(samples_BvsA,"Org+Th17","CTL-org","GroupB_allvsGroupA")
res <- results(dds, alpha=0.05, contrast = c("condition","Org+Th17","CTL-org")) #alpha is the adjusted p value cutoff, default is 0.1 
#heatmap
rld = rlogTransformation(dds)
mat = assay(rld)[order(res$padj),] 
mat = mat - rowMeans(mat) # Subtract the row means from each value
##subset mat to only include genes that are significantly upregulated in GrpBvsGrpA and upregulated in GroupBvsGroupD
mat = subset(mat, rownames(mat) %in% GrpBvsGrpA_subset[GrpBvsGrpA_subset$padj<0.05 & !is.na(GrpBvsGrpA_subset$padj),c('gene_id')])
# Optional, but to make the plot nicer:
df = as.data.frame(colData(rld)[,c("condition")]) # Create a dataframe with a column of the conditions
colnames(df) = "condition" # Rename the column header
df$condition <- paste(df$condition,"  ",sep="") #trying to avoid cropping of legend
rownames(df) = colnames(mat) # add rownames
#change the rownames from ENSG to gene names
rownames_df <- as.data.frame(rownames(mat))
colnames(rownames_df) <- "gene_id"
rownames_df <- merge(rownames_df,gtf,by='gene_id', all.x=TRUE)
rownames_df <- rownames_df[order(rownames_df$gene_id),] #make sure these are both in same order
mat <- mat[order(row.names(mat)), ] #note this means plot is in order of ENSG number
rownames(mat) = rownames_df$gene_name #change from ENSG to gene id
# and plot the actual heatmap
pheatmap(mat, annotation_col=df, show_rownames = T, fontsize_row = 5, legend = F, annotation_names_col = F, height = 10, 
         filename ="/results/GroupBvsGroupA_subset/GroupBvsGroupA_subset_heatmap.png")

##run GO enrichment, using res from above
# needed for goseq
txdb <- GenomicFeatures::makeTxDbFromGFF('gencode.v26.GRCh38.ERCC.genes.fix.gtf', format="gtf")
txsByGene=GenomicFeatures::transcriptsBy(txdb,"gene")
lengthData=median(width(txsByGene))
#DE = genes in GrpBvsGrpA that are SIGNIFICANTLY upregulated, and are upregulated in GroupBvsGroupD, Control is all genes assayed in GrpBvsGrpA
#have to remove version numbers from ENSG which are the rownames
rownames(res) <- gsub("\\..*","",rownames(res))
assayed.genes <- rownames(res)
fdr.threshold <- 0.05
##subset mat to only include genes that are significantly upregulated in GrpBvsGrpA and upregulated in GroupBvsGroupD
sig.genes <- res[!is.na(res$padj) & (res$padj < fdr.threshold) & rownames(res) %in% gsub("\\..*","",GrpBvsGrpA_subset$gene_id),]
de.genes <- rownames(sig.genes)[ which(sig.genes$padj < fdr.threshold) ]

gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

lengthDataDE <- lengthData[intersect(names(lengthData),names(gene.vector))]
gene.vectorDE <- gene.vector[intersect(names(gene.vector),names(lengthData))]

pwf=nullp(gene.vectorDE,"hg38","ensGene", bias.data = lengthDataDE)
GO.wall=goseq(pwf,"hg38","ensGene")
write.table(as.data.frame(GO.wall), file = "/results/GroupBvsGroupA_subset/GroupBvsGroupA_subset_GO.enrichment_results.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
