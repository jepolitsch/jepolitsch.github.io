# load packages ----
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(readr)
library(ggplot2)

# salmon results directory ----
setwd("/home/julian/Desktop/deseq/")
dir("/home/julian/Desktop/deseq/")

# import sample Metadata ----
colData <- 
  read.csv("sampledata.csv", header = FALSE)
colnames(colData) <- c("fastq", "sample")

colData$group <- rep(c("wt_nostess", "wt_stress", "mut_nostress", "mut_stress"), each = 3)
colData$directory <- dir("/salmon_deseq/")

rownames(colData) <- colData$sample

# quant.sf file paths ----
files <- file.path("/salmon_deseq/", colData$directory, "quant.sf")
names(files) <- colData$sample

# create tx2gene object ----
txdb <- GenomicFeatures::makeTxDbFromGFF("/Saccharomyces_cerevisiae.R64-1-1.100.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

# reorder columns of tx2gene ----
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
# check tx2gene
head(tx2gene)

# import salmon quantification data ----
txi.salmon <- tximport(files = files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# check imported data ----
head(txi.salmon$counts)

# check if sample names match between objects ----
identical(x = rownames(colData), y = colnames(txi.salmon$counts))


# DESeq2 pipeline ----
dds <- DESeqDataSetFromTximport(txi = txi.salmon, colData = colData, design = ~group)

# filtering not expressed genes ----
keep <- rowSums(counts(dds)) > 1
table(keep)
dds <- dds[keep, ]
nrow(dds)

dds$group <- relevel(x = dds$group, ref = "WT")
dds <- DESeq(dds)

# save the results of DE analysis ----
res <- results(dds)

# look at the results ----
summary(res)
res

# Visualization of results: MA plot ----
DESeq2::plotMA(res, alpha = 0.05)

# PCA ----
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
  theme_bw()

# Sample distances ----
library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


#remove genes with NA values 
DE <- res[!is.na(res$padj),]
#select genes with adjusted p-values below 0.1
DE <- DE[DE$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

# export data to file ----
write.table(x = res_005, file = "results_yeast.txt", sep = "\t", col.names = NA)