---
layout: alt
title: Julian Politsch
description: projects
---
### DeSeq WorkFlow of Regulation of transcription elongation in response to osmostress:
#### Julian Politsch, Bioinformatics 2021 La Sapienza

RNA-seq data comes from this paper: [Regulation of transcription elongation in response to osmostresslink](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720810/). The data have been deposited in the Gene Expression Omnibus (GEO) database.

#### Abstract: 
Cells trigger massive changes in gene expression upon environmental fluctuations. The Hog1 stress-activated protein kinase (SAPK) is an important regulator of the transcriptional activation program that maximizes cell fitness when yeast cells are exposed to osmostress. Besides being associated with transcription factors bound at target promoters to stimulate transcriptional initiation, activated Hog1 behaves as a transcriptional elongation factor that is selective for stress-responsive genes. Here, we provide insights into how this signaling kinase functions in transcription elongation. Hog1 phosphorylates the Spt4 elongation factor at Thr42 and Ser43 and such phosphorylations are essential for the overall transcriptional response upon osmostress. The phosphorylation of Spt4 by Hog1 regulates RNA polymerase II processivity at stress-responsive genes, which is critical for cell survival under high osmostress conditions. Thus, the direct regulation of Spt4 upon environmental insults serves to stimulate RNA Pol II elongation efficiency.

#### High Level Overview:
1. Create and prepare our linux machine (This will be done locally but for more intensive jobs the workflow will be through a server
2. Download the RNAseq data from NCIB using SRA
3. Quality checking the reads with FastQC
4. Generate a comprehensive report using MultiQC
5. Data quantitification with Salmon
6. Differential expression analysis with DESeq2 in RStudio
7. Gene Ontology (GO) enrichment analysis to find the Up regulated genes
8. Differential Expression analysis of raw counts from paper authers
9. Quality checking the workflow


#### Step 1-Conda Enviroment Setup:

Lets start by creating a new conda environment with only the packages we need and their dependencies in order to avoid future errors:
```console
$ conda create -n deseq
$ conda activate deseq
$ conda install sra-tools=2.11 FastQC MultiQC Salmon
$ conda install -c bioconda trim-galore

```

SRA tools can also be installed with apt-get, but it is an older verson with less features.
    
#### Step 2-Download Data From SRA:

This study contains 12 RNA-seq samples from yeast and there are 3 replicates for 4 conditions.

First download a list of Run IDs in a file named SRR_lst.txt, and then download the SRA files with:
```
$cd Desktop/deseq/data
$cat SRR_Acc_list.txt | xargs -n 1 prefetch ; fastq-dump -0 home/julian/Desktop/deseq/data/fastq_dump *.sra --gzip
```
This step takes a lot of time and storage, each file is almost one GB its better to run it in the background by pressing "Ctrl + Z" followed by "$bg", then the progress is shown by entering "$jobs"


#### Step 3-Check reads quality:

In the same directory, containing the fastq.gz files we run the following command in order to check the quality of our reads:
```
$ cd ~/Desktop/deseq/data/fastq_dump
$ fastqc -o ~/Desktop/deseq/data/fastqc_reports *.fastq.gz
```    

Comparing out samples to the fastqc [good illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html#M8) and [bad illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) we can see a few flags in our data:

##### Duplicate Sequences:
In a diverse library most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication is more likely to indicate some kind of enrichment bias (eg PCR over amplification). 

Since we used rnaSEQ of osmostressed cells this warning message can be ignored for now as we would expect some sequences to be over represented considering the data type, and likely does not show a bias towards an enrichment type.

In order to be able to observe lowly expressed transcripts it is therefore common to greatly over-sequence high expressed transcripts, and this will potentially create large set of duplicates. This will result in high overall duplication in this test, and will often produce peaks in the higher duplication bins.

##### Per base sequence content: 
In a random library you would expect that there would be little to no difference between the different bases of a sequence run, so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in your genome, but in any case they should not be hugely imbalanced from each other.

It's worth noting that some types of library will always produce biased sequence composition, normally at the start of the read. Libraries produced by priming using random hexamers (including nearly all RNA-Seq libraries) and those which were fragmented using transposes inherit an intrinsic bias in the positions at which reads start. This bias does not concern an absolute sequence, but instead provides enrichment of a number of different K-mers at the 5' end of the reads. Whilst this is a true technical bias, it isn't something which can be corrected by trimming and in most cases doesn't seem to adversely affect the downstream analysis. It will however produce a warning or error in this module.

#### Step 4-Compare all fastQC reports:

MultiQC is a tool to aggregate results from bioinformatics analysis logs across multiple samples into a single report, perfect for summarizing the output of numerous bioinformatics tools.

Simply run with:
```
$ multiqc -i multiqc_report_untrimmed -o ~/Desktop/deseq/data .
```
While in the directory you want the tool to search

#### Step 5-Quality and adaptor trimming with Trim Galore

Although the official documentation says the per base sequence content warning can be safely ignored, we should still try to trim the data in order to improve the quality of the model.

```
$conda install -c bioconda trim-galore
$for fastq in *.fastq.gz
$do 
>trim_galore --quality 25 --stringency 5 --length 35 --fastqc $fastq
$done
```
After running the multiqc report on the trimmed samples, we see the same error flags as the untrimmed data, although around 3 million (0.5%) of the base pairs were able to be trimmed. We do see a marked decrease in low base quality, and removing primers is always a good idea.

We will move on with the trimmed files

#### Step 6-Quantifying of samples with salmon

First we need to donwload the reference transcriptome of S. Cerevisiae [Here](http://ftp.ensembl.org/pub/release-104/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz)

Next we create an index for salmon using:
```
salmon index -t Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i S_cere_index
```

And then we can write a quick recursive nano script to run all the samples through salmon in one command

```
for file in *.fq.gz
do
salmon quant -i S_cere_index -l A -r ${file} --validateMappings -o ${file}_quant
done
```
and in the same directory we run
```
$sh [sh file name].sh
```

and the result will be 12 nicely formatted quantified expression files which we can take to deseq analysis in R

#### Step 7-Differential Expression Analysis in R:

First install all the packages:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

packages <- c("tximport", "GenomicFeatures", "DESeq2", "clusterProfiler", "AnnotationDbi", "org.Sc.sgd.db", "org.Mm.eg.db")
BiocManager::install(packages)

install.packages(c("readr", "ggplot2"))
```
Next I will show the step-by-step of my R code, but dont get bogged down on the details and skip to Step 8 as I will attach the scripts in a seperate file

##### Setup
```{r}
files <- file.path("/Users/julian/Desktop/deseq/results", colData$directory, "quant.sf")
names(files) <- colData$sample
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(readr)
library(ggplot2)
library(stringi)
```
##### Initial Statistics for each quante file
```{r counts_summary} 
quant <- 
counts <- read.table("/Users/julian/Desktop/deseq/quants/quant.sf", header = TRUE, row.names = 1)
summary(counts)
```

##### Data import and quantification from salmon
```{r}
# salmon results directory ----
setwd("/Users/julian/Desktop/deseq")
dir("/Users/julian/Desktop/deseq")

# import sample Metadata ----
colData <- 
  read.csv("/Users/julian/Desktop/deseq/sampledata.csv", header = FALSE)
colnames(colData) <- c("fastq", "sample")

colData$group <- rep(c("wt_nostress", "wt_stress", "mut_nostress", "mut_stress"), each = 3)
colData$directory <- dir("/Users/julian/Desktop/deseq/results/")

rownames(colData) <- colData$sample

# quant.sf file paths ----
files <- file.path("/Users/julian/Desktop/deseq/results", colData$directory, "quant.sf")
names(files) <- colData$sample

txi.inf.rep <- tximport::tximport(files, type = "salmon", txOut = TRUE)
names(txi.inf.rep$sample)
```

```{r}
# create tx2gene object ----
txdb <- GenomicFeatures::makeTxDbFromGFF("/Users/julian/Desktop/deseq/Saccharomyces_cerevisiae.R64-1-1.104.gtf")
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

# check if sample names match between objects
identical(x = rownames(colData), y = colnames(txi.salmon$counts))
```

##### Experimental Setup and formula design 
```{r deseq_setup_1}
colData <- read.table("/Users/julian/Desktop/deseq/sampledata.csv", header = FALSE, sep = ',', stringsAsFactors = TRUE)
colData
designFormula <- "~ group"
colData$group <- rep(c("wt_nostress", "wt_stress", "mut_nostress", "mut_stress"), each = 3)
colData$directory <- dir("/Users/julian/Desktop/deseq/results/")
```

##### Filtering non-expressed genes
```{r} 
dds <- DESeqDataSetFromTximport(txi = txi.salmon, colData = colData, design = ~group)
dds$group <- relevel(x = dds$group, ref = "wt_nostress")
dds <- DESeq(dds)

is_expressed <- assay(dds) >= 1
head(is_expressed)
keep <- rowSums(counts(dds)) > 10
table(keep)
dds <- dds[keep, ]
nrow(dds)
```

### DeSEQ pipeline
```{r}
dds <- DESeq(dds)
head(dds)
res <- results(dds)
summary(res)

ws_wns <- results(dds, contrast = c("group", "wt_stress" ,"wt_nostress"))
ms_mns <- results(dds, contrast = c("group", "mut_stress", "mut_nostress"))
ms_ws <- results(dds, contrast = c("group", "mut_stress", "wt_stress"))
mns_wns <- results(dds, contrast = c("group", "mut_nostress", "wt_nostress"))
#### MA Plots
# - An **MA plot** is useful to observe if the data normalization worked well. The **MA plot** is a scatter plot where the x-axis denotes the average of normalized counts across samples and the y-axis denotes the log fold change in the given contrast. Most points are expected to be on the horizontal 0 line (most genes are not expected to be differentially expressed).

ws_wns <- ws_wns[order(ws_wns$pvalue), ]
DESeq2::plotMA(ws_wns, alpha = 0.05)

ms_mns <- ms_mns[order(ms_mns$pvalue), ]
DESeq2::plotMA(ms_mns, alpha = 0.05)

ms_ws <- ms_ws[order(ms_ws$pvalue), ]
DESeq2::plotMA(ms_ws, alpha = 0.05)

mns_wns <- mns_wns[order(mns_wns$pvalue), ]
DESeq2::plotMA(mns_wns, alpha = 0.05)
```

##### Sample Distances and PCA plots

A final diagnosis is to check the biological reproducibility of the sample replicates in a PCA plot or a heatmap. 

To plot the PCA results, we need to extract the **normalized counts** from the DESeqDataSet object. 

It is possible to color the points in the scatter plot by the variable of interest, which helps to see if the replicates cluster well.

```{r}
#### Sample distances
library("pheatmap")
library("RColorBrewer")

#### PCA ----
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
  theme_bw()
# Sample Dist.
sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
```

```{r}
library(DESeq2)
DESeq2::plotMA(object = dds, ylim = c(-5, 5))
```

```{r}
library(ggplot2)
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
   theme_bw()
```

##### finding DE genes
```{r}
#remove genes with NA values 
de1 <- ws_wns[!is.na(ws_wns$padj),]
#select genes with adjusted p-values below 0.1
de1 <- de1[de1$padj < 0.1,]
#select genes with absolute log2 fold change above 1 (two-fold change)
de1 <- de1[abs(de1$log2FoldChange) > 1,]

de2 <- ms_mns[!is.na(ms_mns$padj),]
de2 <- de2[de2$padj < 0.1,]
de2 <- de2[abs(de2$log2FoldChange) > 1,]

de3 <- ms_ws[!is.na(ms_ws$padj),]
de3 <- de3[de3$padj < 0.1,]
de3 <- de3[abs(de3$log2FoldChange) > 1,]

de4 <- mns_wns[!is.na(mns_wns$padj),]
de4 <- de4[de4$padj < 0.1,]
de4 <- de4[abs(de4$log2FoldChange) > 1,]

de1
de2
de3
de4
```
