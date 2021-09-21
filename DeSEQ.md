---
layout: alt
title: Julian Politsch
description: DeSEQ
---

RNA-seq data comes from this paper: [Regulation of transcription elongation in response to osmostresslink](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720810/). The data have been deposited in the Gene Expression Omnibus (GEO) database.
10

11
#### Abstract: 
12
Cells trigger massive changes in gene expression upon environmental fluctuations. The Hog1 stress-activated protein kinase (SAPK) is an important regulator of the transcriptional activation program that maximizes cell fitness when yeast cells are exposed to osmostress. Besides being associated with transcription factors bound at target promoters to stimulate transcriptional initiation, activated Hog1 behaves as a transcriptional elongation factor that is selective for stress-responsive genes. Here, we provide insights into how this signaling kinase functions in transcription elongation. Hog1 phosphorylates the Spt4 elongation factor at Thr42 and Ser43 and such phosphorylations are essential for the overall transcriptional response upon osmostress. The phosphorylation of Spt4 by Hog1 regulates RNA polymerase II processivity at stress-responsive genes, which is critical for cell survival under high osmostress conditions. Thus, the direct regulation of Spt4 upon environmental insults serves to stimulate RNA Pol II elongation efficiency.
13

14
#### High Level Overview:
15
1. Create and prepare our linux machine (This will be done locally but for more intensive jobs the workflow will be through a server
16
2. Download the RNAseq data from NCIB using SRA
17
3. Quality checking the reads with FastQC
18
4. Generate a comprehensive report using MultiQC
19
5. Data quantitification with Salmon
20
6. Differential expression analysis with DESeq2 in RStudio
21
7. Gene Ontology (GO) enrichment analysis to find the Up regulated genes
22
8. Differential Expression analysis of raw counts from paper authers
23
9. Quality checking the workflow
@jepolitsch
Update DeSEQ.md
last month
24

### Step 1-Setup:

First we should make a fresh conda environment with only the nessasary packages in order to avoid future errors:
```  
    $conda create -n deseq
    $conda install fastqc multiqc trimgalore salmon sra-tools=2.11
``` 
#### Step 2-Download Data From SRA:

This study contains 12 RNA-seq samples from yeast and there are 3 replicates for 4 conditions.

First download a list of Run IDs in a file named [SRR_lst.txt](LINK TO HTTP OF LIST), and then download the SRA files with:
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

# Rstudio Portion
## Step 7-Differential Expression Analysis in R:

First install all the packages:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

packages <- c("tximport", "GenomicFeatures", "DESeq2", "clusterProfiler", "AnnotationDbi", "org.Sc.sgd.db", "org.Mm.eg.db")
BiocManager::install(packages)

install.packages(c("readr", "ggplot2"))
```
Next

### Setup
```{r Setup, message=FALSE, warning=FALSE, results='hide'}
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(readr)
library(ggplot2)
library(stringi)
library("AnnotationDbi")
library("org.Sc.sgd.db")
library(DOSE)
library(tximport)
library(clusterProfiler)
```

### Initial Statistics for each quant file
```{r counts_summary} 
counts <- read.table("/Users/julian/Desktop/deseq/quants/quant.sf", header = TRUE, row.names = 1)
summary(counts)
```

### Load sample Metadata
```{r Load Meta Data}
colData <-  read.csv("/Users/julian/Desktop/deseq/sampledata.csv", header = FALSE)
colnames(colData) <- c("fastq", "sample")
files <- file.path("/Users/julian/Desktop/deseq/results", colData$directory, "quant.sf")
```

### Data import and quantification from Salmon
Salmon is a tool for transcript quantification salmon has two phases, indexing to a reference transcriptome, and quantification

```{r Salmon Quantification, message=FALSE, warning=FALSE, results='hide'}
#Labeling
colData$group <- rep(c("wt_nostress", "wt_stress", "mut_nostress", "mut_stress"), each = 3)
colData$directory <- dir("/Users/julian/Desktop/deseq/results/")
rownames(colData) <- colData$sample

#quant.sf file paths containing the quantification results
files <- file.path("/Users/julian/Desktop/deseq/results", colData$directory, "quant.sf")
names(files) <- colData$sample
txi.inf.rep <- tximport::tximport(files, type = "salmon", txOut = TRUE)
```

### Create txtb object 
Object which maps the 5' and 3' UTR and maps the coding sequences, exons, to their associated genome, in this case S. Cerevisiae
```{r TXTB Object, message=FALSE, warning=FALSE, results='hide'}
txdb <- GenomicFeatures::makeTxDbFromGFF("/Users/julian/Desktop/deseq/Saccharomyces_cerevisiae.R64-1-1.104.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

#Reorder columns of tx2gene
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]

#Check tx2gene
head(tx2gene)

#Import salmon quantification data
txi.salmon <- tximport(files = files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# check imported data ----
head(txi.salmon$counts)
```

### Checking quantification correctness
```{r Check correctness}
# check if sample names match between objects
identical(x = rownames(colData), y = colnames(txi.salmon$counts))
```

### Filtering non-expressed genes
```{r Filtering NE Genes, message=FALSE, warning=FALSE} 
dds <- DESeqDataSetFromTximport(txi = txi.salmon, colData = colData, design = ~group)
dds <- DESeq(dds)

keep <- rowSums(counts(dds)) > 1
table(keep)
dds <- dds[keep, ]
dds$group <- relevel(x = dds$group, ref = "wt_nostress")
```

### DeSEQ pipeline
```{r DESEQ, message=FALSE, warning=FALSE, results='hide'}
dds <- DESeq(dds)
res <- results(dds)

ws_wns <- results(dds, contrast = c("group", "wt_stress" ,"wt_nostress"))
ms_mns <- results(dds, contrast = c("group", "mut_stress", "mut_nostress"))
ms_ws <- results(dds, contrast = c("group", "mut_stress", "wt_stress"))
mns_wns <- results(dds, contrast = c("group", "mut_nostress", "wt_nostress"))
```
### Finding DE genes
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

# Now these objects are only those genes with padj < 0.1 and log2fold change >1 as: de1, de2, de3, de4
```

### MA Plots
A **MA plot** is useful to observe if the data normalization worked well. The **MA plot** is a scatter plot where the x-axis denotes the average of normalized counts across samples and the y-axis denotes the log fold change in the given contrast. Most points are expected to be on the horizontal 0 line (most genes are not expected to be differentially expressed).
```{r MA Plots}

ws_wns <- ws_wns[order(ws_wns$pvalue), ]
DESeq2::plotMA(de1, title('WS vs WNS'), alpha = 0.05)

ms_mns <- ms_mns[order(ms_mns$pvalue), ]
DESeq2::plotMA(de2, title('MS vs MNS'),alpha = 0.05)

ms_ws <- ms_ws[order(ms_ws$pvalue), ]
DESeq2::plotMA(de3,title('MS vs WS'), alpha = 0.05)

mns_wns <- mns_wns[order(mns_wns$pvalue), ]
DESeq2::plotMA(de4, title('MNS vs WNS'), alpha = 0.05)
```

### Sample Distances and PCA plots

A final diagnosis is to check the biological reproducibility of the sample replicates in a PCA plot or a heatmap. 

To plot the PCA results, we need to extract the **normalized counts** from the DESeqDataSet object. 

It is possible to color the points in the scatter plot by the variable of interest, which helps to see if the replicates cluster well.

```{r PCA and Sample Distances}
# Pkgs
library("pheatmap")
library("RColorBrewer")

# PCA Plots
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
  theme_bw()

# Sample Dist.
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


```

### Go enrichment analysis (Gene Ontology)
Gene ontology aims to cluster groups of genes that are up or down regulated based on cellular role or specific pathway. Here we only want to look at those upregulated

##### WS vs WNS
```{r WSvsWNS, message=FALSE, warning=FALSE, results='hide'}
de1go <- de1[de1$padj < 0.05,]
de1go <- de1go[(de1go$log2FoldChange) > 1,]

de1go$symbol <- mapIds(org.Sc.sgd.db,
                     keys=rownames(de1go),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

go1 <- clusterProfiler::enrichGO(row.names(de1go), "org.Sc.sgd.db", keyType = "ENSEMBL", ont = "BP")
clusterProfiler::buildGOmap(go1)
```

##### MS vs MNS
```{r MSvsMNS, message=FALSE, warning=FALSE, results='hide'}
de2go <- de2[de2$padj < 0.05,]
de2go <- de2go[(de2go$log2FoldChange) > 1,]
de2go$symbol <- mapIds(org.Sc.sgd.db,
                     keys=rownames(de2go),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")


go2 <- clusterProfiler::enrichGO(row.names(de2go), "org.Sc.sgd.db", keyType = "ENSEMBL", ont = "BP")
clusterProfiler::buildGOmap(go2)
```

##### MS vs WS
```{r MSvsWS, message=FALSE, warning=FALSE, results='hide'}
de3go <- de3[de3$padj < 0.05,]
de3go <- de3go[(de3go$log2FoldChange) > 1,]
de3go$symbol <- mapIds(org.Sc.sgd.db,
                     keys=rownames(de3go),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")


go3 <- clusterProfiler::enrichGO(row.names(de3go), "org.Sc.sgd.db", keyType = "ENSEMBL", ont = "BP")
clusterProfiler::buildGOmap(go3)
```

##### MNS vs WNS
```{r MNSvsWNS, message=FALSE, warning=FALSE, results='hide'}
de4go <- de4[de4$padj < 0.05,]
de4go <- de4go[(de4go$log2FoldChange) > 1,]
de4go$symbol <- mapIds(org.Sc.sgd.db,
                     keys=rownames(de4go),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")


go4 <- clusterProfiler::enrichGO(row.names(de4go), "org.Sc.sgd.db", keyType = "ENSEMBL", ont = "BP")
clusterProfiler::buildGOmap(go4)
```

#### Results
```{r fig.align='left', fig.width=12, error= FALSE, warning=FALSE, message=FALSE}
library(ggplot2)
library(DOSE)
dotplot(go1, showCategory=15, title = "WS vs WNS Up regulated GO pathways")
dotplot(go2, showCategory=15, title = "MS vs MNS Up regulated GO pathways")
dotplot(go3, showCategory=15, title = "MS vs WS Up regulated GO pathways")
dotplot(go4, showCategory=15, title = "MNS vs WNS Up regulated GO pathways")

#since we are looking for only upregulated genes the new MAplots will look like:
de1ex <- de1[(de1$log2FoldChange) > 1,]
DESeq2::plotMA(de1ex, title('WS vs WNS UP only'), alpha = 0.05)
```
## Step 7

```{r Step 7, , message=FALSE, warning=FALSE}
counts <- read_tsv("RawCounts.tsv")
colnames(counts)[1] <- "GENE_ID"
colnames(counts)[2:13] <- rownames(colData)
countData <- as.matrix(subset(counts, select = c(-GENE_ID)))
summary(countData)
# DESeq2 pipeline 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~group)
# filtering not expressed genes 
keep <- rowSums(counts(dds)) > 1
dds <- DESeq(dds)
# PCA 
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
  theme_bw()
# save the results of DE analysis 
WSvsWN2 <- results(dds, contrast = c("group", "wt_stress" ,"wt_nostress"))
MSvsMN2 <- results(dds, contrast = c("group", "mut_stress", "mut_nostress"))
MSvsWS2 <- results(dds, contrast = c("group", "mut_stress", "wt_stress"))
MNvsWN2 <- results(dds, contrast = c("group", "mut_nostress", "wt_nostress"))

WSvsWN2 <- WSvsWN2[!is.na(WSvsWN2$padj),]
WSvsWN2 <- WSvsWN2[WSvsWN2$padj < 0.05,]
WSvsWN2 <- WSvsWN2[abs(WSvsWN2$log2FoldChange) > 1,]
DESeq2::plotMA(WSvsWN2, alpha = 0.05)

MSvsMN2 <- MSvsMN2[!is.na(MSvsMN2$padj),]
MSvsMN2 <- MSvsMN2[MSvsMN2$padj < 0.05,]
MSvsMN2 <- MSvsMN2[abs(MSvsMN2$log2FoldChange) > 1,]
DESeq2::plotMA(MSvsMN2, alpha = 0.05)

MSvsWS2 <- MSvsWS2[!is.na(MSvsWS2$padj),]
MSvsWS2 <- MSvsWS2[MSvsWS2$padj < 0.05,]
MSvsWS2 <- MSvsWS2[abs(MSvsWS2$log2FoldChange) > 1,]
DESeq2::plotMA(MSvsWS2, alpha = 0.05)

MNvsWN2 <- MNvsWN2[!is.na(MNvsWN2$padj),]
MNvsWN2 <- MNvsWN2[MNvsWN2$padj < 0.05,]
MNvsWN2 <- MNvsWN2[abs(MNvsWN2$log2FoldChange) > 1,]
DESeq2::plotMA(MNvsWN2, alpha = 0.05)
```

### Q1 and Q2
How many genes are differentially expressed at the thresholds of padj < 0.05? And padj < 0.1?
```{r Q1}
de1q1 <- ws_wns[!is.na(ws_wns$padj),]
de1q1 <- de1q1[de1q1$padj < 0.1,]
de1q1 <- de1q1[abs(de1q1$log2FoldChange) > 1,]
#thresholds of padj < 0.1
nrow(de1q1)
#thresholds of padj < 0.05
de1q1 <- de1q1[de1q1$padj < 0.05,]
nrow(de1q1)

de2q1 <- ms_mns[!is.na(ms_mns$padj),]
de2q1 <- de2q1[de2q1$padj < 0.1,]
de2q1 <- de2q1[abs(de2q1$log2FoldChange) > 1,]

#thresholds of padj < 0.1
nrow(de2q1)
#thresholds of padj < 0.05
de2q1 <- de2q1[de2q1$padj < 0.05,]
nrow(de2q1)

de3q1 <- ms_ws[!is.na(ms_ws$padj),]
de3q1 <- de3q1[de3q1$padj < 0.1,]
de3q1 <- de3q1[abs(de3q1$log2FoldChange) > 1,]

#thresholds of padj < 0.1
nrow(de3q1)
#thresholds of padj < 0.05
de3q1 <- de3q1[de3q1$padj < 0.05,]
nrow(de3q1)

de4q1 <- mns_wns[!is.na(mns_wns$padj),]
de4q1 <- de4q1[de4$padj < 0.1,]
de4q1 <- de4q1[abs(de4q1$log2FoldChange) > 1,]

#thresholds of padj < 0.1
nrow(de4q1)
#thresholds of padj < 0.05
de4q1 <- de4q1[de4q1$padj < 0.05,]
nrow(de4q1)
```
### Q3
How many genes at the thresholds of padj < 0.05 are upregulated ( > 1)?

```{r Q2}
de1q1 <- de1q1[(de1q1$log2FoldChange) > 1,]
de2q1 <- de2q1[(de2q1$log2FoldChange) > 1,]
de3q1 <- de3q1[(de3q1$log2FoldChange) > 1,]
de4q1 <- de4q1[(de4q1$log2FoldChange) > 1,]
nrow(de1q1)
nrow(de2q1)
nrow(de3q1)
nrow(de4q1)
```
### Q4
Choice one of the GO enrichment results and report how many categories are significant (p.adjust < 0.05). Hint: first convert the object with GO results with as.data.frame function.
```{r Q3}
go1DataFrame <- as.data.frame(go1)
nrow(go1DataFrame)
```

How many genes are present in the most enriched category of the GO enrichment result?

```{r}
max(go1DataFrame$Count)
```
Which corresponds to oxidation reduction processes (logical)

### Q5
Are the numbers of significant differentially expressed genes of Step 7 the same as the results in Step 5?

```{r}
nrow(WSvsWN2)
nrow(de1)
```
The number of significant dif. expressed genes in step 7 and step 5 is very similar, we can account for trimming as the margin of error, and the results from step 5 are likely more accurate.

### Questions?
