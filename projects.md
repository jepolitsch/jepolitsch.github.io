---
layout: alt
title: Julian Politsch
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

#### Step 1: 
Lets start by creating a new conda environment with only the packages we need and their dependencies in order to avoid future errors:
```console
$ conda create -n deseq
$ conda activate deseq
$ conda install sra-tools=2.11 FastQC MultiQC Salmon
$ conda install -c bioconda trim-galore

```

SRA tools can also be installed with apt get, but it is an older verson with less feature 
In linux terniaml we run the following SRA toolkit package commands:
    
### Step 2-Download Data From SRA:

This study contains 12 RNA-seq samples from yeast and there are 3 replicates for 4 conditions.

First download a list of Run IDs in a file named [SRR_lst.txt](LINK TO HTTP OF LIST), and then download the SRA files with:
```
$cd Desktop/deseq/data
$cat SRR_Acc_list.txt | xargs -n 1 prefetch ; fastq-dump -0 home/julian/Desktop/deseq/data/fastq_dump *.sra --gzip
```

This step takes a lot of time and storage, each file is almost one GB its better to run it in the background by pressing "Ctrl + Z" followed by "$bg", then the progress is shown by entering "$jobs"


### Step 3-Check reads quality:

In the same directory, containing the fastq.gz files we run the following command in order to check the quality of our reads:
```    
    $cd ~/Desktop/deseq/data/fastq_dump
    $fastqc -o ~/Desktop/deseq/data/fastqc_reports *.fastq.gz
```

### Step 4-Compare all fastQC reports:

MultiQC is a tool to aggregate results from bioinformatics analysis logs across multiple samples into a single report, perfect for summarizing the output of numerous bioinformatics tools.

Simply run with:
    $multiqc .
While in the directory you want the tool to search    

