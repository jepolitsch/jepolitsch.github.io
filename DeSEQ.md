---
layout: alt
title: Julian Politsch
description: DeSEQ
---


# RNAseq workflow

### Step 1-Setup:

First we should make a fresh conda environment with only the nessasary packages in order to avoid future errors:
    
    $conda create -n deseq
    $conda install fastqc multiqc trimgalore salmon sra-tools=2.11
    $
    
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
Next  
