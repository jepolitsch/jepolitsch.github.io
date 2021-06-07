---
layout: alt
title: Julian Politsch
---
## DeSeq WorkFlow of Regulation of transcription elongation in response to osmostress:
### Julian Politsch, Bioinformatics 2021 La Sapienza

RNA-seq data comes from this paper: [Regulation of transcription elongation in response to osmostresslink.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720810/). The data have been deposited in the Gene Expression Omnibus (GEO) database.

#### Abstract: 
Cells trigger massive changes in gene expression upon environmental fluctuations. The Hog1 stress-activated protein kinase (SAPK) is an important regulator of the transcriptional activation program that maximizes cell fitness when yeast cells are exposed to osmostress. Besides being associated with transcription factors bound at target promoters to stimulate transcriptional initiation, activated Hog1 behaves as a transcriptional elongation factor that is selective for stress-responsive genes. Here, we provide insights into how this signaling kinase functions in transcription elongation. Hog1 phosphorylates the Spt4 elongation factor at Thr42 and Ser43 and such phosphorylations are essential for the overall transcriptional response upon osmostress. The phosphorylation of Spt4 by Hog1 regulates RNA polymerase II processivity at stress-responsive genes, which is critical for cell survival under high osmostress conditions. Thus, the direct regulation of Spt4 upon environmental insults serves to stimulate RNA Pol II elongation efficiency.

#### High Level Overview:
1. Create and prepare our linux machine (This will be done locally but for more intensive jobs the workflow will be through a server
2. Download the RNAseq data from NCIB using SRA
3. Quality checking the reads with FastQC
4. Generate a comprehensive report using MultiQC
5. Data quantitification with Salmon
6. Differential expression analysis with DESeq2 in RStudio
7. 
```
$conda create -n deseq
```

In linux terniaml we run the following SRA toolkit package commands:
  ### Retreive the RNA-seq Data
  ./prefetch SRR6765939 
  ### Convert from SRA to Fastq
  ./fastq-dump SRR6765939.sra
  

