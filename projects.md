---
layout: alt
title: Julian Politsch
---
# RNAseq WorkFlow of Regulation of transcription elongation in response to osmostress:
## Julian Politsch, Bioinformatics 2021 La Sapienza



RNA-seq data comes from this paper: [Regulation of transcription elongation in response to osmostresslink:] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720810/). The data have been deposited in the Gene Expression Omnibus (GEO) database.


In linux terniaml we run the following SRA toolkit package commands:
  ### Retreive the RNA-seq Data
  ./prefetch SRR6765939 
  ### Convert from SRA to Fastq
  ./fastq-dump SRR6765939.sra
  

