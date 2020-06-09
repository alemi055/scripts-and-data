# TITRE
## **Audrée Lemieux, Alexandre J. Poulain, Stéphane Aris-Brosou**
### Department of Biology, University of Ottawa, Ottawa, ON K1N 6N5, Canada.
### Correspondence:


# Abstract


# Introduction


# Methods
## Sequencing and data processing
A first look at the quality of the raw data was done using FastQC (v0.11.8). Forward and reverse reads were trimmed for adapters, low-quality reads, unpaired reads and low quality bases at the end by using Trimmomatic (v0.36) (Bolger *et al.*, 2014) with the following parameters: phred33, ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:20, CROP:105, HEADCROP:15, AVGQUAL:20, MINLEN:36.

# Results


# Discussion


# References
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btu170

