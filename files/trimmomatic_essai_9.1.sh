#!/bin/bash
#SBATCH -c 5                               # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=50gb                         # mem in gb
#SBATCH -t 0-10:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
#SBATCH -J essai_9.1                 	# Name of job

module load trimmomatic

cd /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
HI.4444.003.Index_1.GR_RNA_AS3_R1_fastq.gz HI.4444.003.Index_1.GR_RNA_AS3_R2_fastq.gz \
003.Index_1.GR_RNA_AS3_R1_paired_fastq.gz 003.Index_1.GR_RNA_AS3_R1_unpaired_fastq.gz \
003.Index_1.GR_RNA_AS3_R2_paired_fastq.gz 003.Index_1.GR_RNA_AS3_R2_unpaired_fastq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
HI.4444.003.Index_3.GR_RNA_BS3-3_R1_fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2_fastq.gz \
003.Index_3.GR_RNA_BS3-3_R1_paired_fastq.gz 003.Index_3.GR_RNA_BS3-3_R1_unpaired_fastq.gz \
003.Index_3.GR_RNA_BS3-3_R2_paired_fastq.gz 003.Index_3.GR_RNA_BS3-3_R2_unpaired_fastq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
HI.4444.003.Index_5.GR_RNA_S3-3_R1_fastq.gz HI.4444.003.Index_5.GR_RNA_S3-3_R2_fastq.gz \
003.Index_5.GR_RNA_S3-3_R1_paired_fastq.gz 003.Index_5.GR_RNA_S3-3_R1_unpaired_fastq.gz \
003.Index_5.GR_RNA_S3-3_R2_paired_fastq.gz 003.Index_5.GR_RNA_S3-3_R2_unpaired_fastq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
HI.4444.003.Index_9.GR_RNA_RS3-1_R1_fastq.gz HI.4444.003.Index_9.GR_RNA_RS3-1_R2_fastq.gz \
003.Index_9.GR_RNA_RS3-1_R1_paired_fastq.gz 003.Index_9.GR_RNA_RS3-1_R1_unpaired_fastq.gz \
003.Index_9.GR_RNA_RS3-1_R2_paired_fastq.gz 003.Index_9.GR_RNA_RS3-1_R2_unpaired_fastq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
HI.4444.003.Index_12.GR_RNA_S4-3_R1_fastq.gz HI.4444.003.Index_12.GR_RNA_S4-3_R2_fastq.gz \
003.Index_12.GR_RNA_S4-3_R1_paired_fastq.gz 003.Index_12.GR_RNA_S4-3_R1_unpaired_fastq.gz \
003.Index_12.GR_RNA_S4-3_R2_paired_fastq.gz 003.Index_12.GR_RNA_S4-3_R2_unpaired_fastq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
HI.4444.003.Index_19.GR_RNA_S5-3_R1_fastq.gz HI.4444.003.Index_19.GR_RNA_S5-3_R2_fastq.gz \
003.Index_19.GR_RNA_S5-3_R1_paired_fastq.gz 003.Index_19.GR_RNA_S5-3_R1_unpaired_fastq.gz \
003.Index_19.GR_RNA_S5-3_R2_paired_fastq.gz 003.Index_19.GR_RNA_S5-3_R2_unpaired_fastq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/003* /home/alemi055/scratch/ete2020/RNA_Arctic/quality_control_trimmed

rm /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/003*

cd /home/alemi055/scratch/ete2020/RNA_Arctic/quality_control_trimmed

module load fastqc

fastqc *_paired_*
