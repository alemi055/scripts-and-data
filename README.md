# RNA ARCTIC DATA
Pipeline developed on pair of files "HI.4444.003.Index_3.GR_RNA_BS3-3_R1_fastq.gz" & "HI.4444.003.Index_3.GR_RNA_BS3-3_R2_fastq.gz"

## 1. FastQC

    #!/bin/bash
    #SBATCH -c 5                            # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=50G                       # mem in gb
    #SBATCH -t 0-5:0:0                      # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J fastqc                       # Name of job

    module load fastqc

    fastqc -t 5 *.fastq.gz

    cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*html /home/alemi055/scratch/ete2020/RNA_Arctic/fastqc_1/

    cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*zip /home/alemi055/scratch/ete2020/RNA_Arctic/fastqc_1/

    cd /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data

    rm *html

    rm *zip
    
## 2. Trimmomatic

    #!/bin/bash
    #SBATCH -c 6                               # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=50gb                         # mem in gb
    #SBATCH -t 1-0:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J trimmmomatic                    # Name of job

    module load trimmomatic/0.36

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

    cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*_paired.fastq.gz /home/alemi055/scratch/ete2020/RNA_Arctic/paired_trimmed

    cd /home/alemi055/scratch/ete2020/RNA_Arctic/paired_trimmed

    module load fastqc/0.11.8

    fastqc -t 5 003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz

## 3. FastQC

    #!/bin/bash
    #SBATCH -c 5                            # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=50G                       # mem in gb
    #SBATCH -t 0-5:0:0                      # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J fastqc                       # Name of job

    module load fastqc

    fastqc -t 5 *_paired.fastq.gz

    cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*html /home/alemi055/scratch/ete2020/RNA_Arctic/fastqc_1/

    cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*zip /home/alemi055/scratch/ete2020/RNA_Arctic/fastqc_1/

    cd /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data

    rm *html

    rm *zip
    
## 3. Assembly de novo

#### Trinity

    #!/bin/bash
    #SBATCH -c 10                              # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=160G                         # mem in gb
    #SBATCH -t 13-23:0:0                       # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J trinity                 		   # Name of job

    module load gcc/7.3.0 openmpi/3.1.4 samtools jellyfish salmon trinity/2.9.0
    Trinity --seqType fq --max_memory 120G --CPU 10 --left 003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz --right 003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz

#### SPAdes

    #!/bin/bash
    #SBATCH -c 6                               # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=120G                         # mem in gb
    #SBATCH -t 6-23:0:0                        # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J spades                  	       # Name of job

    module load gcc/7.3.0 spades/3.13.1

    spades.py --rna -1 004.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz -2 004.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz -o /home/alemi055/scratch/ete2020/RNA_Arctic/paired_trimmed/004_spades
    
### 3.1 Assembly - QC

### Trinity



### SPAdes

1. contig_sizes_spades.R

2. split_spades.R
