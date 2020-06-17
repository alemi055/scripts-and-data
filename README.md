# RNA ARCTIC DATA

## 1. Download raw data
    read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://genomequebec.mcgill.ca/nanuqMPS/readsetList?projectId=15225&tech=HiSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

We will need to enter our Nanuq username and password.


## 2. FastQC

    #1. Download the fastqc_script.sh to a directory. 
    scp /drives/c/Users/Audrée/Downloads/fastqc_script.sh alemi055@cedar.computecanada.ca:~/scratch/
    
    #2. Move to the directory.
    cd /home/alemi055/scratch/
    
    #3. Convert the Windows file to a Unix file
    dos2unix fastqc_script.sh
    
    #4. Submit the script as a batch job
    sbatch fastqc_script.sh


## 3. Quality Control with Trimmomatic
   
Code - Try #10

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 CROP:105 HEADCROP:15 AVGQUAL:20 MINLEN:36

GENERIC CODE

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.00!.Index_!.GR_RNA_!_R1_fastq.gz HI.4444.00!.Index_!.GR_RNA_!_R2_fastq.gz \
    00!.Index_!.GR_RNA_!_R1_paired_fastq.gz 00!.Index_!.GR_RNA_!_R1_unpaired_fastq.gz \
    00!.Index_!.GR_RNA_!_R2_paired_fastq.gz 00!.Index_!.GR_RNA_!_R2_unpaired_fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 CROP:105 HEADCROP:15 AVGQUAL:20 MINLEN:36
    
*Replace "!" by appropriate number/letter*


## 4. De novo Assembly

### Trinity

Code - Try #4

    #!/bin/bash
    #SBATCH -c 10                              # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=180G                         # mem in gb
    #SBATCH -t 7-0:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J essai_4                 		   # Name of job
    
    module load gcc/7.3.0 openmpi/3.1.4 samtools jellyfish salmon trinity/2.9.0
    Trinity --seqType fq --max_memory 160G --CPU 10 --left 003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz --right 003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz

### SPAdes

Code - Try #3

    #!/bin/bash
    #SBATCH -c 6                               # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=180G                         # mem in gb
    #SBATCH -t 5-0:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J essai_3_spades                  # Name of job
    
    module load gcc/7.3.0 spades/3.13.1
    spades.py --rna -1 003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz -2 003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz -o /home/alemi055/scratch/ete2020/RNA_Arctic/paired_trimmed/spades
    
## 5. Get rid of all the non-viral contigs

Code - Try #1

    #!/bin/bash
    #SBATCH -c 10                           # Number of CPUS requested. If omitted, the default is 1 CPU.
    #SBATCH --mem=50G                       # mem in gb
    #SBATCH -t 5-0:0:0                      # How long will your job run for? If omitted, the default is 3 hours.
    #SBATCH -J db                           # Name of job
    
    #Build the database
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*tar.gz
    #After that, extract all the files with tar xvzf
