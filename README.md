# 1. Download raw data
    read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://genomequebec.mcgill.ca/nanuqMPS/readsetList?projectId=15225&tech=HiSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

We will need to enter our Nanuq username and password.


# 2. FastQC

### Download the fastqc_script.sh to a directory.
I put it in /scratch. 

```scp /drives/c/Users/Audrée/Downloads/fastqc_script.sh alemi055@cedar.computecanada.ca:~/scratch/```

### Move to the "scratch" directory
```cd /home/alemi055/scratch/```

### Convert the Windows file to a Unix file
```dos2unix fastqc_script.sh```

### Submit the script as a batch job
```sbatch fastqc_script.sh```


# 3. Trimmomatic
   
Code - Try #9

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    Essai9_R1_paired_GR_RNA_BS3-3.fastq.gz Essai9_R1_unpaired_GR_RNA_BS3-3.fastq.gz \
    Essai9_R2_paired_GR_RNA_BS3-3.fastq.gz Essai9_R2_unpaired_GR_RNA_BS3-3.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20

GENERIC CODE

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.00!.Index_!.GR_RNA_!_R1_fastq.gz HI.4444.00!.Index_!.GR_RNA_!_R2_fastq.gz \
    00!.Index_!.GR_RNA_!_R1_paired_fastq.gz 00!.Index_!.GR_RNA_!_R1_unpaired_fastq.gz \
    00!.Index_!.GR_RNA_!_R2_paired_fastq.gz 00!.Index_!.GR_RNA_!_R2_unpaired_fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20
    
*Replace "!" by appropriate number/letter*
