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

### Submit the script as a sbatch job
```sbatch fastqc_script.sh```


# 3. Trimmomatic

~~Code - Try #1~~

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    GR_RNA_BS3-3_R1_paired.fastq.gz GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    GR_RNA_BS3-3_R2_paired.fastq.gz GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

Code - Try #2

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    GR_RNA_BS3-3_R1_paired.fastq.gz GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    GR_RNA_BS3-3_R2_paired.fastq.gz GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:20 AVGQUAL:20
