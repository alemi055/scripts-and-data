Code - Try #1

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

Code - Try #3

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    GR_RNA_BS3-3_R1_paired.fastq.gz GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    GR_RNA_BS3-3_R2_paired.fastq.gz GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36 CROP:120 HEADCROP:20 AVGQUAL:20
    
Code - Try #4

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    GR_RNA_BS3-3_R1_paired.fastq.gz GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    GR_RNA_BS3-3_R2_paired.fastq.gz GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:100 HEADCROP:20 AVGQUAL:20
    
Code - Try #5

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    Essai5_R1_paired_GR_RNA_BS3-3.fastq.gz Essai5_R1_unpaired_GR_RNA_BS3-3.fastq.gz \
    Essai5_R2_paired_GR_RNA_BS3-3.fastq.gz Essai5_R2_unpaired_GR_RNA_BS3-3.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:120 HEADCROP:8 AVGQUAL:20
    
Code - Try #6

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    Essai6_R1_paired_GR_RNA_BS3-3.fastq.gz Essai6_R1_unpaired_GR_RNA_BS3-3.fastq.gz \
    Essai6_R2_paired_GR_RNA_BS3-3.fastq.gz Essai6_R2_unpaired_GR_RNA_BS3-3.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:100 HEADCROP:8 AVGQUAL:20
    
Code - Try #7

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    Essai7_R1_paired_GR_RNA_BS3-3.fastq.gz Essai7_R1_unpaired_GR_RNA_BS3-3.fastq.gz \
    Essai7_R2_paired_GR_RNA_BS3-3.fastq.gz Essai7_R2_unpaired_GR_RNA_BS3-3.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20
    
Code - Try #8

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    Essai8_R1_paired_GR_RNA_BS3-3.fastq.gz Essai8_R1_unpaired_GR_RNA_BS3-3.fastq.gz \
    Essai8_R2_paired_GR_RNA_BS3-3.fastq.gz Essai8_R2_unpaired_GR_RNA_BS3-3.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:110 HEADCROP:20 AVGQUAL:20
    
Code - Try #9

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
    HI.4444.003.Index_3.GR_RNA_BS3-3_R1.fastq.gz HI.4444.003.Index_3.GR_RNA_BS3-3_R2.fastq.gz \
    003.Index_3.GR_RNA_BS3-3_R1_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R1_unpaired.fastq.gz \
    003.Index_3.GR_RNA_BS3-3_R2_paired.fastq.gz 003.Index_3.GR_RNA_BS3-3_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:3:26:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:105 HEADCROP:15 AVGQUAL:20
