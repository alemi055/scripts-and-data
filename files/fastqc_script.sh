#!/bin/bash
#SBATCH -c 5                               # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=50gb                         # mem in gb
#SBATCH -t 0-4:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
#SBATCH -J fastqc                 	# Name of job

cd /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data

module load fastqc

fastqc -t 5 *.fastq.gz

cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*html /home/alemi055/scratch/ete2020/RNA_Arctic/quality_control/

cp /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data/*zip /home/alemi055/scratch/ete2020/RNA_Arctic/quality_control/

cd /home/alemi055/scratch/ete2020/RNA_Arctic/raw_data

rm *html

rm *zip
