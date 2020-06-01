# FastQC

### Download the fastqc_script.sh to a directory.
I put it in /scratch. 

```scp /drives/c/Users/Audrée/Downloads/fastqc_script.sh alemi055@cedar.computecanada.ca:~/scratch/```

### Move to the "scratch" directory
```cd /home/alemi055/scratch/```

### Convert the Windows file to a Unix file
```dos2unix fastqc_script.sh```

### Submit the script as a sbatch job
```sbatch fastqc_script.sh```
