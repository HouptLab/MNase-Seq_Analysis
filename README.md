# MNase-Seq Analysis

## download sequence files from RCC

```

> ssh thoupt@hpc-login.rcc.fsu.edu
> cd /gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/
> rsync -avP *.fastq.gz USERNAME@pauper.bio.fsu.edu:~/FOLDERNAMEOFCHOICE
```

or, copy to local directory
```
rsync -avP thoupt@hpc-login.rcc.fsu.edu:/gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/Thomas_Houpt_11-19-2025_SN_Medull ./
```

# view first 10 lines

```
gzcat SN_Medulla_10U_S1_L008_R1_001.fastq.gz | head -n 10
```

## Check quality with fastqc 

install ```brew install fastqc```
install ```brew install parallel```

Script runs against all fastq.gz files in source directory, uses ```parallel``` for speed up, logs fastqc messages to fastqc_raw.log

```
./do_fastqc.sh <source_directory> <fastqc_output_directory>

```

## Trim reads with Trimmomatic



## Align with bowtie2

(bowtie2)[https://bowtie-bio.sourceforge.net/bowtie2/index.shtml]

install ```brew install bowtie2```

- need to get bowtie2 indices for rat genome (GRCr8 for rat)

download indexes from https://benlangmead.github.io/aws-indexes/bowtie (bowtie2 maintainer site), put into sequences/bowtie2_indexes/GRCr8/

to copy bowtie2 indexes to pauper, use curl:

```
curl -L -o /path/to/destination/GRCr8.zip https://genome-idx.s3.amazonaws.com/bt/GRCr8.zip
```



run bowtie2:

```
bowtie2 -x bowtie2_indexes/GRCr8 -1 reads_R1.fq.gz -2 reads_R2.fq.gz -S GRCr8.1_alignment.sam -p 8
```


* -x: Specifies the index (use the prefix).
* -1, -2: Your forward and reverse read files (can be gzipped).
* -S: Output SAM file.
* -p 8: Uses 8 threads for faster alignment (adjust as needed). 
