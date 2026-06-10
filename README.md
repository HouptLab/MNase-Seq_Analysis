# MNase-Seq Analysis

## download sequence files from RCC

```
ssh thoupt@hpc-login.rcc.fsu.edu
cd /gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/
rsync -avP *.fastq.gz USERNAME@pauper.bio.fsu.edu:~/FOLDERNAMEOFCHOICE
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

### Pauper

```
mkdir fastqc_raw 
for i in *fastq*; do fastqc $i -t 15 -o fastqc_raw/; done &> fastqc_raw.log 
```

### Macos

install ```brew install fastqc```
install ```brew install parallel```

Script runs against all fastq.gz files in source directory, uses ```parallel``` for speed up, logs fastqc messages to fastqc_raw.log

```
./do_fastqc.sh <source_directory> <fastqc_output_directory>

```

On MacStudio for 2 samples with R1 and R2 (so 4 fastq files) about 1 hour

## Trim reads with Trimmomatic

https://pmc.ncbi.nlm.nih.gov/articles/PMC4103590/
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
http://www.usadellab.org/cms/index.php?page=trimmomatic

### Pauper

run in same directory as fastq.gz files
```
cd ~/medulla_analysis2/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla
for i in *_R1*; do java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 $i ${i/R1/R2} ${i/R1/R1_paired} ${i/R1/R1_unpaired} ${i/R1/R2_paired} ${i/R1/R2_unpaired} ILLUMINACLIP:/home/sn23h/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:1:TRUE MINLEN:25 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15; done &> trimming_B.log & 
```

You can monitor progress with ```tail -f trimming_B.log```.

### Macos

TODO: put command into shell script, with timestamps

Download jar from [Trimmomatic releases](https://github.com/usadellab/Trimmomatic/releases): version 0.40 has parallel unzipping.

adapter files are in ```Trimmomatic-0.40/adapters```, and Trimmomatic looks there automatically.

Query: which are appropriate adapters? NEBNext_PE from the library kit?


NB: add trailing ```&``` o run in background.

```

./do_trimming.sh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla &

```

You can monitor progress with ```tail -f trimming_B.log```.

## Align with bowtie2

[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

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
