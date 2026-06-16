# MNase-Seq Analysis

## Programs and Dependencies

- Trimmomatic 0.40 
  
  download java jar from [Trimmomatic releases](https://github.com/usadellab/Trimmomatic/releases), place in /Applications
  
- fastqc  0.12.1 

 `brew install fastqc`
  
- bowtie2 2.5.5

  `brew install bowtie2`
  
- parallel GNU parallel 20260522 

  `brew install parallel`
  
- samtools 1.23 

  `brew install samtools`
  

## Reference sequence files

- bowtie2 indexes for GRCr8

  download indexes from https://benlangmead.github.io/aws-indexes/bowtie (bowtie2 maintainer site), put into `sequences/bowtie2_indexes/GRCr8/`

```bash
curl -L -o sequences/bowtie2_indexes/GRCr8.zip \
           https://genome-idx.s3.amazonaws.com/bt/GRCr8.zip
```

- GRCr8_TSS_1kb.bed

  see `build_TSS_1kb.md` in this repository. The 1kb flanking TSS regions file is made by `make_tss_bed.zsh` from GRCr8 annotation files.

## download sequence files from RCC

```bash
ssh thoupt@hpc-login.rcc.fsu.edu
cd /gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/
rsync -avP *.fastq.gz USERNAME@pauper.bio.fsu.edu:~/FOLDERNAMEOFCHOICE
```

or, copy to local directory
```bash
rsync -avP thoupt@hpc-login.rcc.fsu.edu:/gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/Thomas_Houpt_11-19-2025_SN_Medull ./
```

view first 10 lines

```bash
gzcat SN_Medulla_10U_S1_L008_R1_001.fastq.gz | head -n 10
```

## Check quality with fastqc 

### Pauper

```bash
mkdir fastqc_raw 
for i in *fastq*; do fastqc $i -t 15 -o fastqc_raw/; done &> fastqc_raw.log 
```

To download and view the fastqc.html files, use ```rsync```

```bash
rsync -avc sn23h@pauper.bio.fsu.edu:~/medulla_analysis2/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/fastqc_raw/*.html .
```


### Macos

* install ```brew install fastqc```
* install ```brew install parallel```

Script runs against all fastq.gz files in source directory, uses ```parallel``` for speed up, logs fastqc messages to fastqc_raw.log

```bash
./do_fastqc.sh <source_directory> <fastqc_output_directory>
```

On MacStudio for 2 samples with R1 and R2 (so 4 fastq files) about 1 hour

## Trim reads with Trimmomatic

https://pmc.ncbi.nlm.nih.gov/articles/PMC4103590/
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
http://www.usadellab.org/cms/index.php?page=trimmomatic

### Pauper

run in same directory as fastq.gz files

```bash
cd ~/medulla_analysis2/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla
nohup bash -c 'for i in *_R1*; do java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 -phred33 "$i" "${i/R1/R2}" "${i/R1/R1_paired}" "${i/R1/R1_unpaired}" "${i/R1/R2_paired}" "${i/R1/R2_unpaired}" ILLUMINACLIP:/home/sn23h/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:1:TRUE MINLEN:25 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 < /dev/null; done' > trimming_B.log 2>&1 &
```

You can monitor progress with ```tail -f trimming_B.log```.

### Macos


Download jar from [Trimmomatic releases](https://github.com/usadellab/Trimmomatic/releases): version 0.40 has parallel unzipping.

adapter files are in ```Trimmomatic-0.40/adapters```, and Trimmomatic looks there automatically.

Query: which are appropriate adapters? NEBNext_PE from the library kit?


```bash
./trim_pe.zsh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla

```


You can monitor progress with ```tail -f trimming.log```.

**put the original sequencing files in ```/raw```, and the paired/unpaired files in ```/trimmed```

## Check quality of trimmed with fastqc 

### Macos

```bash
./do_fastqc.sh ./do_fastqc.sh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/trimmed ./fastqc_trimmed

```

## Align with bowtie2

[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

We align to current (2024) reference genome assembly GRCr8 for rat 
* [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610589/) 
* [NCBI site](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036323735.1/)
* [assembly report (txt download)](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_assembly_report.txt
) -- gives ascension numbers for each chromosome and mitochondrion

### Macos

install `brew install bowtie2`

- need to get bowtie2 indices for rat genome (GRCr8 for rat)

download indexes from https://benlangmead.github.io/aws-indexes/bowtie (bowtie2 maintainer site), put into sequences/bowtie2_indexes/GRCr8/

to copy bowtie2 indexes to pauper, use curl:

```bash
curl -L -o sequences/bowtie2_indexes/GRCr8.zip https://genome-idx.s3.amazonaws.com/bt/GRCr8.zip
```



run bowtie2 and pipe through samtools to get BAM files:

```bash
nohup ./run_bowtie2.zsh <source_directory> <destination_directory> &
```

e.g.

```bash
nohup ./run_bowtie2.zsh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/trimmed ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/aligned &
 ```

Outputs BAM files to the destination directory. Logs to bowtie.log (and  bowtie2 logs per-sample bowtie2.log). 

* specify non-discordant and no-mixed
* -x: Specifies the index (use the prefix).
* -1, -2: Your forward and reverse read files (can be gzipped).
* -S: Output SAM file. (not used, because outpur piped through samtools to make BAM file)
* -p 8: Uses 8 threads for faster alignment (adjust as needed). 



