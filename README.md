# MNase-Seq Analysis

## Programs and Dependencies

- Trimmomatic 0.40
- fastqc  0.12.1
- bowtie2 2.5.5
- parallel GNU parallel 20260522
- samtools 1.23
- datasets 18.30.1 (need to place in /bin or similiar)
- dataformat ?
- bedtools  v2.31.1

## Reference sequence files

- bowtie2 indexes for GRCr8
- GRCr8_TSS_1kb.bed

### for making GRCr8_TSS.bed and GRCr8_TSS_1kb.bed
- ncbi/ncbi_dataset/data/GCF_036323735.1/genomic.gtf
- ./ncbi/ncbi_dataset/data/GCF_036323735.1GCF_036323735.1_GRCr8_genomic.fna
- ncbi/ncbi_dataset/data/GCF_036323735.1/GRCr8.chrom.sizes

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

install ```brew install bowtie2```

- need to get bowtie2 indices for rat genome (GRCr8 for rat)

download indexes from https://benlangmead.github.io/aws-indexes/bowtie (bowtie2 maintainer site), put into sequences/bowtie2_indexes/GRCr8/

to copy bowtie2 indexes to pauper, use curl:

```bash
curl -L -o /path/to/destination/GRCr8.zip https://genome-idx.s3.amazonaws.com/bt/GRCr8.zip
```



run bowtie2 and pipe through samtools to get BAM files:

```bash
nohup ./run_bowtie2.zsh <source_directory> <destination_directory> &
```

Outputs BAM files to the destination directory. Logs to bowtie.log (and  bowtie2 logs per-sample bowtie2.log). 

* -x: Specifies the index (use the prefix).
* -1, -2: Your forward and reverse read files (can be gzipped).
* -S: Output SAM file. (not used, because outpur piped through samtools to make BAM file)
* -p 8: Uses 8 threads for faster alignment (adjust as needed). 


## Construct TSS bed file (`GRCr8_TSS_1kb.bed`)

We need to filter reads to only those that are ±1000bp of TSS sites of rat genes in GRCr8.

1. get gene transcription start sites for GRCr8 from NCBI as GTF file
2. use bedtools to make a TSS bed file with ±1000bp spans

* install NCBI [`datasets` and `dataformat`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).

```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x datasets dataformat
mv datasets usr/local/bin
mv dataformat usr/local/bin
```

* install bedtools with homebrew: `brew install bedtools`

* download the genome sequence records for Rattus norvegicus RefSeq assembly GCF_036323735.1 (GRCr8) as annotated by the NCBI Eukaryotic Genome Annotation Pipeline; this annotation should be referred to as "GCF_036323735.1-RS_2024_02".

```bash
mkdir ncbi
cd ncbi
./datasets download genome accession GCF_036323735.1 --include gtf,genome
unzip ncbi_dataset.zip
# GTF lands at: ncbi/ncbi_dataset/data/GCF_036323735.1/genomic.gtf
```

Note: the chrom.sizes contig names must match column 1 of the BED (and your BAM). The TSS BED inherits the GTF's seqnames, so for an NCBI GRCr8 GTF those are RefSeq accessions — generate the matching chrom.sizes from the same NCBI FASTA (e.g. `samtools faidx GRCr8.fa` then `cut -f1,2 GRCr8.fa.fai > GRCr8.chrom.sizes`) so all three agree.

GRCr8.fa is same as ./ncbi/ncbi_dataset/data/GCF_036323735.1GCF_036323735.1_GRCr8_genomic.fna, so make a symlink rather than renaming original, then:

* run `samtools faidx GRCr8.fa` to build an index of GRCr8.fa, puts it in GRCr8.fa.fai
* run `cut -f1,2 GRCr8.fa.fai > GRCr8.chrom.sizes` takes the first 2 columns of the index (name of region and length), which will make a valid bed format file of chromosome sizes.

```bash
samtools faidx GRCr8.fa 
cut -f1,2 GRCr8.fa.fai > GRCr8.chrom.sizes
```

* derive strand-aware TSS BED from genomic.gtf, using `awk`. Specify 4th argument as "transcript" for one TSS per transcript (captures alternative TSSs; multiple per gene), or as "gene" if you want a single TSS per gene instead.

```bash
./make_tss_bed.zsh <input.gtf[.gz]> <output.bed> [chrom.sizes] [gene|transcript] 
# e.g. ./make_tss_bed.zsh ncbi/ncbi_dataset/data/GCF_036323735.1/genomic.gtf GRCr8_TSS.bed   ncbi/ncbi_dataset/data/GCF_036323735.1/GRCr8.chrom.sizes
```
produces 2 files: `GRCr8_TSS.bed` with just TSS of each gene (i.e. 1 base wide), and `GRCr8_TSS_1kb.bed`, with span from 1kb upstream to 1kb downstream of TSS. Use `GRCr8_TSS_1kb.bed` to filter aligned reads to within TSS span.

