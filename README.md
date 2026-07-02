# MNase-Seq Analysis

Based on notes from Dr. Jane Benoit (alum of Dennis lab) and scripted with assistance from Claude Opus 4.8.

Documentation of workflow for TSS analysis of MNase-seq data, with the goal of mapping changes in nucleosome occupancy and sensitivity in the rat genome, e.g. in the adrenal medulla after stress induction.

This repository organizes and codifies the set-up and processing steps, with a series of shell scripts containing the command line tools and logging results. Ultimately we will have a makefile to specify the build and automate the data pipeline from raw sequencing files to heatmaps and gene ontology mapping. To ensure reproducibility, tool versions (and RNG seeds) are recorded.

---

* TODO: add citations for all programs used
* TODO: add links to manuals/websites for all programs used
* TODO: add a flowchart of dataflow with data files, transforms, and tools used to transforms
* TODO: make sure each step has a functional description
* TODO: add time estimate for each step
* TODO: add approximate file sizes for each stem
* TODO: describe how to get estimated read count from sequencing quality files

---

Contents of README:

- [Overview of workflow](#overview-of-workflow)
- [Set up](#set-up)
- [Download sequence files from RCC](#download-sequence-files-from-RCC)
- [Check quality with fastqc](#check-quality-with-fastqc)
- [Trim reads with Trimmomatic](#trim-reads-with-trimmomatic)
- [Check quality of trimmed reads](#check-quality-of-trimmed-reads)
- [Align reads with bowtie2](#align-reads-with-bowtie2)
- [Remove duplicate reads](#remove-duplicate-reads)
- [Downsample BAM files](#downsample-bam-files)
- [Run QC with deepTools](#run-qc-with-deeptools)
- [Call nucleosome positions](#call-nucleosome-positions)

---

## Overview of workflow

1. download raw sequencing files -> `/raw/*.fasta.gz`
2. check quality with `fastqc` -> `/fastqc_raw/`
3. trim sequencing adapters and barcodes with `Trimmomatic` -> `/trimmed/\*.fasta.gz`
4. check quality with `fastqc`  -> `/fastqc_trimmed/`
5. align reads with `bowtie2` -> `/aligned/\*.bam`
6. remove duplicate reads with `samtools` -> `/nodups/\*.bam`
7. downsample with `Picard` -> `/downsampled/\*.bam`
8. merge heavy & light files for total occupancy -> `/downsampled/\*.bam`
9. check quality with `deepTools`
10. call nucleosome positions with `DANPOS3`, filtered by TSS -> `/nucleosomes/\*.bw`
11. get log2ratio of light/heavy with `bigwigCompare` -> `/nucleosomes/\*.bw`
12. get deepTools matrix file with `computeMatrix` -> `/nucleosomes/\*.matrix.gz`
13. get clusters, heatmaps and profiles with `plotHeatmap`/`plotProfiles` -> `/results/*.png`
14. get GO plot with `clusterProfiler` -> `/results/`


## Set up

### Shell Scripts

From this repository:

- do_fastqc.sh
- trim_pe.zsh
- run_bowtie2.zsh
- markdup_dedup.zsh
- downsample_bams_picard.zsh


### Programs and Dependencies

- Trimmomatic 0.40 
  
    download java jar from [Trimmomatic releases](https://github.com/usadellab/Trimmomatic/releases), place in `/Applications`; version 0.40 has parallel unzipping.
  
- fastqc  0.12.1 

    `brew install fastqc`
  
- bowtie2 2.5.5

    `brew install bowtie2`
  
- parallel GNU parallel 20260522 

    `brew install parallel`
  
- samtools 1.23 

    `brew install samtools`
    
- Picard 3.4.0

    download java jar from [Picard releases](https://github.com/broadinstitute/picard/releases/tag/3.4.0), place in `/Applications`
    
- deepTools

   install [deepTools with Anaconda](https://deeptools.readthedocs.io/en/latest/content/installation.html)
   
- DANPOS3

  install as [a bioconda package](https://github.com/boenc28-cmyk/DANPOS)
  
- clusterProfiler

  https://guangchuangyu.github.io/software/clusterProfiler/
---

### Reference sequence files

To align reads with the rat reference genome, and to limit analysis to transcription start site (TSS) of protein-encoding genes, we need to provide some reference files.

We align to current (2024) reference genome assembly GRCr8 for rat 
* [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610589/) 
* [NCBI site](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036323735.1/)
* [assembly report (txt download)](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_assembly_report.txt) -- gives ascension numbers for each chromosome and mitochondrion


- `bowtie2` indexes for GRCr8

    download indexes from [bowtie2 maintainer site]( https://benlangmead.github.io/aws-indexes/bowtie), put into `sequences/bowtie2_indexes/GRCr8/`

    ```bash
    curl -L -o sequences/bowtie2_indexes/GRCr8.zip \
               https://genome-idx.s3.amazonaws.com/bt/GRCr8.zip
    ```

- GRCr8_TSS_pc_1kb.bed

  see `build_TSS_1kb.md` in this repository. The file of 1kb flanking TSS regions of protein-coding genes is made by `make_tss_bed.zsh` from GRCr8 annotation files.

### Directory organization

```
project
    ├── aligned # BAM files aligned by bowtie2
    ├── downsampled # BAM files normalized by Picard
    ├── fastqc_raw
    ├── fastqc_trimmed
    ├── nodup # duplicates removed by samtools
    ├── nucleosomes # peaks called by DANPOS3 and deepTools matrix
    ├── raw # sequencing files
    ├── results  # heatmaps, profiles, clusters, GO output  
    └── trimmed # fasta.gz files from Trimmomatic
```


---

## Download sequence files from RCC

The sequencing core will tell you where your data is located, e.g.

```bash
ssh USERNAME@hpc-login.rcc.fsu.edu
cd /gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/
rsync -avP *.fastq.gz USERNAME@pauper.bio.fsu.edu:~/FOLDERNAMEOFCHOICE
```

or, copy to local directory

```bash
rsync -avP thoupt@hpc-login.rcc.fsu.edu:/gpfs/research/medicine/sequencer/NovaSeqXPlus/Outputs_XP/2025_Outputs_XP/Thomas_Houpt_11-19-2025_SN_Medull ./
```

view first 10 lines of a FASTQ sequencing file

```bash
gzcat SN_Medulla_10U_S1_L008_R1_001.fastq.gz | head -n 10
```

*put the original sequencing files in `/raw`*

---

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

Script runs against all fastq.gz files in source directory, uses ```parallel``` for speed up, logs fastqc messages to fastqc_raw.log

```bash
./do_fastqc.sh <source_directory> <fastqc_output_directory>
```

On MacStudio for 2 samples with R1 and R2 (so 4 fastq files) about 1 hour

*put into /fastqc_raw directory*

---

## Trim reads with Trimmomatic

*TODO: a little discussion of what gets trimmed (adapters and barcodes).*


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

Query: which are appropriate adapters? NEBNext_PE from the library kit?

`./trim_pe.zsh` script runs Trimmomatic with `-phred33` and
 `ILLUMINACLIP:${ADAPTERS}:2:30:10:1:TRUE MINLEN:25 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15`. The adapters are hardcoded in the script as` TruSeq3-PE.fa`. The adapter files are in ```Trimmomatic-0.40/adapters```, and Trimmomatic looks there automatically.

```bash
./trim_pe.zsh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla
```


You can monitor progress with `tail -f trimming.log`.

*put the trimmed paired/unpaired files in `/trimmed`*


---

## Check quality of trimmed reads

### Macos

```bash
./do_fastqc.sh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/trimmed ./fastqc_trimmed
```

*put results in `./fastqc_trimmed` directory*

---

## Align reads with bowtie2

### Macos

to copy `bowtie2` indexes to pauper, use curl:

```bash
curl -L -o sequences/bowtie2_indexes/GRCr8.zip https://genome-idx.s3.amazonaws.com/bt/GRCr8.zip
```

and place in `./bowtie2_indexes/GRCr8`


The script `run_bowtie2.zsh` runs `bowtie2` and pipes SAM output through `samtools` to get BAM files:

```bash
nohup ./run_bowtie2.zsh <source_directory> <destination_directory> &
```
Note that because the alignment can take dozens of hours, we use `nohup` and `&` to run in the background (`&`)  and continue running even if we hangup by closing terminal (`nohup`).

e.g.

```bash
nohup ./run_bowtie2.zsh ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/trimmed ./sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/aligned &
```

Outputs BAM files to the destination directory. Logs to `bowtie.log` (and  `bowtie2` itself logs into per-sample `bowtie2.log`s). 

The `bowtie2` invocation specifies:
* non-discordant and no-mixed
* -x to specify the index (use the prefix). Index location currently hardcoded to "$SCRIPT_DIR/bowtie2_indexes/GRCr8"
* -1, -2: Your forward and reverse read files (can be gzipped).
* -p : Uses number of cores for threads for faster alignment, or adjust as needed with THREADS env variable. 

The `bowtie2` alignment results are piped to `samtools` to directly produce sorted BAM files, with reads with quality less than 10 dropped (`-q 10`). For each generated BAM file, `samtools index` is  called to generate bam.bai index files, and `samtools flagstat` is called to provide summary statistics. Downstream tools will use the bam.bai index files to speed up random-access into the BAM files during processing.

To view BAM file contents:

```bash
samtools view input.bam | head -10        # first 10 alignment records
samtools view -h input.bam | head -10     # include header lines (@HD, @SQ, etc.)
samtools head input.bam                    # header only
```

To copy to pauper:

```
rsync -avP -c  houpt@bio-k2067c-mac.bio.fsu.edu:/Users/houpt/Programming_Github/MNase-Seq_Analysis/sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/aligned/ ./aligned 2>&1 | grep -i -E 'error|denied|failed|permission'
```

---

## Remove duplicate reads

*TODO: a little discussion of what duplicate reads are and how they are identified.*

Use [`samtools markdup`] (https://www.htslib.org/doc/samtools-markdup.html) to remove identical  duplicate reads (PCR artifacts?)

```bash
./markdup_dedup.zsh /path/to/bams /path/to/dedup
```

e.g. 

```bash
./markdup_dedup.zsh /sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/aligned /sequences/Thomas_Houpt_05-29-2026_Houpt_SN_Medulla/Houpt_SN_Medulla/dedup
```


Defaults: INPUT_DIR=current dir, OUTPUT_DIR=./dedup. A log file is written at `markdup_$(date +%Y%m%d_%H%M%S).log`

Override thread count with the THREADS environment variable (or modify script to call `hw.perflevel0.physicalcpu` to get count of performance cores only, if you want to avoid loading efficiency cores ).

```
THREADS=8 ./markdup_dedup.zsh /path/to/bams /path/to/output
```

This script runs over all BAM files in the source directory and applies:

```bash
samtools fixmate -m Sorted_names.bam Fixmate.bam
samtools sort -o Sorted.bam Fixmate.bam
samtools markdup -r -s Sorted.bam Final_File.bam
```

`samtools fixmate -m` requires name-collated (name-sorted) input, which is why the script runs `samtools sort -n` first. The manual states `fixmate` should be run on a name-sorted/name-collated file and that the `-m `option adds the mate score tags needed by `markdup`. `samtools markdup` requires position-sorted input with the `fixmate` tags present, which the script produces with the second `samtools sort` before `markdup`. The `-r` flag removes duplicates and `-s` prints statistics.

---

## Downsample BAM files

*TODO: a little discussion of why we downsample.*


Normalize number of paired reads in all BAM files to the number of paired reads in the smallest BAM file, using Picard:

```bash
./downsample_bams_picard.zsh /path/to/bams /path/to/downsampled

# or reset some script values
PICARD_JAR=/path/picard.jar SEED=42 STRATEGY=HighAccuracy \
      ./downsample_bams_picard.zsh /path/to/bams /path/to/downsampled
```

Defaults:
- INPUT_DIR  = current dir
- OUTPUT_DIR = ./downsampled
- SEED       = 42
- STRATEGY   = HighAccuracy   (better adherence to target proportion)
- ACCURACY   = 0.0001
- PICARD_JAR = /Applications/picard.jar

The script counts the number of paired reads in each BAM file using `samtools`:

```bash
samtools view -c -f 0x40 -F 0x90C -@ "$THREADS" "$bam"
```

This uses `-f 0x40 -F 0x90C` to count templates (first-in-pair, primary, both mates mapped) rather than reads, so the denominator matches what Picard samples. Picard's `PROBABILITY` is a per-template keep probability — its docs state the goal is retaining reads from `PROBABILITY` × (input templates), so P = `target_min / this_file's_template_count`.

For the Picard invocation, `STRATEGY=HighAccuracy` is the default because the Picard docs recommend it for smaller inputs: `ConstantMemory` should be accurate 99.9% of the time when the input contains ≥ 50,000 templates; for smaller inputs `HighAccuracy` is recommended instead. Override with `STRATEGY=ConstantMemory` if your libraries are large and memory is a concern. `RANDOM_SEED` is set for reproducibility. 

A final pass with `samtools index` creates `downsampled.bam.bai` index files.

*Note: can't confirm Picard's memory needs for particular BAM sizes under `HighAccuracy`; if you hit Java heap errors, add `-Xmx` (e.g. `java -Xmx8g -jar` ...) so we can debug.*


---


## Run QC with deepTools


---


## Call Nucleosome Positions

*TODO: include some suggestions on how to characterize the position data, e.g. volcano plots, histograms, etc.*


---


## Clustering, Heatmaps, and Profiles

run `plotHeatmap --kmeans <N> ---outFileSortedRegions my_clustered_regions.bed` to get list of clustered regions. Need to run `bedtools intersect` to get gene names? What does `clusterProfiler` want for input?

### Notes on reproducibility of clustering:

Because plotHeatmap uses scikit-learn's k-means clustering under the hood, the cluster initialization centroids are chosen randomly each time (there is no way to specify a seed for the randomization). If you are running `plotHeatmap --kmeans <N>`, your rows will likely be grouped and ordered slightly differently every single time you execute the command.

To fix clustering, either use hierarchical clustering `--hclust <N>`, which is deterministic, or generate clusters once and reuse it by giving computeMatrix the regions BED file, then re-plot without clustering:

```bash
plotHeatmap -m matrix.gz -out temporary_heatmap.png --kmeans 4 --outFileSortedRegions my_clustered_regions.bed
computeMatrix reference-point -S signal.bigWig -R my_clustered_regions.bed -o reproducible_matrix.gz
plotHeatmap -m reproducible_matrix.gz -out reproducible_heatmap.png
```

---



## Gene Ontology mapping
