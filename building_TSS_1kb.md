# Building TSS_1kb.bed file

MNase-seq of nucleosome sensitivity in the regions flanking transcription start sites (TSS), or TSS-seq, requires filtering aligned reads to those flanking regions, which requires a bed file with the location of flanking regions for every gene in the reference assembly.

We are using rat tissue, using assembly GRCr8, so need a GRCr8 specific listing. These are the steps for producing  `GRCr8_TSS_1kb.bed`:

1. get gene transcription start sites for GRCr8 from NCBI as GTF file
2. make a file containing the chromosome sizes from GRCr8.fa
3. use bedtools to make a TSS bed file with ±1000bp spans from the TSS file and the chromosome sizes file

## Dependencies

1. install NCBI [`datasets` and `dataformat`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).

```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x datasets dataformat
mv datasets usr/local/bin
mv dataformat usr/local/bin
```

## Get list of gene_biotypes and transcript_biotypes

```bash
grep -v "^#" genomic.gtf | awk -F\t '$3=="gene"' | grep -o 'gene_biotype "[^"]*"' | sort | uniq -c
```
yields the list in column 9; only select the "protein_coding" genes
```
 1 gene_biotype "antisense_RNA"
  10 gene_biotype "C_region"
10879 gene_biotype "lncRNA"
 466 gene_biotype "miRNA"
  13 gene_biotype "misc_RNA"
  34 gene_biotype "ncRNA"
23154 gene_biotype "protein_coding"
8687 gene_biotype "pseudogene"
   1 gene_biotype "RNase_MRP_RNA"
 159 gene_biotype "rRNA"
1564 gene_biotype "snoRNA"
1013 gene_biotype "snRNA"
   1 gene_biotype "SRP_RNA"
   1 gene_biotype "telomerase_RNA"
 117 gene_biotype "transcribed_pseudogene"
 771 gene_biotype "tRNA"
 491 gene_biotype "V_segment"
```
 
```bash
grep -v "^#" genomic.gtf | awk -F\t '$3=="transcript"' | grep -o 'transcript_biotype "[^"]*"' | sort | uniq -c
```

yields the list in column 9; not how we want to define "all transcripts"; maybe "mRNA"?
```
1 transcript_biotype "antisense_RNA"
  10 transcript_biotype "C_gene_segment"
21024 transcript_biotype "lnc_RNA"
 796 transcript_biotype "miRNA"
85576 transcript_biotype "mRNA"
  34 transcript_biotype "ncRNA"
 466 transcript_biotype "primary_transcript"
   1 transcript_biotype "RNase_MRP_RNA"
 159 transcript_biotype "rRNA"
1564 transcript_biotype "snoRNA"
1013 transcript_biotype "snRNA"
   1 transcript_biotype "SRP_RNA"
   1 transcript_biotype "telomerase_RNA"
4885 transcript_biotype "transcript"
 771 transcript_biotype "tRNA"
 491 transcript_biotype "V_gene_segment"
 ```
 
2. install bedtools with homebrew: `brew install bedtools`

3.  download the genome sequence records for Rattus norvegicus RefSeq assembly GCF_036323735.1 (GRCr8) as annotated by the NCBI Eukaryotic Genome Annotation Pipeline; this annotation should be referred to as "GCF_036323735.1-RS_2024_02".

```bash
mkdir ncbi
cd ncbi
./datasets download genome accession GCF_036323735.1 --include gtf,genome
unzip ncbi_dataset.zip
# GTF lands at: ncbi/ncbi_dataset/data/GCF_036323735.1/genomic.gtf
```

## GRCr8.chrom.sizes file

The chrom.sizes contig names must match column 1 of the TSS BED (and the BAM alignment files). The TSS BED inherits the GTF's seqnames, so for an NCBI GRCr8 GTF those are RefSeq accessions. The aliases for the chromosome names are tabled at [UCSC GCF_036323735.1.chromAlias.txt](https://hgdownload.soe.ucsc.edu/hubs/GCF/036/323/735/GCF_036323735.1/GCF_036323735.1.chromAlias.txt); a download is included in this repository.

We can extract the chromsome sizes from `GRCr8.fa`. Make a symlink named `GRCr8.fa` to `./ncbi/ncbi_dataset/data/GCF_036323735.1GCF_036323735.1_GRCr8_genomic.fna`, just to save on typing. 

* run `samtools faidx GRCr8.fa` to build an index of GRCr8.fa, puts it in GRCr8.fa.fai
* run `cut -f1,2 GRCr8.fa.fai > GRCr8.chrom.sizes`, which takes the first 2 columns of the index (name of region and length), which will make a valid BED format file of chromosome sizes.

```bash
samtools faidx GRCr8.fa 
cut -f1,2 GRCr8.fa.fai > GRCr8.chrom.sizes
```

## TSS BED files 

Run `make_tss_bed.zsh` to derive strand-aware TSS BED from `genomic.gtf`,

```bash
./make_tss_bed.zsh [--no-drop] <input.gtf[.gz]> <output_TSS.bed> [chrom.sizes] [gene|transcript] 
# generates "<output>_TSS.bed" and "<output>_TSS_1kb.bed"
```

TODO: script undercounts number of mRNA transcripts? some bug in awk program? by "gene" "protein_coding" is ok.

`make_tss_bed.zsh` takes 4 arguments and one flag:

1. the path to GRCr8 gtf file.
2. the output BED file (e.g. `GRCr8_TSS.bed`), which will contain just the first base of each gene/transcript from the gtf file.
3. the path to the chromosome sizes file (e.g. `GRCr8.chrom.sizes`), required to make valid +/- 1kb flanking regions. If the chromosome sizes file is provided, then the +/-1kb TSS bed file (e.g. `GRCr8_TSS_1kb.bed`) is generated. If no chromosome sizes file is provided, then only the TSS bed file is generated. 
4. An optional 4th argument as "transcript" for one TSS per transcript (captures alternative TSSs; multiple per gene), or as "gene" if you want a single TSS per gene instead. The default if unspecified is "gene".
5. a `--no-drop` flag that will retain TSSs that are within 1kb of the chromosome ends (which appear as truncated flanking regions in the +/-1kb TSS bed file.) The default is that TSSs too close to the ends are dropped from the final BED files.

The `make_tss_bed.zsh` script uses `awk` to extract the TSS for each gene from `genomic.gtf`. If the chromosome sizes file is provided, the script then calls `bedtools slop` to produce a bed file with 1kb flanking spans. (Currently, the size of the flank is hardcoded in the script as `FLANK=1000`.)

```bash
./make_tss_bed.zsh ncbi/ncbi_dataset/data/GCF_036323735.1/genomic.gtf GRCr8_TSS.bed   GRCr8.chrom.sizes
```

THe script produces 2 files: `GRCr8_TSS.bed` with just TSS of each gene (i.e. 1 base wide), and `GRCr8_TSS_1kb.bed`, with span from 1kb upstream to 1kb downstream of TSS. Use `GRCr8_TSS_1kb.bed` with `samtools` to filter aligned BAM reads to just those with the +/- 1kb TSS span.

## Notes on numbering

* The GTF source file is 1 indexed, while BED files are zero-indexed. So, the `awk` routine in the `make_tss_bed.zsh` script moves the base-numbering down by one.

* BED files specify the end of a region with the index of the first base *after* the last base in the region; i.e., the start of the range is *inclusive*, the end of the range is *non-inclusive*, or  `[chromStart, chromEnd)`. 

  So for a gene like `Npy`, the BED line for the entire transcript is:

```bed
NC_086022.1	80212110	80219309	Npy	.	+
```

  Transcription start (the first base of the transcript) is at base 80212110 of chromosome 4 (the chromStart field). The last base of the transcript is base 80212110 (1 base before the BED chromEnd field). `Npy` is on the + strand, so it runs from lower base index to higher base index.

  For a gene like `Th`, the BED line for the entire transcript is:

```bed
NC_086019.1	207500958	207509276	Th	.	-
```

  Transcription start (the first base of the transcript, "+1 == TSS") is at 207509275 (because it is on the - strand, so it runs from higher index to lower index, but because range is non-incluse it starts at chromEnd-1).

* The +/- 1kb TSS region contains 2001 bases: 1000 bases upstream from the transcription start base, the  first transcribed base itself ("+1 == TSS"), and 2000 bases upstream from the transcription start base. So the BED range looks like chromStart = TSS - 1000, chromEnd = TSS + 1000. In eukaryotic transcription start site numbering convention, the region spans −1000 TSS to +1001 TSS, inclusive.

  So for `Npy`, with transcription start base = 80212110, the `GRCr8_TSS_1kb.bed` entry will be.

```bed
NC_086022.1	80211110	80213111	Npy	.	+
```

* By default, if a TSS is too close to the ends of a chromosome (e.g., with 1 kb of either 0 or chromosome size), then the TSS is dropped from the TSS BED file and the +/- 1kb BED file.  If the `--no-drop` flag is set, then all TSS are retained, and  `bedtools slop`  generates an entry for the TSS but truncates the flanking region. So a gene with TSS at 500, with +/- 1kb from -500 to 1501, will be truncated to 0 to 1501. 
