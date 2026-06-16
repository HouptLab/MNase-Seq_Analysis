#!/bin/zsh
#
# make_tss_bed.zsh  (macOS / zsh; also works on Linux)
# Derive a strand-aware, 1-bp TSS BED from a GTF annotation, then optionally
# expand each TSS into a promoter window with bedtools slop.
#
# For each selected feature, the TSS is the 5' end of the feature: the GTF
# start (column 4) on the + strand, or the GTF end (column 5) on the - strand.
# GTF is 1-based; BED is 0-based half-open, so a single-base TSS at position p
# is emitted as [p-1, p). Output is BED6, coordinate-sorted and de-duplicated.
#
# The feature type selects granularity and which ID becomes the BED name:
#   gene        -> one TSS per gene       (name = gene_id)        [default]
#   transcript  -> one TSS per transcript (name = transcript_id)
#
# If a chrom.sizes file is available, the script then runs:
#   bedtools slop -b $FLANK -i <output.bed> -g <chrom.sizes> > <output>_<N>.bed
# producing windows of +/- FLANK bp around each TSS (clipped to chromosome
# ends). FLANK is configured below and defaults to 1000 (-> "..._1kb.bed").
# If no chrom.sizes is found, the slop step is skipped and only the raw TSS
# BED is written.
#
# The seqname (column 1) is passed through unchanged, so the BED's chromosome
# names match those in the GTF (for NCBI GRCr8 GTFs that means RefSeq
# accessions). Those names must match both your chrom.sizes and your BAM.
#
# By default, TSSs within FLANK bp of either chromosome end are dropped before
# slop, so every emitted window is a full +/- FLANK bp with no clamping. Pass
# --no-drop to disable this and keep all TSSs (slop will then clamp edge windows
# to chromosome bounds, as before).
#
# Usage:
#   ./make_tss_bed.zsh [--no-drop] <input.gtf[.gz]> <output.bed> [chrom.sizes] [gene|transcript]

setopt PIPE_FAIL       # a pipeline fails if any stage fails

# ---- Temp-file cleanup -------------------------------------------------------
# Remove scratch files on any exit: normal completion, error, or interruption.
# Each temp file is appended to TMP_FILES as it is created; rm -f is a no-op for
# files that were already mv'd into place on success.
typeset -a TMP_FILES
cleanup() { (( ${#TMP_FILES} )) && rm -f -- "${TMP_FILES[@]}"; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM HUP

# ---- Configuration -----------------------------------------------------------
SCRIPT_DIR=${0:A:h}
FLANK=1000                                    # bp added on each side by slop (-b)
CHROM_SIZES="$SCRIPT_DIR/GRCr8.chrom.sizes"   # default; override with 3rd argument

# ---- Argument handling -------------------------------------------------------
# Optional flags precede the positional arguments. --drop (default) removes
# TSSs too close to a chromosome end before slop; --no-drop keeps them.
DROP_EDGE=1
while [[ $# -gt 0 && $1 == --* ]]; do
    case $1 in
        --no-drop) DROP_EDGE=0; shift ;;
        --drop)    DROP_EDGE=1; shift ;;
        --)        shift; break ;;
        *)
            print -u2 -- "Error: unknown option '$1'."
            print -u2 -- "Usage: $0 [--no-drop] <input.gtf[.gz]> <output.bed> [chrom.sizes] [gene|transcript]"
            exit 1
            ;;
    esac
done

if [[ $# -lt 2 || $# -gt 4 ]]; then
    print -u2 -- "Usage: $0 [--no-drop] <input.gtf[.gz]> <output.bed> [chrom.sizes] [gene|transcript]"
    exit 1
fi

GTF=$1
BED=$2
[[ -n $3 ]] && CHROM_SIZES=$3 # optional chrom.sizes override
FEATURE=${4:-gene}            # default to one TSS per gene

if [[ $FEATURE != "gene" && $FEATURE != "transcript" ]]; then
    print -u2 -- "Error: feature must be 'gene' or 'transcript' (got '$FEATURE')."
    exit 1
fi

if [[ ! -e $GTF ]]; then
    print -u2 -- "Error: input file '$GTF' does not exist."
    exit 1
fi

# Create the output directory if needed.
out_dir=${BED:h}
if [[ -n $out_dir && ! -d $out_dir ]]; then
    mkdir -p -- "$out_dir" || {
        print -u2 -- "Error: could not create output directory '$out_dir'."
        exit 1
    }
fi

# ---- Reader (transparently handle gzip) --------------------------------------
# gzip -cd works on both macOS and Linux (zcat/gzcat differ between them).
if [[ $GTF == *.gz ]]; then
    reader=(gzip -cd --)
else
    reader=(cat --)
fi

# ---- Extract TSS -> sorted, de-duplicated BED --------------------------------
# NOTE: uses the POSIX 2-arg match()/substr()/sub() form so it runs under macOS
# (BSD) awk as well as GNU awk. The 3-arg match($0, /re/, arr) form is a
# gawk-only extension and is NOT used here. The feature type is passed in via
# -v; the ID attribute (gene_id vs transcript_id) is chosen to match it.
tmp="$BED.tmp.$$"
TMP_FILES+=("$tmp")

# -F'\t'

"${reader[@]}" "$GTF" \
  | awk  -v feature="$FEATURE" 'BEGIN { OFS = "\t" } 
         $3 == feature {
             if ($7 == "+") { s = $4 - 1; e = $4 }
             else           { s = $5 - 1; e = $5 }
             id = "."
             if (feature == "gene") {
                 if (match($0, /gene_id "[^"]+"/)) {
                     id = substr($0, RSTART, RLENGTH)
                     sub(/^gene_id "/, "", id)
                     sub(/"$/, "", id)
                 }
             } else {
                 if (match($0, /transcript_id "[^"]+"/)) {
                     id = substr($0, RSTART, RLENGTH)
                     sub(/^transcript_id "/, "", id)
                     sub(/"$/, "", id)
                 }
             }
             print $1, s, e, id, ".", $7
         }' \
  | LC_ALL=C sort -k1,1 -k2,2n -u > "$tmp"
rc=$?                              # PIPE_FAIL: nonzero if any stage failed

if (( rc != 0 )); then
    print -u2 -- "Error: TSS extraction failed (exit status $rc)."
    exit $rc                       # trap removes the temp file
fi

mv -- "$tmp" "$BED"

n=$(wc -l < "$BED" | tr -d ' ')
print -- "Wrote $n TSS record(s) ($FEATURE-level) to '$BED'."

# ---- Optional: expand TSS into promoter windows (bedtools slop) --------------
if [[ -e $CHROM_SIZES ]]; then
    if ! command -v bedtools >/dev/null 2>&1; then
        print -u2 -- "Error: bedtools not found in PATH; cannot run slop."
        print -u2 -- "       TSS BED '$BED' was written; install bedtools to enable slop."
        exit 1
    fi

    # Build a readable suffix from the flank size (1000 -> 1kb, 500 -> 500bp).
    if (( FLANK % 1000 == 0 )); then
        flank_label="$((FLANK / 1000))kb"
    else
        flank_label="${FLANK}bp"
    fi
    SLOP_BED="${BED%.bed}_${flank_label}.bed"
    slop_tmp="$SLOP_BED.tmp.$$"
    TMP_FILES+=("$slop_tmp")

    # Choose slop's input: either the full TSS BED, or a filtered copy with
    # edge-proximal TSSs removed (default). Filtering keeps a TSS only if
    # (start - FLANK) >= 0 AND (end + FLANK) <= chrom length, reading lengths
    # from $CHROM_SIZES. Records on chromosomes absent from chrom.sizes are
    # dropped (no length to validate; slop would otherwise error on them).
    slop_input="$BED"
    if (( DROP_EDGE )); then
        inbounds_tmp="$BED.inbounds.$$"
        TMP_FILES+=("$inbounds_tmp")

        awk -v f="$FLANK" 'BEGIN { OFS = "\t" }
             FNR == NR { len[$1] = $2; next }
             ($1 in len) && ($2 - f) >= 0 && ($3 + f) <= len[$1]' \
             "$CHROM_SIZES" "$BED" > "$inbounds_tmp"
        rc=$?
        if (( rc != 0 )); then
            print -u2 -- "Error: in-bounds TSS filtering failed (exit status $rc)."
            exit $rc               # trap removes the temp file
        fi

        kept=$(wc -l < "$inbounds_tmp" | tr -d ' ')
        dropped=$(( n - kept ))
        print -- "Filtered TSSs: kept $kept, dropped $dropped within ${FLANK} bp of a chromosome end (or on an unknown chromosome)."
        slop_input="$inbounds_tmp"
    else
        print -- "Edge filtering disabled (--no-drop); slop will clamp edge windows to chromosome bounds."
    fi

    bedtools slop -b "$FLANK" -i "$slop_input" -g "$CHROM_SIZES" > "$slop_tmp"
    rc=$?
    if (( rc != 0 )); then
        print -u2 -- "Error: bedtools slop failed (exit status $rc)."
        exit $rc                   # trap removes the temp file
    fi
    mv -- "$slop_tmp" "$SLOP_BED"

    m=$(wc -l < "$SLOP_BED" | tr -d ' ')
    print -- "Wrote $m window(s) (+/- ${FLANK} bp) to '$SLOP_BED'."
else
    print -u2 -- "Note: chrom.sizes '$CHROM_SIZES' not found; skipping bedtools slop."
    print -u2 -- "      Pass a chrom.sizes as the 3rd argument (or set CHROM_SIZES) to enable it."
fi
