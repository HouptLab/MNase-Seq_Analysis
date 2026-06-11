#!/bin/zsh
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <source_directory> <output directory>"
    exit 1
fi


SRCDIR="$1"

OUTDIR="${2:-bowtie}" 

# Check directory exists
if [[ ! -d "$SRCDIR" ]]; then
    echo "Error: Directory '$SRCDIR' not found."
    exit 1
fi

# Put all fastq.gz files into array
files=("$SRCDIR"/*.fastq.gz)
printf '%s\n' "${files[@]}"

mkdir -p "$OUTDIR"



LOG="$OUTDIR/bowtie.log"

echo "$LOG"

# Backup existing log if present
if [[ -f "$LOG" ]]; then
    mv "$LOG" "$OUTDIR/$LOG_$(date '+%Y%m%d_%H%M%S').log"
    echo "Existing log backed up."
fi