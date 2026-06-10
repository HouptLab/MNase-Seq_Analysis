#!/bin/zsh


# brew install fastqc
# brew install parallel
# One thing to be aware of — the first time you run it, GNU parallel shows a citation # notice asking you to run: parallel --citation. Type will cite to dismiss it, otherwise it # prints the notice every run. 
# You can also suppress it with --will-cite flag in your script:

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <source_directory> <output directory>"
    exit 1
fi


SRCDIR="$1"

OUTDIR="${2:-fastqc_raw}" 

# Check directory exists
if [[ ! -d "$SRCDIR" ]]; then
    echo "Error: Directory '$SRCDIR' not found."
    exit 1
fi

# Put all fastq.gz files into array
files=("$SRCDIR"/*.fastq.gz)
printf '%s\n' "${files[@]}"

mkdir -p "$OUTDIR"



LOG="$OUTDIR/fastqc_raw.log"

echo "$LOG"

# Backup existing log if present
if [[ -f "$LOG" ]]; then
    mv "$LOG" "$OUTDIR/$LOG_$(date '+%Y%m%d_%H%M%S').log"
    echo "Existing log backed up."
fi




# # Iterate and run fastqc
# for i in "${files[@]}"; do
#     echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastqc Processing: $i" | tee -a "$LOG"
#     fastqc "$i" -t 15 -o "$OUTDIR"/ >> "$LOG" 2>&1
# done 

# use parallel


# Auto-detect cores and divide among jobs
CORES=$(sysctl -n hw.logicalcpu)
JOBS=4
THREADS=$(( CORES / JOBS ))


echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running $JOBS parallel jobs, $THREADS threads each (of $CORES total cores)" | tee -a "$LOG"


FASTQC_CMD="fastqc {} -t $THREADS -o $OUTDIR/ 2>&1 | tee -a $LOG"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running $JOBS parallel jobs, $THREADS threads each (of $CORES total cores)" | tee -a "$LOG"

printf '%s\n' "${files[@]}" | parallel --will-cite --group -j "$JOBS" "$FASTQC_CMD"
    

echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastqc Done." | tee -a "$LOG"