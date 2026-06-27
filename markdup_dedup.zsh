#!/usr/bin/env zsh
set -euo pipefail

# ===========================================================================
# markdup_dedup.zsh
# Remove PCR duplicates from BAM files using samtools.
# Designed for macOS (Apple silicon, e.g. Mac Studio).
#
# Usage:
#   ./markdup_dedup.zsh [INPUT_DIR] [OUTPUT_DIR]
#   THREADS=8 ./markdup_dedup.zsh /path/to/bams /path/to/output
#
# Defaults: INPUT_DIR=current dir, OUTPUT_DIR=./dedup
# Override thread count with the THREADS environment variable.
# ===========================================================================

# ---------------------------------------------------------------------------
# Directories / arguments
# ---------------------------------------------------------------------------
INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-./dedup}"
LOGFILE="${OUTPUT_DIR}/markdup_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# Logging helper: write to both stdout and the logfile, with timestamp
# ---------------------------------------------------------------------------
log() {
    print -r -- "[$(date '+%Y-%m-%d %H:%M:%S %Z')] $*" | tee -a "$LOGFILE"
}

# ---------------------------------------------------------------------------
# 1. Check samtools is installed
# ---------------------------------------------------------------------------
if ! command -v samtools >/dev/null 2>&1; then
    print -r -- "ERROR: samtools not found in PATH. Install it (e.g. 'brew install samtools') and retry." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Threads from core count (Apple silicon, macOS)
#    sysctl -n hw.ncpu reports logical cores on macOS
# ---------------------------------------------------------------------------
CORES=$(sysctl -n hw.ncpu)
THREADS="${THREADS:-$CORES}"

# ---------------------------------------------------------------------------
# 5. Record samtools version  +  4. start timestamp
# ---------------------------------------------------------------------------
log "=== markdup pipeline START ==="
log "Host: $(hostname)  |  macOS $(sw_vers -productVersion)  |  arch: $(uname -m)"
log "Logical cores detected: ${CORES}  |  THREADS set to: ${THREADS}"
log "Input dir:  ${INPUT_DIR}"
log "Output dir: ${OUTPUT_DIR}"
log "samtools path: $(command -v samtools)"
log "samtools version:"
samtools --version | tee -a "$LOGFILE"

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for bam in "$INPUT_DIR"/*.bam; do
    [[ -e "$bam" ]] || { log "ERROR: No BAM files found in ${INPUT_DIR}"; exit 1; }

    sample="${bam:t:r}"
    log "--- Processing sample: ${sample} ---"

    tmp_namesort="$OUTPUT_DIR/${sample}.namesort.bam"
    tmp_fixmate="$OUTPUT_DIR/${sample}.fixmate.bam"
    tmp_possort="$OUTPUT_DIR/${sample}.possort.bam"
    final="$OUTPUT_DIR/${sample}.dedup.bam"

    # markdup requires fixmate (-m) on name-sorted input, then position sort
    samtools sort  -n -@ "$THREADS" -o "$tmp_namesort" "$bam"        2>>"$LOGFILE"
    samtools fixmate -m -@ "$THREADS"   "$tmp_namesort" "$tmp_fixmate" 2>>"$LOGFILE"
    samtools sort     -@ "$THREADS" -o "$tmp_possort"   "$tmp_fixmate" 2>>"$LOGFILE"

    # -r remove duplicates, -s print stats (stats captured to logfile)
    samtools markdup -r -s -@ "$THREADS" "$tmp_possort" "$final"    2>>"$LOGFILE"
    samtools index "$final" 2>>"$LOGFILE"

    rm -f "$tmp_namesort" "$tmp_fixmate" "$tmp_possort"
    log "Completed: ${final}"
done

# ---------------------------------------------------------------------------
# 4. Completion timestamp
# ---------------------------------------------------------------------------
log "=== markdup pipeline COMPLETE ==="
