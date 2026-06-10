#!/bin/zsh
#
# do_trimming.zsh — run Trimmomatic PE on every *_R1* read file in a directory.
#               Always runs in the background; the prompt returns immediately.
#
# Usage:
#   ./do_trimming.sh <source_directory> [log_file]
#
# Examples:
#   ./trim_pe.zsh /data/run_B
#   ./trim_pe.zsh /data/run_B /data/run_B/trimming_B.log
#
# Follow progress with:
#   tail -f /data/run_B/trimming_B.log
#

set -u  # treat unset variables as an error

# --- configuration -------------------------------------------------------
TRIMMOMATIC_JAR="/Applications/Trimmomatic-0.40/trimmomatic-0.40.jar"
ADAPTERS="TruSeq3-PE.fa"   # resolved relative to the source dir (see note below)
THREADS=20

# --- argument handling (runs in the foreground, so errors are visible) ---
if (( $# < 1 )); then
    print -u2 "Usage: $0 <source_directory> [log_file]"
    exit 1
fi

SRC_DIR="$1"
LOG_FILE="${2:-$SRC_DIR/trimming_B.log}"

if [[ ! -d "$SRC_DIR" ]]; then
    print -u2 "Error: '$SRC_DIR' is not a directory."
    exit 1
fi

if [[ ! -f "$TRIMMOMATIC_JAR" ]]; then
    print -u2 "Error: Trimmomatic jar not found at $TRIMMOMATIC_JAR"
    exit 1
fi

# Make the log path absolute so it is unaffected by the cd below.
LOG_FILE="${LOG_FILE:A}"

# --- self-background -----------------------------------------------------
# On first invocation, relaunch detached (nohup survives terminal/SSH close),
# report the PID + log path, and return the prompt. The relaunched copy carries
# TRIM_PE_BACKGROUND=1 so it skips this block and does the actual work.
if [[ -z "${TRIM_PE_BACKGROUND:-}" ]]; then
    TRIM_PE_BACKGROUND=1 nohup "$0" "$@" >/dev/null 2>&1 &
    print "trim_pe: running in background (PID $!)."
    print "trim_pe: follow progress with ->  tail -f \"$LOG_FILE\""
    exit 0
fi

# --- main work (background instance only) --------------------------------
cd "$SRC_DIR" || exit 1

# zsh errors on a non-matching glob by default; NULL_GLOB makes it expand
# to nothing instead so an empty directory doesn't abort the script.
setopt NULL_GLOB

# Timestamped line appended to the log.
log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "$LOG_FILE" }

log "=== Trimmomatic run started in $PWD (PID $$) ==="

for i in *_R1*; do
    log "START  $i"
    java -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" -phred33 \
        "$i" "${i/R1/R2}" \
        "${i/R1/R1_paired}" "${i/R1/R1_unpaired}" \
        "${i/R1/R2_paired}" "${i/R1/R2_unpaired}" \
        ILLUMINACLIP:${ADAPTERS}:2:30:10:1:TRUE \
        MINLEN:25 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
        >> "$LOG_FILE" 2>&1
    rc=$?
    log "FINISH $i (exit $rc)"
done

log "=== All files processed ==="
