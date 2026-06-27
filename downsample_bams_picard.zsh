#!/usr/bin/env zsh
set -euo pipefail

# ===========================================================================
# downsample_bams_picard.zsh
# Downsample a directory of paired-end BAM files to the lowest template
# (fragment) count among them, using Picard DownsampleSam. Picard keeps
# all reads of a template (both mates + secondary/supplementary) together,
# so mate pairs are never split. Designed for macOS (Apple silicon).
#
# Usage:
#   ./downsample_bams_picard.zsh [INPUT_DIR] [OUTPUT_DIR]
#   PICARD_JAR=/path/picard.jar SEED=42 STRATEGY=HighAccuracy \
#       ./downsample_bams_picard.zsh /path/to/bams /path/to/output
#
# Defaults:
#   INPUT_DIR  = current dir
#   OUTPUT_DIR = ./downsampled
#   SEED       = 42
#   STRATEGY   = HighAccuracy   (better adherence to target proportion)
#   ACCURACY   = 0.0001
#   PICARD_JAR = /Applications/picard.jar
#
# Fragment (template) counting uses:
#   samtools view -c -f 0x40 -F 0x90C
#     0x40  = first in pair        (count each template once)
#     0x900 = exclude secondary(0x100)+supplementary(0x800)
#     0x00C = exclude read-unmapped(0x4)+mate-unmapped(0x8)
#   -> 0x90C = 0x900 + 0x00C
# ===========================================================================

# ---------------------------------------------------------------------------
# Directories / arguments
# ---------------------------------------------------------------------------
INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-./downsampled}"
SEED="${SEED:-42}"
STRATEGY="${STRATEGY:-HighAccuracy}"
ACCURACY="${ACCURACY:-0.0001}"
PICARD_JAR="${PICARD_JAR:-/Applications/picard.jar}"
LOGFILE="${OUTPUT_DIR}/downsample_picard_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# Logging helper: write to both stdout and the logfile, with timestamp
# ---------------------------------------------------------------------------
log() {
    print -r -- "[$(date '+%Y-%m-%d %H:%M:%S %Z')] $*" | tee -a "$LOGFILE"
}

# ---------------------------------------------------------------------------
# 1. Check dependencies: samtools, java, picard.jar
# ---------------------------------------------------------------------------
if ! command -v samtools >/dev/null 2>&1; then
    print -r -- "ERROR: samtools not found in PATH. Install it (e.g. 'brew install samtools') and retry." >&2
    exit 1
fi

if ! command -v java >/dev/null 2>&1; then
    print -r -- "ERROR: java not found in PATH. Install a JDK/JRE (e.g. 'brew install temurin') and retry." >&2
    exit 1
fi

if [[ ! -f "$PICARD_JAR" ]]; then
    print -r -- "ERROR: picard.jar not found at '${PICARD_JAR}'. Set PICARD_JAR to its path and retry." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Threads from core count (Apple silicon, macOS)
#    Picard/DownsampleSam is largely single-threaded; THREADS is used for the
#    samtools counting/indexing steps.
# ---------------------------------------------------------------------------
CORES=$(sysctl -n hw.ncpu)
THREADS="${THREADS:-$CORES}"

# ---------------------------------------------------------------------------
# 5. Record versions  +  4. start timestamp
# ---------------------------------------------------------------------------
log "=== Picard downsample pipeline START ==="
log "Host: $(hostname)  |  macOS $(sw_vers -productVersion)  |  arch: $(uname -m)"
log "Logical cores detected: ${CORES}  |  THREADS set to: ${THREADS}"
log "Input dir:  ${INPUT_DIR}"
log "Output dir: ${OUTPUT_DIR}"
log "RNG seed:   ${SEED}   |   Strategy: ${STRATEGY}   |   Accuracy: ${ACCURACY}"
log "picard.jar: ${PICARD_JAR}"
log "samtools path: $(command -v samtools)"
log "samtools version:"
samtools --version | tee -a "$LOGFILE"
log "java version:"
java -version 2>&1 | tee -a "$LOGFILE"
log "picard version:"
java -jar "$PICARD_JAR" DownsampleSam --version 2>&1 | tee -a "$LOGFILE" || true

# ---------------------------------------------------------------------------
# Gather BAM files
# ---------------------------------------------------------------------------
bams=("$INPUT_DIR"/*.bam(N))   # (N) = nullglob: empty array if no matches
if (( ${#bams[@]} == 0 )); then
    log "ERROR: No BAM files found in ${INPUT_DIR}"
    exit 1
fi

# ---------------------------------------------------------------------------
# 6. Count templates (fragments) in each BAM; find target_min
#    -f 0x40 -F 0x90C  -> first-in-pair, primary, both mates mapped
# ---------------------------------------------------------------------------
log "--- Counting templates (samtools view -c -f 0x40 -F 0x90C) ---"

typeset -A counts
target_min=""
for bam in "${bams[@]}"; do
    c=$(samtools view -c -f 0x40 -F 0x90C -@ "$THREADS" "$bam")
    counts[$bam]=$c
    log "  ${bam:t}: ${c} templates"
    if [[ -z "$target_min" || $c -lt $target_min ]]; then
        target_min=$c
    fi
done
log "target_min = ${target_min} templates"

if [[ "$target_min" -eq 0 ]]; then
    log "ERROR: target_min is 0 templates; check input BAMs and flag filters."
    exit 1
fi

# ---------------------------------------------------------------------------
# 7. Downsample each BAM to target_min templates with Picard DownsampleSam.
#    PROBABILITY (P) = target_min / this_file's_template_count, in [0,1].
#    Reads of a template travel together, so mates stay intact.
# ---------------------------------------------------------------------------
log "--- Downsampling (Picard DownsampleSam) ---"

for bam in "${bams[@]}"; do
    c=${counts[$bam]}
    out="${OUTPUT_DIR}/${bam:t:r}.downsampled.bam"

    if [[ $c -eq $target_min ]]; then
        # Already at the minimum; copy through unchanged (no subsampling).
        log "  ${bam:t}: at target_min, copying without subsampling"
        samtools view -b -@ "$THREADS" "$bam" > "$out" 2>>"$LOGFILE"
    else
        prob=$(awk -v t="$target_min" -v c="$c" 'BEGIN { printf "%.6f", t / c }')
        log "  ${bam:t}: ${c} -> ~${target_min} templates (P=${prob})"
        java -jar "$PICARD_JAR" DownsampleSam \
            I="$bam" \
            O="$out" \
            STRATEGY="$STRATEGY" \
            PROBABILITY="$prob" \
            ACCURACY="$ACCURACY" \
            RANDOM_SEED="$SEED" \
            >>"$LOGFILE" 2>&1
    fi

    # Report resulting template count for verification
    rc=$(samtools view -c -f 0x40 -F 0x90C -@ "$THREADS" "$out")
    log "    -> wrote ${out:t} (${rc} templates)"
done

# ---------------------------------------------------------------------------
# Final step: index all downsampled BAM files
# ---------------------------------------------------------------------------
log "--- Indexing downsampled BAM files (samtools index) ---"
for out in "$OUTPUT_DIR"/*.downsampled.bam(N); do
    samtools index -@ "$THREADS" "$out" 2>>"$LOGFILE"
    log "  indexed: ${out:t} -> ${out:t}.bai"
done

# ---------------------------------------------------------------------------
# 4. Completion timestamp
# ---------------------------------------------------------------------------
log "=== Picard downsample pipeline COMPLETE ==="
