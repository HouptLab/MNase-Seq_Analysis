#!/usr/bin/env zsh
#
# trim_pe.zsh — run Trimmomatic PE on every raw R1 read file in a directory.
#               Always runs in the background; the prompt returns immediately.
#
# Usage:
#   ./trim_pe.zsh <source_directory> [log_file]
#
# Follow progress with:
#   tail -f <source_directory>/trimming_B.log
#

set -u  # treat unset variables as an error

# --- configuration -------------------------------------------------------
TRIMMOMATIC_JAR="/Applications/Trimmomatic-0.40/trimmomatic-0.40.jar"
ADAPTERS="TruSeq3-PE.fa"        # resolved relative to the source dir
THREADS=20
R1_GLOB="*_R1_001.fastq.gz"    # RAW R1 inputs only — excludes _paired_/_unpaired_ outputs

# Resolve a real Java runtime explicitly. A non-interactive / nohup'd shell does
# NOT source ~/.zshrc, so a JDK added to PATH there is invisible and a bare
# `java` falls back to the macOS /usr/bin/java stub (silent exit 1). Prefer
# java_home, which returns a real JDK or nothing. Override with: JAVA=/path ./trim_pe.zsh ...
JAVA="${JAVA:-}"
if [[ -z "$JAVA" ]]; then
    if [[ -x /usr/libexec/java_home ]] && _jh="$(/usr/libexec/java_home 2>/dev/null)" && [[ -n "$_jh" ]]; then
        JAVA="$_jh/bin/java"
    elif command -v java >/dev/null 2>&1; then
        JAVA="$(command -v java)"
    fi
fi

# --- argument handling (foreground, so errors are visible immediately) ---
if (( $# < 1 )); then
    print -u2 "Usage: $0 <source_directory> [log_file]"
    exit 1
fi

SRC_DIR="$1"
LOG_FILE="${2:-$SRC_DIR/trimming.log}"

if [[ ! -d "$SRC_DIR" ]]; then
    print -u2 "Error: '$SRC_DIR' is not a directory."
    exit 1
fi

if [[ -z "$JAVA" || ! -x "$JAVA" ]]; then
    print -u2 "Error: no Java runtime found. Install a JDK, or run with an explicit path:"
    print -u2 "       JAVA=/opt/homebrew/opt/openjdk/bin/java $0 $*"
    exit 1
fi

if [[ ! -f "$TRIMMOMATIC_JAR" ]]; then
    print -u2 "Error: Trimmomatic jar not found at $TRIMMOMATIC_JAR"
    exit 1
fi

# Make the log path absolute so it is unaffected by the cd below.
LOG_FILE="${LOG_FILE:A}"

# --- self-background -----------------------------------------------------
# On first invocation, relaunch detached (nohup survives terminal/SSH close).
# The relaunched copy carries TRIM_PE_BACKGROUND=1 (and the resolved JAVA path,
# so the background shell doesn't have to re-find it) and skips this block.
if [[ -z "${TRIM_PE_BACKGROUND:-}" ]]; then
    # Rotate an existing log to a timestamped backup so the previous run's
    # log is preserved and this run starts with a fresh trimming.log.
    if [[ -f "$LOG_FILE" ]]; then
        BACKUP="${LOG_FILE:r}_$(date '+%Y%m%d-%H%M%S').${LOG_FILE:e}"
        mv "$LOG_FILE" "$BACKUP"
        print "trim_pe: backed up previous log -> $BACKUP"
    fi
    TRIM_PE_BACKGROUND=1 JAVA="$JAVA" nohup "$0" "$@" </dev/null >/dev/null 2>&1 &
    print "trim_pe: running in background (PID $!)."
    print "trim_pe: using java -> $JAVA"
    print "trim_pe: follow progress with ->  tail -f \"$LOG_FILE\""
    exit 0
fi

# --- main work (background instance only) --------------------------------
cd "$SRC_DIR" || exit 1
setopt NULL_GLOB   # empty glob -> empty list instead of an error

log() { print -r -- "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >> "$LOG_FILE" }

log "=== Trimmomatic run started in $PWD (PID $$) ==="
log "java: $JAVA"

for i in ${~R1_GLOB}; do
    log "START  $i"
    "$JAVA" -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" -phred33 \
        "$i" "${i/_R1_/_R2_}" \
        "${i/_R1_/_R1_paired_}" "${i/_R1_/_R1_unpaired_}" \
        "${i/_R1_/_R2_paired_}" "${i/_R1_/_R2_unpaired_}" \
        ILLUMINACLIP:${ADAPTERS}:2:30:10:1:TRUE \
        MINLEN:25 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
        </dev/null >> "$LOG_FILE" 2>&1
    rc=$?
    log "FINISH $i (exit $rc)"
done

log "=== All files processed ==="
