#!/bin/zsh
#
# run_bowtie2.zsh  (macOS / zsh)
# Align paired-end FASTQ files to GRCr8 with bowtie2.
#
# Usage:
#   nohup ./run_bowtie2.zsh <source_directory> <destination_directory> &
#
# Outputs (SAM files, per-sample bowtie2 logs, and bowtie.log) are written
# to <destination_directory> (created if it does not exist). The bowtie2
# index path, thread count, and concurrency cap are configured below.

setopt PIPE_FAIL       # a pipeline fails if any stage fails
setopt NULL_GLOB       # a non-matching glob expands to nothing (no error)

# ---- Argument handling -------------------------------------------------------
if [[ $# -ne 2 ]]; then
    print -u2 -- "Usage: $0 <source_directory> <destination_directory>"
    exit 1
fi

SRC_DIR=$1
DEST_DIR=$2

if [[ ! -d $SRC_DIR ]]; then
    print -u2 -- "Error: source '$SRC_DIR' is not a directory."
    exit 1
fi

# Create the destination directory if needed.
if [[ ! -d $DEST_DIR ]]; then
    mkdir -p -- "$DEST_DIR" || {
        print -u2 -- "Error: could not create destination '$DEST_DIR'."
        exit 1
    }
fi

# ---- Configuration -----------------------------------------------------------
# Directory containing this script (absolute, symlinks resolved). The bowtie2
# index is looked up relative to here, not the current working directory.
SCRIPT_DIR=${0:A:h}

INDEX="$SCRIPT_DIR/bowtie2_indexes/GRCr8"
THREADS=8
MAX_JOBS=4        # max bowtie2 jobs to run concurrently (MAX_JOBS * THREADS <= cores)
PATTERN="*_R1_paired*.fastq.gz"
LOGFILE="$DEST_DIR/bowtie.log"

# ---- Log file rotation -------------------------------------------------------
# If bowtie.log already exists, rename it with a timestamp and start fresh.
if [[ -e $LOGFILE ]]; then
    mv -- "$LOGFILE" "${LOGFILE%.log}.$(date '+%Y%m%d_%H%M%S').log"
fi
: > "$LOGFILE"

# ---- Logging helper ----------------------------------------------------------
log() {
    print -r -- "$(date '+%Y-%m-%d %H:%M:%S') $*" >> "$LOGFILE"
}

# ---- Index resolution --------------------------------------------------------
# Accept either a bowtie2 basename prefix (…/GRCr8 with …/GRCr8.1.bt2 alongside)
# or a directory that contains the index files, and print the usable basename.
# Returns non-zero if no .bt2/.bt2l index files can be found.
resolve_index() {
    local base=$1
    while [[ $base == */ && $base != / ]]; do base=${base%/}; done   # trim trailing /

    # Case 1: already a valid basename prefix.
    if [[ -e ${base}.1.bt2 || -e ${base}.1.bt2l ]]; then
        print -r -- "$base"
        return 0
    fi

    # Case 2: a directory holding the index files -> derive the prefix.
    if [[ -d $base ]]; then
        local hits=( ${base}/*.1.bt2(N) ${base}/*.1.bt2l(N) )
        if (( ${#hits[@]} > 0 )); then
            local f=${hits[1]}
            f=${f%.1.bt2l}
            f=${f%.1.bt2}
            print -r -- "$f"
            return 0
        fi
    fi

    return 1
}

# ---- Main --------------------------------------------------------------------
log "=== Run started ==="
log "Source directory: $SRC_DIR"
log "Destination directory: $DEST_DIR"

# Resolve the bowtie2 index basename before launching anything.
INDEX_INPUT=$INDEX
if ! INDEX=$(resolve_index "$INDEX"); then
    msg="Error: no bowtie2 index (.bt2/.bt2l) found at or under '$INDEX_INPUT'."
    print -u2 -- "$msg"
    log "$msg"
    log "=== Run aborted ==="
    exit 1
fi

log "Index: $INDEX  Threads: $THREADS  Max concurrent: $MAX_JOBS  Pattern: $PATTERN"

# Collect matching R1 files. ${~PATTERN} forces glob expansion of the variable;
# NULL_GLOB makes an empty match yield an empty array instead of an error.
files=("$SRC_DIR"/${~PATTERN})

# Track jobs in flight (FIFO queue) and remember each PID's sample name.
typeset -A sample_of
running=()        # PIDs currently executing
fail=0
launched=0

# Wait for the oldest running job, record its status, drop it from the queue.
reap_one() {
    local pid=${running[1]}        # zsh arrays are 1-indexed
    running=(${running[2,-1]})     # remove the first element
    if wait $pid; then
        log "PID $pid (sample '${sample_of[$pid]}') completed successfully."
    else
        local status=$?
        log "PID $pid (sample '${sample_of[$pid]}') FAILED with exit status $status."
        (( fail++ ))
    fi
}

for i in "${files[@]}"; do
    r2=${i/R1/R2}

    if [[ ! -e $r2 ]]; then
        log "WARNING: mate file not found for '$i' (expected '$r2') -- skipping."
        continue
    fi

    # Derive a per-sample name so concurrent jobs don't clobber one SAM file.
    sample=${i:t}                  # :t = tail (basename)
    sample=${sample%%_R1_paired*}
    sam="$DEST_DIR/${sample}_GRCr8.1_alignment.sam"
    sample_log="$DEST_DIR/${sample}.bowtie2.log"

    # Throttle: if at capacity, wait for the oldest job before launching more.
    while (( ${#running[@]} >= MAX_JOBS )); do
        reap_one
    done

    log "Launching alignment for sample '$sample' -> $sam (bowtie2 output: $sample_log)"

    nohup bowtie2 -x "$INDEX" -1 "$i" -2 "$r2" -S "$sam" -p "$THREADS" \
        > "$sample_log" 2>&1 &

    pid=$!
    running+=($pid)
    sample_of[$pid]=$sample
    (( launched++ ))
    log "  PID $pid : R1='$i'  R2='$r2'  (in flight: ${#running[@]}/$MAX_JOBS)"
done

if (( launched == 0 )); then
    log "No files matching '$PATTERN' found in '$SRC_DIR'. Nothing to do."
    log "=== Run finished ==="
    exit 0
fi

log "All $launched job(s) launched. Waiting for remaining jobs to finish..."

# Drain any jobs still running.
while (( ${#running[@]} > 0 )); do
    reap_one
done

log "All jobs finished. Failures: $fail"
log "=== Run finished ==="

exit $(( fail > 0 ? 1 : 0 ))
