#!/usr/bin/env bash
#
# run_trf.sh — Run Tandem Repeat Finder (TRF) on all .fa files in a directory
#
# Usage:
#   run_trf.sh -s SEQUENCE_DIR -o OUTPUT_DIR -l LOG_DIR
#
# Options:
#   -s  Directory containing .fa sequence files.
#   -o  Directory for TRF output files.
#   -l  Directory for log files.
#   -h  Show this help message and exit.
#
# Default TRF parameters:
#   MATCH=2, MISMATCH=7, DEL_INS=7, MATCH_SCORE=80
#   MIN_REPEAT=10, MAX_PERIOD=50, HIT_THRESHOLD=500
#

usage(){
  cat <<EOF
Usage: $(basename "$0") -s SEQUENCE_DIR -o OUTPUT_DIR -l LOG_DIR

Options:
  -s  Directory containing .fa sequence files.
  -o  Directory for TRF output files.
  -l  Directory for log files.
  -h  Show this help message and exit.

Default TRF parameters:
  MATCH=2, MISMATCH=7, DEL_INS=7, MATCH_SCORE=80
  MIN_REPEAT=10, MAX_PERIOD=50, HIT_THRESHOLD=500
EOF
}

#load required modules
module load trf/4.09.1

# Parse CLI args
while getopts ":s:o:l:h" opt; do
  case $opt in
    s) SEQUENCE_DIR=$OPTARG ;;
    o) OUTPUT_DIR=$OPTARG ;;
    l) LOG_DIR=$OPTARG ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done

# Check required args
if [[ -z ${SEQUENCE_DIR:-} || -z ${OUTPUT_DIR:-} || -z ${LOG_DIR:-} ]]; then
  echo "Error: missing required arguments." >&2
  usage
  exit 1
fi

# Verify TRF is available
if ! command -v trf &>/dev/null; then
  echo "Error: 'trf' not found in PATH. Install TRF or load its module." >&2
  exit 1
fi

# TRF parameters
MATCH=2
MISMATCH=7
DEL_INS=7
MATCH_SCORE=80
MIN_REPEAT=10
MAX_PERIOD=50
HIT_THRESHOLD=500

# Make output dirs
mkdir -p "$LOG_DIR" "$OUTPUT_DIR"

# Loop over FASTA files
shopt -s nullglob
files=("$SEQUENCE_DIR"/*.fa)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "Error: no .fa files found in $SEQUENCE_DIR" >&2
  exit 1
fi

for SEQ in "${files[@]}"; do
  BASE=$(basename "$SEQ" .fa)
  LOG_FILE="$LOG_DIR/${BASE}_trf.log"
  OUT_FILE="$OUTPUT_DIR/${BASE}_trf.txt"

  echo "Running TRF on $SEQ" | tee -a "$LOG_FILE"
  trf "$SEQ" \
      $MATCH $MISMATCH $DEL_INS $MATCH_SCORE $MIN_REPEAT $MAX_PERIOD $HIT_THRESHOLD \
      -f -d -m \
    > "$OUT_FILE" 2>> "$LOG_FILE"

  if [[ $? -ne 0 ]]; then
    echo "ERROR: TRF failed on $SEQ" | tee -a "$LOG_FILE" >&2
    continue
  fi

  echo "Completed TRF for $SEQ" | tee -a "$LOG_FILE"
done

echo "All done."
