#!/usr/bin/env bash
#
# run_rrna_blast.sh — BLAST rRNA query against all genome FASTAs
#
# Usage:
#   run_rrna_blast.sh -q QUERY_FASTA -g GENOME_DIR -o OUT_DIR -s SUMMARY_FILE [-l LOG_DIR] [-c CPUS]
#
# Options:
#   -q QUERY_FASTA   Path to rRNA query FASTA (e.g. rRNA_seq.fasta)
#   -g GENOME_DIR    Directory containing genome FASTAs (*.fa)
#   -o OUT_DIR       Directory for per-species BLAST outputs
#   -s SUMMARY_FILE  TSV file to append summary results
#   -l LOG_DIR       Directory for log files (default: ./logs)
#   -c CPUS          Number of threads for blastn (default: detected cores)
#   -h               Show this help message and exit

usage(){
  cat <<EOF
Usage: $(basename "$0") -q QUERY_FASTA -g GENOME_DIR -o OUT_DIR -s SUMMARY_FILE [-l LOG_DIR] [-c CPUS]

Options:
  -q QUERY_FASTA   Path to rRNA query FASTA (required)
  -g GENOME_DIR    Directory containing genome FASTAs (*.fa) (required)
  -o OUT_DIR       Directory for BLAST outputs (required)
  -s SUMMARY_FILE  TSV file to append summary (required)
  -l LOG_DIR       Directory for logs (default: ./logs)
  -c CPUS          Threads for blastn (default: detected cores)
  -h               Show this help and exit
EOF
}

# Defaults
LOG_DIR="./logs"
CPUS=""

# Parse arguments
while getopts ":q:g:o:s:l:c:h" opt; do
  case $opt in
    q) QUERY_FASTA=$OPTARG ;;
    g) GENOME_DIR=$OPTARG ;;
    o) OUT_DIR=$OPTARG ;;
    s) SUMMARY_FILE=$OPTARG ;;
    l) LOG_DIR=$OPTARG ;;
    c) CPUS=$OPTARG ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done

# Load required modules
module load gcc/8.2.0
module load perl/5.28.0
module load blast+

# Check required
if [[ -z "${QUERY_FASTA:-}" || -z "${GENOME_DIR:-}" || -z "${OUT_DIR:-}" || -z "${SUMMARY_FILE:-}" ]]; then
  echo "Error: -q, -g, -o, and -s are required." >&2
  usage
  exit 1
fi

# Verify files/directories
if [[ ! -f "$QUERY_FASTA" ]]; then
  echo "Error: query FASTA not found: $QUERY_FASTA" >&2
  exit 1
fi
if [[ ! -d "$GENOME_DIR" ]]; then
  echo "Error: genome directory not found: $GENOME_DIR" >&2
  exit 1
fi

# Determine CPUS
if [[ -z "$CPUS" ]]; then
  if command -v nproc &>/dev/null; then
    CPUS=$(nproc)
  else
    CPUS=1
  fi
fi

# Create directories
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Initialize summary file with header if missing
if [[ ! -f "$SUMMARY_FILE" ]]; then
  printf "species\tqseqid\tsseqid\tsstart\tsend\tpident\tlength\n" > "$SUMMARY_FILE"
fi

# Loop through genome FASTAs
shopt -s nullglob
for GENOME in "$GENOME_DIR"/*.fa; do
  species=$(basename "$GENOME" .fa)
  LOG_FILE="$LOG_DIR/${species}_rrna_blast.log"
  DB_PREFIX="$OUT_DIR/${species}_db"
  BLAST_OUT="$OUT_DIR/${species}_rrna_hits.blastn"
  TMP_HITS="$OUT_DIR/${species}_rrna_tmp_hits.tsv"

  echo "=== [$species] $(date +'%Y-%m-%d %H:%M:%S') ===" | tee "$LOG_FILE"

  # Build BLAST DB
  echo "Building BLAST DB..." | tee -a "$LOG_FILE"
  makeblastdb -in "$GENOME" -dbtype nucl -out "$DB_PREFIX" \
    >> "$LOG_FILE" 2>&1

  # Run BLASTN
  echo "Running blastn of $QUERY_FASTA vs. $species..." | tee -a "$LOG_FILE"
  blastn \
    -query "$QUERY_FASTA" \
    -db "$DB_PREFIX" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen sstart send evalue bitscore" \
    -evalue 1e-10 \
    -num_threads "$CPUS" \
    -out "$BLAST_OUT" \
    >> "$LOG_FILE" 2>&1

  # Extract contigs with hits
  cut -f2 "$BLAST_OUT" | sort -u > "$TMP_HITS"
  if [[ ! -s "$TMP_HITS" ]]; then
    echo "No hits found for $species" | tee -a "$LOG_FILE"
    continue
  fi

  hit_count=$(wc -l < "$TMP_HITS")
  echo "Found $hit_count contigs with hits" | tee -a "$LOG_FILE"

  # Append to summary
  awk -v sp="$species" 'BEGIN {OFS="\t"}
    { print sp, $1, $2, $7, $8, $3, $4 }
  ' "$BLAST_OUT" >> "$SUMMARY_FILE"

  echo "Summary updated: $SUMMARY_FILE" | tee -a "$LOG_FILE"
done

echo "All done! Results in $OUT_DIR/, logs in $LOG_DIR/, summary in $SUMMARY_FILE."
