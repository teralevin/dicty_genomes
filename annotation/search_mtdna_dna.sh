#!/usr/bin/env bash
#
# search_mtdna_dna.sh — BLAST reference mtDNA against species contigs
#
# Usage:
#   search_mtdna_dna.sh -d MT_DNA -s SPECIES_DIR -o OUT_DIR [-l LOG_DIR] [-c CPUS]
#
# Options:
#   -d MT_DNA       Reference mitochondrial DNA FASTA (required)
#   -s SPECIES_DIR  Directory containing per-species contig FASTAs (*.fa) (required)
#   -o OUT_DIR      Directory for BLAST DBs & blastn results (required)
#   -l LOG_DIR      Directory for logs (default: ./logs)
#   -c CPUS         Threads for blastn (default: detected cores)
#   -h              Show this help and exit
#

usage() {
  cat <<EOF
Usage: $(basename "$0") -d MT_DNA -s SPECIES_DIR -o OUT_DIR [-l LOG_DIR] [-c CPUS]

Options:
  -d MT_DNA       Reference mitochondrial DNA FASTA (required)
  -s SPECIES_DIR  Directory containing per-species contig FASTAs (*.fa) (required)
  -o OUT_DIR      Directory for BLAST DBs & blastn results (required)
  -l LOG_DIR      Directory for logs (default: ./logs)
  -c CPUS         Threads for blastn (default: detected cores)
  -h              Show this help and exit
EOF
}

# Load required modules
module load gcc/8.2.0
module load perl/5.28.0
module load blast+

# Defaults
LOG_DIR="./logs"
CPUS=""

# Parse arguments
while getopts ":d:s:o:l:c:h" opt; do
  case $opt in
    d) MT_DNA=$OPTARG ;;
    s) SPECIES_DIR=$OPTARG ;;
    o) OUT_DIR=$OPTARG ;;
    l) LOG_DIR=$OPTARG ;;
    c) CPUS=$OPTARG ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done

# Check required args
if [[ -z "${MT_DNA:-}" || -z "${SPECIES_DIR:-}" || -z "${OUT_DIR:-}" ]]; then
  echo "Error: -d, -s, and -o are required." >&2
  usage
  exit 1
fi

# Verify inputs
if [[ ! -f "$MT_DNA" ]]; then
  echo "Error: mtDNA FASTA not found: $MT_DNA" >&2
  exit 1
fi
if [[ ! -d "$SPECIES_DIR" ]]; then
  echo "Error: species directory not found: $SPECIES_DIR" >&2
  exit 1
fi

# Determine CPUS if not set
if [[ -z "$CPUS" ]]; then
  CPUS=$(command -v nproc &>/dev/null && nproc || echo 1)
fi

# Prepare directories
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Loop over FASTAs (exclude ax4.fa by default)
shopt -s nullglob
for FASTA in "$SPECIES_DIR"/*.fa; do
  BASENAME=$(basename "$FASTA")
  [[ "$BASENAME" == "ax4.fa" ]] && {
    echo "Skipping excluded file: $BASENAME"
    continue
  }
  SPECIES=${BASENAME%.fa}
  LOG_FILE="$LOG_DIR/${SPECIES}_mtblastn.log"
  DB_PREFIX="$OUT_DIR/${SPECIES}_db"
  OUT_FILE="$OUT_DIR/${SPECIES}_vs_mtDNA.blastn"

  echo " === $SPECIES ===" | tee "$LOG_FILE"

  echo "Building BLAST DB for $FASTA..." | tee -a "$LOG_FILE"
  makeblastdb -in "$FASTA" -dbtype nucl -out "$DB_PREFIX" \
    >> "$LOG_FILE" 2>&1

  echo "Running blastn of $MT_DNA vs. $SPECIES contigs..." | tee -a "$LOG_FILE"
  blastn \
    -query "$MT_DNA" \
    -db "$DB_PREFIX" \
    -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend sstart send' \
    -evalue 1e-10 \
    -num_threads "$CPUS" \
    -out "$OUT_FILE" \
    >> "$LOG_FILE" 2>&1

  if [[ $? -ne 0 ]]; then
    echo "ERROR: blastn failed for $SPECIES" | tee -a "$LOG_FILE" >&2
  else
    echo "Completed: results in $OUT_FILE" | tee -a "$LOG_FILE"
  fi
done

echo "All done! Logs => $LOG_DIR/    Results => $OUT_DIR/"
