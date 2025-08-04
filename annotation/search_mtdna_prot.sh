#!/usr/bin/env bash
#
# search_mtdna_prot.sh — Build BLAST DBs and run tBLASTn of mitochondrial proteins vs. contigs
#
# Usage:
#   search_mtdna_prot.sh -m MT_PROT -s SPECIES_DIR -o OUT_DIR [-l LOG_DIR] [-c CPUS]
#
# Options:
#   -m MT_PROT       FASTA of mitochondrial proteins (e.g. mt_proteins.faa)
#   -s SPECIES_DIR   Directory containing contig FASTAs (*.fa)
#   -o OUT_DIR       Directory to write BLAST DBs and tBLASTn results
#   -l LOG_DIR       Directory for log files (default: ./logs)
#   -c CPUS          Number of threads for tblastn (default: detected cores)
#   -h               Show this help message and exit
#

usage(){
  cat <<EOF
Usage: $(basename "$0") -m MT_PROT -s SPECIES_DIR -o OUT_DIR [-l LOG_DIR] [-c CPUS]

Options:
  -m MT_PROT       FASTA of mitochondrial proteins (required)
  -s SPECIES_DIR   Directory containing contig FASTAs (*.fa) (required)
  -o OUT_DIR       Directory to write BLAST DBs and results (required)
  -l LOG_DIR       Directory for logs (default: ./logs)
  -c CPUS          Threads for tblastn (default: detected cores)
  -h               Show this help and exit
EOF
}

#load required modules
module load gcc/8.2.0
module load perl/5.28.0
module load blast+

# defaults
LOG_DIR="./logs"
CPUS=""
EXCLUDE_BASENAME="ax4.fa"

# parse options
while getopts ":m:s:o:l:c:h" opt; do
  case $opt in
    m) MT_PROT=$OPTARG ;;
    s) SPECIES_DIR=$OPTARG ;;
    o) OUT_DIR=$OPTARG ;;
    l) LOG_DIR=$OPTARG ;;
    c) CPUS=$OPTARG ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done

# check required args
if [[ -z "${MT_PROT:-}" || -z "${SPECIES_DIR:-}" || -z "${OUT_DIR:-}" ]]; then
  echo "Error: -m, -s and -o are required." >&2
  usage
  exit 1
fi

# detect CPUS if not provided
if [[ -z "$CPUS" ]]; then
  if command -v nproc &>/dev/null; then
    CPUS=$(nproc)
  else
    CPUS=1
  fi
fi

# validations
if [[ ! -f "$MT_PROT" ]]; then
  echo "Error: mt proteins FASTA not found: $MT_PROT" >&2
  exit 1
fi
if [[ ! -d "$SPECIES_DIR" ]]; then
  echo "Error: species directory not found: $SPECIES_DIR" >&2
  exit 1
fi

# prepare directories
mkdir -p "$OUT_DIR" "$LOG_DIR"

# loop over contig FASTAs
shopt -s nullglob
for FASTA in "$SPECIES_DIR"/*.fa; do
  BASENAME=$(basename "$FASTA")
  [[ "$BASENAME" == "$EXCLUDE_BASENAME" ]] && {
    echo "Skipping excluded file: $BASENAME"
    continue
  }
  SPECIES=${BASENAME%.fa}
  LOG_FILE="$LOG_DIR/${SPECIES}_mtblast.log"
  DB_PREFIX="$OUT_DIR/${SPECIES}_contigs_db"
  OUT_TBL="$OUT_DIR/${SPECIES}_vs_mtprot.tblastn"

  echo " === Processing $SPECIES ===" | tee "$LOG_FILE"

  echo "Building BLAST DB for $FASTA..." | tee -a "$LOG_FILE"
  makeblastdb -in "$FASTA" -dbtype nucl -out "$DB_PREFIX" \
    >> "$LOG_FILE" 2>&1

  echo "Running tBLASTn of $MT_PROT vs. $SPECIES contigs..." | tee -a "$LOG_FILE"
  tblastn \
    -query "$MT_PROT" \
    -db "$DB_PREFIX" \
    -evalue 1e-5 \
    -outfmt 6 \
    -num_threads "$CPUS" \
    -out "$OUT_TBL" \
    >> "$LOG_FILE" 2>&1

  if [[ $? -ne 0 ]]; then
    echo "ERROR: tBLASTn failed for $SPECIES" | tee -a "$LOG_FILE" >&2
  else
    echo "Done with $SPECIES: results in $OUT_TBL" | tee -a "$LOG_FILE"
  fi
done

echo "All species processed. Logs in $LOG_DIR, results in $OUT_DIR."
