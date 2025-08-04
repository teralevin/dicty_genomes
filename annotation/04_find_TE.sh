#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 04_find_TE.sh
#
# Usage:
#   04_find_TE.sh <contigs_dir> [threads]
#
#   <contigs_dir> : path to directory containing your .fa contig files
#   [threads]     : number of CPU threads to use (default: all available)
# -----------------------------------------------------------------------------

usage() {
  cat <<EOF >&2
Usage: $0 <contigs_dir> [threads]
  <contigs_dir> : Path to directory with .fa files
  [threads]     : Number of CPU threads (default: $(nproc))
EOF
  exit 1
}

# ─── parse args ────────────────────────────────────────────────────────────────
if (( $# < 1 || $# > 2 )); then
  usage
fi

CONTIG_DIR="$1"
THREADS="${2:-$(nproc)}"

if [[ ! -d "$CONTIG_DIR" ]]; then
  echo "Error: '$CONTIG_DIR' is not a directory." >&2
  exit 1
fi

echo "→ Contigs dir: $CONTIG_DIR"
echo "→ Threads:     $THREADS"
echo

# ─── load modules ─────────────────────────────────────────────────────────────
module load gcc/8.2.0
module load perl/5.28.0
module load blast+

# ─── download TE file ─────────────────────────────────────────────────────────
TE_FILE="TE_amoebozoa_500_ddis.fa"
if [[ ! -f "$TE_FILE" ]]; then
  echo "Downloading transposable elements to '$TE_FILE'…"
  wget -q -O "$TE_FILE" \
    https://raw.githubusercontent.com/Bart-Edelbroek/firmibasis_genome/main/transposable_elements/TE_amoebozoa_500_ddis.fa
else
  echo "Found existing '$TE_FILE', skipping download."
fi
echo

# ─── combine contigs ──────────────────────────────────────────────────────────
COMBINED_FASTA="all_species.fa"
echo "Combining all .fa from '$CONTIG_DIR' into '$COMBINED_FASTA'…"
cat "$CONTIG_DIR"/*.fa > "$COMBINED_FASTA"
echo

# ─── make BLAST DB ────────────────────────────────────────────────────────────
DB_NAME="dicty_db"
echo "Creating BLAST database '$DB_NAME' from '$COMBINED_FASTA'…"
makeblastdb -in "$COMBINED_FASTA" -dbtype nucl -out "$DB_NAME"
echo

# ─── run tblastx ──────────────────────────────────────────────────────────────
OUT_CSV="tblastx_hits_dicty.csv"
echo "Running tblastx (evalue=1e-10, sorthits=1)…"
tblastx \
  -db "$DB_NAME" \
  -query "$TE_FILE" \
  -evalue 1e-10 \
  -sorthits 1 \
  -num_threads "$THREADS" \
  -outfmt "10 delim=, sseqid bitscore evalue qlen qstart qend slen sstart send sseq length qseqid" \
  -out "$OUT_CSV"

echo
echo "Done. Results written to '$OUT_CSV'."
