#!/usr/bin/env bash
set -euo pipefail


# -------------------------------------------------------------------
# 02_extract_telomere_motif.sh
#  - takes a single TRF .dat file
#  - extracts repeats → BED
#  - writes consensus.fa
#  - runs blastn of telomeric motifs vs. consensus
#  - outputs BLAST hits and telomeric coords BED
# -------------------------------------------------------------------

# usage check
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <trf_output.dat>"
  exit 1
fi

#load the necessary modules
module load gcc/8.2.0
module load perl/5.28.0
module load blast+

dat_file="$1"
if [[ ! -f "$dat_file" ]]; then
  echo "Error: File '$dat_file' not found."
  exit 1
fi

# derive base name (strip extension)
fname=$(basename "$dat_file")
base="${fname%%.*}"

# create output directory
out_dir="./telomere_out_${base}"
mkdir -p "$out_dir"

# 1) Convert .dat → BED-like (seq name, 0-based start, end, consensus)
bed_file="${out_dir}/${base}_repeats.bed"
awk '
  /^Sequence:/      { seq=$2; next }
  /^[0-9]/          { printf "%s\t%d\t%d\t%s\n", seq, $1-1, $2, $(NF-1) }
' "$dat_file" > "$bed_file"

# 2) Write consensus FASTA
consensus_fa="${out_dir}/${base}_consensus.fa"
awk 'BEGIN{OFS="\n"} {
      chrom=$1; start=$2; end=$3; seq=$4;
      printf(">%s:%d-%d\n%s\n", chrom, start, end, seq)
    }' "$bed_file" > "$consensus_fa"

# 3) Ensure telomeric motifs FASTA exists
motif_fa="telomeric_motifs.fa"
if [[ ! -f "$motif_fa" ]]; then
  cat << 'EOF' > "$motif_fa"
>motif1
GAGGAGAGAGTCCCTTTTTTT
>motif2
GGGGAGAGACAGGGGAGAGACA
EOF
  echo "Created $motif_fa"
fi

# 4) Run blastn-short of motifs vs. consensus
out_tbl="${out_dir}/${base}_telo_matches.tsv"
blastn \
  -task blastn-short \
  -query "$motif_fa" \
  -subject "$consensus_fa" \
  -penalty -1 \
  -perc_identity 70 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  > "$out_tbl"

# 5) Extract matching sequence IDs
match_ids="${out_dir}/${base}_matching_ids.txt"
awk '{ print $2 }' "$out_tbl" | sort -u > "$match_ids"

# 6) Convert matching IDs to BED of coords (chrom, start, end)
coords_bed="${out_dir}/${base}_telomeric_coords.bed"
awk -F'[:-]' '{ print $1 "\t" $2 "\t" $3 }' "$match_ids" > "$coords_bed"

# Done
echo "All done for ${base}:"
echo "  BED of repeats:      $bed_file"
echo "  Consensus FASTA:     $consensus_fa"
echo "  BLAST hits table:    $out_tbl"
echo "  Matching IDs list:   $match_ids"
echo "  Telomeric coords BED:$coords_bed"