#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# run_craq_single.sh
#
# Usage:
#   run_craq_single.sh -n <nanopore.{fastq.gz|bam}> \
#                     -a <assembly.fa> \
#                     -i <illumina.{fastq.gz|bam}[,<illumina2.{fastq.gz|bam}>] \
#                     [-t threads] \
#                     [-o outdir]
# -----------------------------------------------------------------------------

usage() {
  echo "Usage: $0 -n <nanopore.{fastq.gz|bam}> -a <assembly.fa> -i <illumina.{fastq.gz|bam}[,<illumina2.{fastq.gz|bam}>] [-t threads] [-o outdir]"
  exit 1
}

# Defaults
THREADS="$(nproc)"
OUTDIR="craq_output"

# Parse args
while getopts "n:a:i:t:o:h" opt; do
  case "$opt" in
    n) NANO_FILE="$OPTARG" ;;
    a) ASSEMBLY="$OPTARG" ;;
    i) ILLUMINA_INPUT="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    *) usage ;;
  esac
done

# Check required
if [ -z "${NANO_FILE:-}" ] || [ -z "${ASSEMBLY:-}" ] || [ -z "${ILLUMINA_INPUT:-}" ]; then
  usage
fi

# Derive species name
SPECIES="$(basename "$ASSEMBLY" .fa)"
SPECIES_LOWER="$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')"

echo "→ Species:      $SPECIES_LOWER"
echo "→ Assembly:     $ASSEMBLY"
echo "→ Threads:      $THREADS"
echo

# -----------------------------------------------------------------------------
# 1) Prepare Nanopore BAM
# -----------------------------------------------------------------------------
case "$NANO_FILE" in
  *.fastq|*.fastq.gz)
    echo "Mapping Nanopore reads to BAM..."
    NANO_BAM="nanopore.sorted.bam"
    minimap2 -a "$ASSEMBLY" "$NANO_FILE" \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$NANO_BAM"
    samtools index "$NANO_BAM"
    ;;
  *.bam)
    echo "Using provided Nanopore BAM..."
    NANO_BAM="$NANO_FILE"
    [ -f "${NANO_BAM}.bai" ] || samtools index "$NANO_BAM"
    ;;
  *)
    echo "Error: Nanopore input must be .fastq, .fastq.gz, or .bam" >&2
    exit 1
    ;;
esac

# -----------------------------------------------------------------------------
# 2) Prepare Illumina BAM (supports 1 or 2 inputs via comma)
# -----------------------------------------------------------------------------
# split on comma
IFS=',' read -r -a ILL_FILES <<< "$ILLUMINA_INPUT"
NUM_ILL=${#ILL_FILES[@]}

if (( NUM_ILL == 1 )); then
  I1="${ILL_FILES[0]}"
  ext="${I1##*.}"
  case "$ext" in
    bam)
      echo "Using single Illumina BAM..."
      ILLUMINA_BAM="$I1"
      [ -f "${ILLUMINA_BAM}.bai" ] || samtools index "$ILLUMINA_BAM"
      ;;
    fastq|gz)
      echo "Mapping single Illumina FASTQ..."
      ILLUMINA_BAM="illumina.sorted.bam"
      minimap2 -a -x sr "$ASSEMBLY" "$I1" \
        | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$ILLUMINA_BAM"
      samtools index "$ILLUMINA_BAM"
      ;;
    *)
      echo "Error: Illumina input must be .fastq, .fastq.gz, or .bam" >&2
      exit 1
      ;;
  esac

elif (( NUM_ILL == 2 )); then
  I1="${ILL_FILES[0]}"
  I2="${ILL_FILES[1]}"
  ext1="${I1##*.}"
  ext2="${I2##*.}"

  if [[ "$ext1" == "bam" && "$ext2" == "bam" ]]; then
    echo "Merging two Illumina BAMs..."
    ILLUMINA_BAM="illumina.merged.bam"
    samtools merge -@ "$THREADS" "$ILLUMINA_BAM" "$I1" "$I2"
    samtools index "$ILLUMINA_BAM"

  elif ([[ "$ext1" == "fastq" || "$ext1" == "gz" ]] && [[ "$ext2" == "fastq" || "$ext2" == "gz" ]]); then
    echo "Mapping paired Illumina FASTQs..."
    ILLUMINA_BAM="illumina.sorted.bam"
    minimap2 -a -x sr "$ASSEMBLY" "$I1" "$I2" \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$ILLUMINA_BAM"
    samtools index "$ILLUMINA_BAM"

  else
    echo "Error: For two inputs, both must be FASTQs or both must be BAMs." >&2
    exit 1
  fi

else
  echo "Error: Too many Illumina inputs (max 2 allowed)." >&2
  exit 1
fi

# -----------------------------------------------------------------------------
# 3) Run CRAQ
# -----------------------------------------------------------------------------
mkdir -p "$OUTDIR"
echo "Running CRAQ…"
craq \
  -g "$ASSEMBLY" \
  -sms "$NANO_BAM" \
  -ngs "$ILLUMINA_BAM" \
  -D "$OUTDIR"

echo
echo "CRAQ complete. Results in: $OUTDIR"
