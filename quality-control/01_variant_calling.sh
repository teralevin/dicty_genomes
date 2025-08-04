#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# 01_variant_calling.sh
#
# Usage:
#   # with defaults (PL=1, CUTOFF=0.05)
#   ./01_variant_calling.sh <assembly.fasta> <reads_R1.fastq> [reads_R2.fastq]
#
#   # override PL & CUTOFF
#   ./01_variant_calling.sh <ploidy> <cutoff> <assembly.fasta> <reads_R1.fastq> [reads_R2.fastq]
#
# Examples:
#   ./01_variant_calling.sh ref.fa R1.fastq R2.fastq
#   ./01_variant_calling.sh 2 0.10 ref.fa R1.fq
# -----------------------------------------------------------------------------

usage() {
  echo "Usage:"
  echo "  $0 [<ploidy> <cutoff>] <assembly.fasta> <reads_R1.fastq> [reads_R2.fastq]"
  exit 1
}

# ---- parse PL & CUTOFF if given (must have at least 4 args to override) ----
if (( $# >= 4 )); then
  PL="$1"           # ploidy
  CUTOFF="$2"       # allele‐freq cutoff
  REF="$3"          # assembly
  R1="$4"           # read1
  R2="${5-}"        # read2 (optional)
elif (( $# >= 2 )); then
  PL=1
  CUTOFF=0.05
  REF="$1"
  R1="$2"
  R2="${3-}"
else
  usage
fi

# ---- sanity check ----
if [[ -z "$REF" || -z "$R1" ]]; then
  usage
fi

# derive basename
BASE="$(basename "$REF")"
BASE="${BASE%.*}"

echo "Starting variant calling"
echo "   Ploidy:          $PL"
echo "   AF cutoff:       $CUTOFF"
echo "   Assembly:        $REF"
echo "   Reads R1:        $R1"
if [[ -n "$R2" ]]; then echo "   Reads R2:        $R2"; fi
echo "   Output prefix:   $BASE"
echo

# ---- load modules ----
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.9
module load bcftools/1.20

# ---- ensure BWA index ----
if [[ ! -f "${REF}.bwt" ]]; then
  echo "Building BWA index for $REF"
  bwa index "$REF"
fi

# ---- alignment & sorting ----
echo "Aligning reads with BWA MEM"
if [[ -n "$R2" ]]; then
  bwa mem "$REF" "$R1" "$R2" \
    | samtools view -b - \
    | samtools sort -o "${BASE}.bam"
else
  bwa mem "$REF" "$R1" \
    | samtools view -b - \
    | samtools sort -o "${BASE}.bam"
fi
samtools index "${BASE}.bam"

# ---- variant calling ----
echo "Calling variants (bcftools mpileup + call)"
bcftools mpileup -Ou -f "${REF}" "${BASE}.bam" \
  | bcftools call -mv -Ob --ploidy "${PL}" -o "${BASE}.bcf"

# ---- allele‐frequency extraction ----
echo "Extracting allele frequencies"
bcftools query -f '%CHROM\t%POS\t%AN\t%AC{0}\n' "${BASE}.bcf" \
  | awk -v CUT="$CUTOFF" '{
      af = $4/$3;
      printf("%s\t%s\t%d\t%d\t%.6f\n",$1,$2,$3,$4,af)
    }' > "${BASE}.variant_freq.txt"

# ---- report mixed‐allele sites ----
echo "Reporting mixed‐allele sites"
awk -v CUT="$CUTOFF" '$5 < 1.0 && $5 >= CUT' "${BASE}.variant_freq.txt" \
  > "${BASE}.mixed_variants.txt"

# ---- all‐sites pileup & call ----
echo "Calling ALL sites for error‐rate estimation"
bcftools mpileup --fasta-ref "${REF}" --annotate FORMAT/DP,AD "${BASE}.bam" \
  | bcftools call --multiallelic-caller --output-type b -o "${BASE}.all_sites.bcf"
bcftools index "${BASE}.all_sites.bcf"

# ---- extract all‐sites allele frequencies ----
echo "Extracting all‐sites AF"
bcftools query -f '%CHROM\t%POS\t[%DP]\t[%AD]\n' "${BASE}.all_sites.bcf" \
  | awk '{
      split($4,ad,",");
      tot = ad[1]+ad[2];
      if (tot>0) {
        af = ad[2]/tot;
        print $1, $2, tot, ad[1], ad[2], af;
      }
    }' > "${BASE}.all_sites_af.txt"

# ---- summarize error vs. mixed‐allele sites ----
echo "Summarizing error rates (cutoff = $CUTOFF)"
F="${BASE}.all_sites_af.txt"
TOTAL=$(wc -l < "$F")
ERR=$(awk -v t="$CUTOFF" '$6 < t' "$F" | wc -l)
VAR=$(awk -v t="$CUTOFF" '$6 >= t' "$F" | wc -l)
PCT_ERR=$(awk 'BEGIN{printf "%.2f", '"$ERR"'/'"$TOTAL"'*100}')
PCT_VAR=$(awk 'BEGIN{printf "%.2f", '"$VAR"'/'"$TOTAL"'*100}')

cat <<EOF > "${BASE}.error_summary.txt"
Total sites surveyed:        $TOTAL
Sites with AF < $CUTOFF:      $ERR  (${PCT_ERR}%)
Sites with AF ≥ $CUTOFF:      $VAR  (${PCT_VAR}%)
EOF

echo "Done. Summary in ${BASE}.error_summary.txt"
