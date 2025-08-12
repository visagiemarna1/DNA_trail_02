#!/usr/bin/env bash
set -euo pipefail

ROOT="$HOME/documents/DNA_trail_02"
RAW="$ROOT/Fastq"
OUT="$ROOT/output"
TMP="$OUT/tmp"
MERGED="$TMP/merged"
mkdir -p "$OUT" "$TMP" "$MERGED"

# ---------- tweak if needed ----------
STRIPLEFT_R1=0
STRIPLEFT_R2=0
TRUNCLEN_R1=0
TRUNCLEN_R2=0
MIN_OVERLAP=20
MAXEE=1.0
TRUNCLEN_MERGED=0
# -------------------------------------

: > "$TMP/merged_list.txt"
: > "$OUT/merge_summary.tsv"   # sampleID\tmerged_fastq_size_bytes

shopt -s nullglob
for R1 in "$RAW"/*_R1_001.fastq "$RAW"/*_R1_001.fastq.gz; do
  [[ -e "$R1" ]] || continue
  base=$(basename "$R1")
  sample_full=${base%%_L001_*}         # e.g., 293_S71
  sample_id=${sample_full%%_*}         # e.g., 293

  # find matching R2 preserving extension
  if [[ "$R1" == *.gz ]]; then
    R2="$RAW/${base/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  else
    R2="$RAW/${base/_R1_001.fastq/_R2_001.fastq}"
  fi
  if [[ ! -f "$R2" ]]; then
    echo "ERROR: Missing R2 for sam    echo "ERROR: Missing R2    exit 1
  fi

  # Trim
  tR1="$TMP/${sample_id}_R1.trim.fastq"
  tR2="$TMP/${sample_id}_R2.trim.fastq"
  args1=(-fastq_filter "$R1" -fastqout "$tR1"); ((STRIPLEFT_R1>0))&&args1+=(-fastq_stripleft $STRIPLEFT_R1); ((TRUNCLEN_R1>0))&&args1+=(-fastq_trunclen $TRUNCLEN_R1)
  args2=(-fastq_filter "$R2" -fastqout "$tR2"); ((STRIPLEFT_R2>0))&&args2+=(-fastq_stripleft $STRIPLEFT_R2); ((TRUNCLEN_R2>0))&&args2+=(-fastq_trunclen $TRUNCLEN_R2)
  usearch "${args1[@]}"
  usearch "${args2[@]}"

  # Merge (fail if empty)
  report="$MERGED/${sample_id}_merge_report.txt"
  m="$MERGED/${sample_id}_merged.fastq"
  usearch -fastq_mergepairs "$tR1" \
          -reverse "$tR2" \
          -fastq_minovlen $MIN_OVERLAP \
          -fastqout "$m" \
          -report "$report"

  if [[ ! -s "$m" ]]; then
    echo "ERROR: Merge produced no reads for sample ${sample_id}." >&2
    echo "See report: $report" >&2
    exit 1
  fi

  echo -e "${sample_id}\t$(stat -f%z "$m" 2>/dev/null || stat -c%s "$m")" >> "$OUT/merge_summary.tsv"
  echo "$m" >> "$TMP/merged_list.txt"
done

if [[ ! -s "$TMP/merged_list.txt" ]]; then
  echo "ERROR: No merged FASTQs produced (unexpected in strict mode)." >&2
  exit 1
fi

# Quality filter all merged and combine
: > "$TMP/all_filtered.fasta"
while IFS= read -r mfq; do
  [[ -n "$mfq" ]] || continue
  in="$mfq"
  if (( TRUNCLEN_MERGED > 0 )); then
    trm="${mfq%.fastq}.trunc.fastq"
    usearch -fastq_filter "$mfq" -fastq_trunclen $TRUNCLEN_MERGED -fastqout "$trm"
    in="$trm"
  fi
  ffa="${mfq%.fastq}_filtered.fasta"
  usearch -fastq_filter "$in" -fastq_maxee $MAXEE -fastaout "$ffa"
  cat "$ffa" >> "$TMP/all_filtered.fasta"
done < "$TMP/merged_list.txt"

# UPARSE clustering (removes singletons + chimeras)
usearch -fastx_uniques "$TMP/all_filtered.fasta" -sizeout -relabel Uniq -fastaout "$TMP/uniques.fasta"
usearch -cluster_otus "$TMP/uniques.fasta" -minsize 2 -otus "$OUT/OTU_sequences_for_BLAST.fasta"

# Build list of merged files sorted numerically by sample ID
MERGED_FILES=$(while read -r f; do
  [[ -s "$f" ]] || { echo "ERROR: Empty merged file encountered unexpectedly: $f" >&2; exit 1; }
  id=$(basename "$f" | cut -d'_' -f1)
  printf "%s\t%s\n" "$id" "$f"
done < "$TMP/merged_list.txt" | sort -n -k1,1 | cut -f2-)

# OTU table (numeric IDs, sorted columns)
usearch -otutab $MERGED_FILES \
        -otus "$OUT/OTU_sequences_for_BLAST.fasta" \
        -otutabout "$OUT/otutable.txt" \
        -sample_delim "_"

# Cleanup intermediates
rm -rf "$TMP"

echo "Strict run complete.
Main outputs:
  - $OUT/OTU_sequences_for_BLAST.fasta
  - $OUT/otutable.txt
Quick check:
  - $OUT/merge_summary.tsv (size of each merged FASTQ)"
