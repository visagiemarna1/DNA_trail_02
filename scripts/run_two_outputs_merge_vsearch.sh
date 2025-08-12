#!/usr/bin/env bash
set -euo pipefail

ROOT="$HOME/documents/DNA_trail_02"
RAW="$ROOT/Fastq"
OUT="$ROOT/output"
TMP="$OUT/tmp"
MERGED="$TMP/merged"
mkdir -p "$OUT" "$TMP" "$MERGED"

# ---- params ----
MIN_OVERLAP=20       # for merging (vsearch)
MAXEE=1.0            # expected errors on merged reads (usearch)
TRUNCLEN_MERGED=0    # 0 = no truncation, else e.g. 220
# -----------------

# sanity: require vsearch + usearch
command -v vsearch >/dev/null || { echo "ERROR: vsearch not in PATH"; exit 1; }
command -v usearch >/dev/null || { echo "ERROR: usearch not in PATH"; exit 1; }

: > "$TMP/merged_list.txt"

shopt -s nullglob
for R1 in "$RAW"/*_R1_001.fastq "$RAW"/*_R1_001.fastq.gz; do
  [[ -e "$R1" ]] || continue
  base=$(basename "$R1")
  sample_full=${base%%_L001_*}      # e.g., 293_S71
  sample_id=${sample_full%%_*}      # e.g., 293

  # match R2 by extension
  if [[ "$R1" == *.gz ]]; then
    R2="$RAW/${base/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  else
    R2="$RAW/${base/_R1_001.fastq/_R2_001.fastq}"
  fi
  [[ -f "$R2" ]] || { echo "ERROR: Missing R2 for $sample_id"; exit 1; }

  m="$MERGED/${sample_id}_merged.fastq"
  # --- MERGE with VSEARCH ---
  vsearch --fastq_mergepairs "$R1" \
          --reverse "$R2" \
          --fastq_minovlen $MIN_OVERLAP \
          --fastqout "$m" \
          --fastq_eeout

  [[ -s "$m" ]] || { echo "ERROR: Merge produced no reads for $sample_id"; exit 1; }
  echo "$m" >> "$TMP/merged_list.txt"
done

# QUALITY FILTER merged reads (USEARCH) and combine to one FASTA
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

# UPARSE clustering (USEARCH) -> OTUs (singleton + chimera filtered)
usearch -fastx_uniques "$TMP/all_filtered.fasta" -sizeout -relabel Uniq -fastaout "$TMP/uniques.fasta"
usearch -cluster_otus "$TMP/uniques.fasta" -minsize 2 -otus "$OUT/OTU_sequences_for_BLAST.fasta"

# Build OTU table (samples = numeric IDs, sorted left→right)
MERGED_FILES=$(while read -r f; do
  [[ -s "$f" ]] || { echo "ERROR: Empty merged file: $f"; exit 1; }
  id=$(basename "$f" | cut -d'_' -f1)
  printf "%s\t%s\n" "$id" "$f"
done < "$TMP/merged_list.txt" | sort -n -k1,1 | cut -f2-)

usearch -otutab $MERGED_FILES \
        -otus "$OUT/OTU_sequences_for_BLAST.fasta" \
        -otutabout "$OUT/otutable.txt" \
        -sample_delim "_"

# Cleanup
rm -rf "$TMP"

echo "Done.
Main outputs:
  - $OUT/OTU_sequences_for_BLAST.fasta
  - $OUT/otutable.txt (numeric IDs, sorted)"
