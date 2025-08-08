#!/usr/bin/env bash
set -euo pipefail
ROOT="$HOME/documents/DNA_trail_02"
RAW="$ROOT/Fastq"
OUT="$ROOT/output"

mkdir -p "$OUT"

# Example for one pair; adapt or loop as needed
R1="$RAW/sample_R1.fastq"
R2="$RAW/sample_R2.fastq"

usearch -fastq_mergepairs "$R1" -reverse "$R2" -fastqout "$OUT/merged.fastq"
usearch -fastq_filter "$OUT/merged.fastq" -fastq_maxee 1.0 -fastaout "$OUT/filtered.fasta"
usearch -fastx_uniques "$OUT/filtered.fasta" -fastaout "$OUT/uniques.fasta" -sizeout
usearch -unoise3 "$OUT/uniques.fasta" -zotus "$OUT/zotus.fasta"
usearch -otutab "$OUT/merged.fastq" -zotus "$OUT/zotus.fasta" -otutabout "$OUT/otutable.txt"

echo "Done. Results in $OUT"
