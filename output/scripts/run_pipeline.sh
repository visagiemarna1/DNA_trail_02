#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# DNA Trail Pipeline (merged+final only in MultiQC)
# - Merge (USEARCH) -> Trim (cutadapt) -> QC filter (USEARCH, maxEE=1.0)
# - Stats CSV: counts + meanQ at: before-merge, after-merge(raw), after-cleaning
# - FastQC runs ONLY on merged (raw) and final cleaned files
# - MultiQC shows plots for merged & final only; CSV is embedded as a table
# - Outputs organized into subfolders
# ==========================================================

# ---------- CONFIG ----------
FWD_PRIMER="TAGAACAGGCTCCTCTAG"
REV_PRIMER="TTAGATACCCCACTATGC"
REV_PRIMER_RC="GCATAGTGGGGTATCTAA"   # RC of reverse primer for 3' trimming on merged reads

# numeric formatting: force decimal DOT
export LC_ALL=C
export LC_NUMERIC=C

# Dry-run support
DRY_RUN="${1:-false}"
[[ "$DRY_RUN" == "--dry-run" ]] && DRY_RUN="true"

# ---------- PATHS ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_DIR="${REPO_ROOT}/Fastq"

OUTPUT_DIR="${REPO_ROOT}/output"
MERGED_DIR="${OUTPUT_DIR}/merged"
TRIMMED_DIR="${OUTPUT_DIR}/trimmed"
CLEAN_DIR="${OUTPUT_DIR}/clean"
FASTQC_DIR="${OUTPUT_DIR}/fastqc"         # we will clear and then write ONLY merged+final FastQC here
MQC_DIR="${OUTPUT_DIR}/multiqc"
MQC_CUSTOM_DIR="${OUTPUT_DIR}/multiqc_custom"
LOGS_DIR="${REPO_ROOT}/logs"

STATS_FILE="${OUTPUT_DIR}/merge_stats.csv"
MQC_TABLE="${MQC_CUSTOM_DIR}/merge_stats_mqc.csv"

mkdir -p "$MERGED_DIR" "$TRIMMED_DIR" "$CLEAN_DIR" "$FASTQC_DIR" "$MQC_DIR" "$MQC_CUSTOM_DIR" "$LOGS_DIR"

# ---------- TOOLS ----------
USEARCH_BIN="$(command -v usearch || true)"
CUTADAPT_BIN="$(command -v cutadapt || true)"
FASTQC_BIN="$(command -v fastqc || true)"   # optional
MULTIQC_BIN="$(command -v multiqc || true)"

[[ -z "$USEARCH_BIN"  ]] && { echo "ERROR: usearch not found"; exit 1; }
[[ -z "$CUTADAPT_BIN" ]] && { echo "ERROR: cutadapt not found"; exit 1; }
[[ -z "$MULTIQC_BIN"  ]] && { echo "ERROR: multiqc not found"; exit 1; }

echo "Repo      : $REPO_ROOT"
echo "Input     : $INPUT_DIR"
echo "Outputs   : $OUTPUT_DIR"
echo "Tools     :"
echo "  usearch : $USEARCH_BIN"
echo "  cutadapt: $CUTADAPT_BIN"
echo "  fastqc  : ${FASTQC_BIN:-SKIP}"
echo "  multiqc : $MULTIQC_BIN"
echo "Dry-run   : $DRY_RUN"
echo

# ---------- HELPERS ----------
run_or_echo(){ if [[ "$DRY_RUN" == "true" ]]; then echo "[DRY-RUN] $*"; else eval "$@"; fi; }

# Count reads in FASTQ (lines/4)
count_reads(){ local f="$1"; echo $(( $(wc -l < "$f") / 4 )); }

# Mean Phred Q across all bases (Phred+33)
mean_q(){
  local f="$1"
  perl -ne '$.%4==0 or next; chomp; for (split //){$s+=ord($_)-33;$n++} END{printf "%.2f\n",$n?$s/$n:0}' "$f"
}

# ---------- CSV HEADER (no GC%) ----------
echo "sample,reads_before_merge,meanQ_before,reads_after_merge_raw,meanQ_after_merge_raw,reads_after_cleaning,meanQ_after_cleaning" > "$STATS_FILE"

# ---------- ensure FastQC plots include ONLY merged & final ----------
# (Clear previous FastQC outputs so MultiQC won't pick up extra files)
if [[ "$DRY_RUN" != "true" ]]; then
  rm -rf "${FASTQC_DIR:?}/"* || true
fi

# ---------- PROCESS ----------
shopt -s nullglob
R1_FILES=( "${INPUT_DIR}"/*_R1_*.fastq )
shopt -u nullglob
(( ${#R1_FILES[@]} > 0 )) || { echo "No *_R1_*.fastq in ${INPUT_DIR}"; exit 1; }

for r1 in "${R1_FILES[@]}"; do
  r2="${r1/_R1_/_R2_}"
  sample="$(basename "$r1" | sed 's/_R1_.*//')"

  merged="${MERGED_DIR}/${sample}_merged.fastq"             # included in MultiQC (FastQC)
  trimmed="${TRIMMED_DIR}/${sample}_trimmed.fastq"          # NOT included in FastQC plots
  clean="${CLEAN_DIR}/${sample}_clean_FINAL.fastq"          # included in MultiQC (FastQC)

  if [[ ! -f "$r2" ]]; then
    echo "WARNING: missing R2 for $(basename "$r1") — skipping $sample."
    continue
  fi

  # BEFORE MERGE (use R1 as proxy for pairs)
  rb=$(count_reads "$r1")
  qb=$(mean_q "$r1")

  echo "Merging $sample -> $(basename "$merged")"
  run_or_echo "$USEARCH_BIN -fastq_mergepairs '$r1' -reverse '$r2' -fastqout '$merged' 2>&1 | tee '${LOGS_DIR}/${sample}_usearch_merge.log' >/dev/null"

  if [[ "$DRY_RUN" == "true" ]]; then
    echo "${sample},${rb},${qb},NA,NA,NA,NA" >> "$STATS_FILE"
    continue
  fi

  if [[ ! -s "$merged" ]]; then
    echo "ERROR: merged file empty for $sample"
    echo "${sample},${rb},${qb},0,0,0,0" >> "$STATS_FILE"
    continue
  fi

  # AFTER MERGE (raw, untrimmed)
  ram=$(count_reads "$merged")
  qam=$(mean_q "$merged")

  # Run FastQC on MERGED (raw) so MultiQC shows this set
  if [[ -n "${FASTQC_BIN}" ]]; then
    run_or_echo "$FASTQC_BIN -q -o '${FASTQC_DIR}' '$merged'"
  fi

  echo "Trim primers $sample -> $(basename "$trimmed")"
  run_or_echo "$CUTADAPT_BIN -e 0.1 -m 50 -g ^${FWD_PRIMER} -a ${REV_PRIMER_RC} -o '$trimmed' '$merged' > '${LOGS_DIR}/${sample}_cutadapt.txt' 2>&1"

  echo "QC filter (maxEE=1.0) $sample -> $(basename "$clean")"
  run_or_echo "$USEARCH_BIN -fastq_filter '$trimmed' -fastq_maxee 1.0 -fastqout '$clean' 2>&1 | tee '${LOGS_DIR}/${sample}_usearch_filter.log' >/dev/null"

  # AFTER CLEANING (final)
  rac=0; qac=0
  if [[ -s "$clean" ]]; then
    rac=$(count_reads "$clean")
    qac=$(mean_q "$clean")
  fi

  # Run FastQC on FINAL only (so plots show merged + final; no R1/R2 or trimmed)
  if [[ -n "${FASTQC_BIN}" && -s "$clean" ]]; then
    run_or_echo "$FASTQC_BIN -q -o '${FASTQC_DIR}' '$clean'"
  fi

  echo "${sample},${rb},${qb},${ram},${qam},${rac},${qac}" >> "$STATS_FILE"
done

# ---------- MEAN ROW ----------
if [[ "$DRY_RUN" != "true" ]]; then
  awk -F',' '
    NR==1{next} $1=="MEAN"{next}
    {n++; s2+=$2; s3+=$3; s4+=$4; s5+=$5; s6+=$6; s7+=$7}
    END{
      if(n>0) printf "MEAN,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
        s2/n,s3/n,s4/n,s5/n,s6/n,s7/n
    }' "$STATS_FILE" >> "$STATS_FILE"
fi

# ---------- MULTIQC CUSTOM TABLE (embed CSV without GC%) ----------
if [[ "$DRY_RUN" != "true" ]]; then
  {
    echo "# plot_type: 'table'"
    echo "# section_name: 'Merge & final QC summary (pipeline)'"
    echo "# description: 'Reads and mean Q per sample at: before merge, after merge (raw), after cleaning (final)'"
    echo "Sample,Reads before merge,MeanQ before,Reads after merge (raw),MeanQ after merge (raw),Reads after cleaning,MeanQ after cleaning"
    tail -n +2 "$STATS_FILE"  # includes MEAN row
  } > "$MQC_TABLE"
fi

# ---------- RUN MULTIQC ----------
if [[ "$DRY_RUN" != "true" ]]; then
  echo "Running MultiQC…"
  # Only fastqc (merged + final) + logs + custom table are scanned
  run_or_echo "$MULTIQC_BIN '$LOGS_DIR' '$FASTQC_DIR' '$MQC_CUSTOM_DIR' --outdir '$MQC_DIR' --force"
  echo "MultiQC report: ${MQC_DIR}/multiqc_report.html"
fi

echo "Done. Stats CSV: $STATS_FILE"

