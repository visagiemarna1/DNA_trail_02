
# DNA Trail Pipeline

This repository contains scripts for processing DNA sequencing data, specifically merging paired-end reads using [USEARCH](https://drive5.com/usearch/).

## Overview
The `run_pipeline.sh` script:
- Scans the `Fastq/` folder for Illumina-style paired files (`*_R1_*.fastq` and `*_R2_*.fastq`).
- Merges each pair using USEARCH.
- Saves merged FASTQ files into the `output/` folder.
- Can run in **dry-run mode** to preview actions before executing.

---

## Requirements
- [USEARCH](https://drive5.com/usearch/) installed and available in your `PATH`.
- Input files in `Fastq/` with names matching:

- Output directory (`output/`) will be created automatically.

---

## Running the Pipeline

### Preview (Dry-Run)
To see what the script will do without actually merging:
```bash
chmod +x scripts/run_pipeline.sh   # Only needed once
./scripts/run_pipeline.sh --dry-run

./scripts/run_pipeline.sh


output/

## Example Dry-Run output
./scripts/run_pipeline.sh --dry-run


### Example output
Repository root: /Users/username/Documents/DNA_trail_02
Input FASTQs   : /Users/username/Documents/DNA_trail_02/Fastq
Output folder  : /Users/username/Documents/DNA_trail_02/output
Dry-run mode   : true
Merging: 293_S71_L001_R1_001.fastq + 293_S71_L001_R2_001.fastq -> 293_S71_L001_merged.fastq
[DRY-RUN] usearch -fastq_mergepairs '/Users/username/Documents/DNA_trail_02/Fastq/293_S71_L001_R1_001.fastq' -reverse '/Users/username/Documents/DNA_trail_02/Fastq/293_S71_L001_R2_001.fastq' -fastqout '/Users/username/Documents/DNA_trail_02/output/293_S71_L001_merged.fastq'
...
Done. Merged: 0, Skipped: 0
Results: /Users/username/Documents/DNA_trail_02/output


### Real run output
./scripts/run_pipeline.sh


Repository root: /Users/username/Documents/DNA_trail_02
Input FASTQs   : /Users/username/Documents/DNA_trail_02/Fastq
Output folder  : /Users/username/Documents/DNA_trail_02/output
Dry-run mode   : false
Merging: 293_S71_L001_R1_001.fastq + 293_S71_L001_R2_001.fastq -> 293_S71_L001_merged.fastq
...
Done. Merged: 12, Skipped: 0
Results: /Users/username/Documents/DNA_trail_02/output
