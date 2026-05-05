#!/bin/bash
# ============================================================
# Script: run_extract_target_orfs_with_assignment.sh
# Purpose: Batch extraction of target ORFs and family assignment
# Stage: Annotation
# Input:
#   - /home/chenlu/skx_cd/diamond_result/*.tsv
#   - /home/chenlu/skx_cd/prodigal/*.genes.fna
# Output:
#   - /home/chenlu/skx_cd/target_orfs/*.target_orfs.fna
#   - /home/chenlu/skx_cd/target_orfs/*.confident_target_orfs.fna
#   - /home/chenlu/skx_cd/target_orfs/*.target_orfs.all_hits.tsv
#   - /home/chenlu/skx_cd/target_orfs/*.target_orfs.assignment.tsv
# ============================================================
set -euo pipefail

DIAMOND_DIR=/home/chenlu/skx_cd/diamond_result
PRODIGAL_DIR=/home/chenlu/skx_cd/prodigal
OUT_DIR=/home/chenlu/skx_cd/target_orfs
SCRIPT=/home/chenlu/skx_cd/scripts/extract_target_orfs_with_assignment.py

mkdir -p "$OUT_DIR"

for tsv in "${DIAMOND_DIR}"/*.tsv; do
    sample=$(basename "$tsv" .tsv)
    genes="${PRODIGAL_DIR}/${sample}.genes.fna"

    if [[ ! -f "$genes" ]]; then
        echo "[WARN] genes.fna not found for sample: $sample"
        continue
    fi

    echo "[INFO] processing: $sample"

    python "$SCRIPT" \
      -d "$tsv" \
      -g "$genes" \
      --out-all-fna "${OUT_DIR}/${sample}.target_orfs.fna" \
      --out-confident-fna "${OUT_DIR}/${sample}.confident_target_orfs.fna" \
      --out-all-hits "${OUT_DIR}/${sample}.target_orfs.all_hits.tsv" \
      --out-assignment "${OUT_DIR}/${sample}.target_orfs.assignment.tsv" \
      --max-evalue 1e-5 \
      --min-bitscore 50 \
      --min-aln-len 60 \
      --bitscore-gap 10
done
