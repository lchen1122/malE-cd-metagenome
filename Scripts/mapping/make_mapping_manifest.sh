#!/bin/bash
# ============================================================
# Script: make_mapping_manifest.sh
# Purpose: Build mapping manifest linking each sample to its
#          confident target ORF reference and raw FASTQ files
# Stage: Mapping preparation
# ============================================================
set -euo pipefail

OUT=/home/chenlu/skx_cd/mapping_manifest.tsv
: > "$OUT"

for ref in /home/chenlu/skx_cd/target_orfs/*.confident_target_orfs.fna; do
    sample=$(basename "$ref" .confident_target_orfs.fna)

    [[ -s "$ref" ]] || continue

    raw=${sample#control_}
    raw=${raw#pollution_}

    if [[ "$sample" == control_* ]]; then
        group=control
        r1=/home/chenlu/skx_cd/control/data/${raw}_1.fastq
        r2=/home/chenlu/skx_cd/control/data/${raw}_2.fastq
    else
        group=pollution
        r1=/home/chenlu/skx_cd/cd_pollution/data/${raw}_1.fastq
        r2=/home/chenlu/skx_cd/cd_pollution/data/${raw}_2.fastq
    fi

    [[ -f "$r1" && -f "$r2" ]] || {
        echo "[WARN] missing raw fastq for $sample"
        continue
    }

    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$group" "$ref" "$r1" "$r2" >> "$OUT"
done

echo "Mapping manifest written:"
wc -l "$OUT"
head "$OUT"

