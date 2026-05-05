#!/bin/bash
# ============================================================
# Script: map_back_array.sh
# Purpose: Map raw reads back to sample-specific confident target
#          ORF references and calculate coverage
# Stage: Mapping / Quantification
# ============================================================
#SBATCH --job-name=map_back
#SBATCH --partition=Par-2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --error=/home/chenlu/skx_cd/log/%A_%a.err
#SBATCH --output=/home/chenlu/skx_cd/log/%A_%a.out

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cd_meta

MANIFEST=/home/chenlu/skx_cd/mapping_manifest.tsv
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")

sample=$(echo "$line" | cut -f1)
group=$(echo "$line" | cut -f2)
ref=$(echo "$line" | cut -f3)
r1=$(echo "$line" | cut -f4)
r2=$(echo "$line" | cut -f5)

mkdir -p /home/chenlu/skx_cd/{mapback_idx,mapback_bam,mapback_cov,tmp,log}

prefix=/home/chenlu/skx_cd/mapback_idx/${sample}

if [[ ! -s ${prefix}.1.bt2 ]]; then
    bowtie2-build "$ref" "$prefix"
fi

bowtie2 --very-sensitive-local \
    -p "${SLURM_CPUS_PER_TASK}" \
    -x "$prefix" \
    -1 "$r1" \
    -2 "$r2" \
    2> /home/chenlu/skx_cd/log/${sample}.mapback.bowtie2.log \
| samtools sort -@ 4 -O BAM \
    -T /home/chenlu/skx_cd/tmp/${sample} \
    -o /home/chenlu/skx_cd/mapback_bam/${sample}.sorted.bam -

samtools index /home/chenlu/skx_cd/mapback_bam/${sample}.sorted.bam

samtools coverage /home/chenlu/skx_cd/mapback_bam/${sample}.sorted.bam \
    > /home/chenlu/skx_cd/mapback_cov/${sample}.coverage.tsv
