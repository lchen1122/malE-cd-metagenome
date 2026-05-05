#!/bin/bash
#SBATCH --job-name=meta_asm
#SBATCH --partition=Par-2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --error=/home/chenlu/skx_cd/log/%A_%a.err
#SBATCH --output=/home/chenlu/skx_cd/log/%A_%a.out

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate meta_asm

manifest=/home/chenlu/skx_cd/sample_manifest.tsv
db=/home/chenlu/skx_cd/diamond_db/targets.final.dmnd

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")

sample=$(echo "$line" | cut -f1)
group=$(echo "$line" | cut -f2)
r1=$(echo "$line" | cut -f3)
r2=$(echo "$line" | cut -f4)

echo "Sample: $sample"
echo "Group: $group"
echo "R1: $r1"
echo "R2: $r2"

mkdir -p /home/chenlu/skx_cd/{assembly,prodigal,diamond_result,log,tmp}

# 1. assembly
if [[ ! -s /home/chenlu/skx_cd/assembly/${sample}/final.contigs.fa ]]; then
    megahit \
      -1 "$r1" \
      -2 "$r2" \
      -o /home/chenlu/skx_cd/assembly/${sample} \
      -t "${SLURM_CPUS_PER_TASK}"
else
    echo "[INFO] Assembly exists, skip: $sample"
fi

# 2. prodigal
if [[ ! -s /home/chenlu/skx_cd/prodigal/${sample}.proteins.faa ]]; then
    prodigal \
      -i /home/chenlu/skx_cd/assembly/${sample}/final.contigs.fa \
      -a /home/chenlu/skx_cd/prodigal/${sample}.proteins.faa \
      -d /home/chenlu/skx_cd/prodigal/${sample}.genes.fna \
      -o /home/chenlu/skx_cd/prodigal/${sample}.prodigal.gff \
      -p meta
else
    echo "[INFO] Prodigal result exists, skip: $sample"
fi

# 3. diamond blastp
if [[ ! -s /home/chenlu/skx_cd/diamond_result/${sample}.tsv ]]; then
    diamond blastp \
      -q /home/chenlu/skx_cd/prodigal/${sample}.proteins.faa \
      -d "$db" \
      -o /home/chenlu/skx_cd/diamond_result/${sample}.tsv \
      -p "${SLURM_CPUS_PER_TASK}" \
      --sensitive \
      -e 1e-5 \
      -k 1 \
      --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen
else
    echo "[INFO] DIAMOND result exists, skip: $sample"
fi
