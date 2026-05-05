#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import sys
from pathlib import Path
from collections import defaultdict


TARGET_DIR = Path("/home/chenlu/skx_cd/target_orfs")
COV_DIR = Path("/home/chenlu/skx_cd/mapback_cov")
CONTROL_DIR = Path("/home/chenlu/skx_cd/control/data")
POLLUTION_DIR = Path("/home/chenlu/skx_cd/cd_pollution/data")
OUT_DIR = Path("/home/chenlu/skx_cd/family_abundance")

OUT_DIR.mkdir(parents=True, exist_ok=True)


def count_reads_from_fastq(fq_path: Path) -> int:
    n = 0
    with open(fq_path) as f:
        for i, _ in enumerate(f, 1):
            pass
    return i // 4 if 'i' in locals() else 0


def get_raw_fastqs(sample: str):
    raw = sample
    if sample.startswith("control_"):
        raw = sample[len("control_"):]
        r1 = CONTROL_DIR / f"{raw}_1.fastq"
        r2 = CONTROL_DIR / f"{raw}_2.fastq"
    elif sample.startswith("pollution_"):
        raw = sample[len("pollution_"):]
        r1 = POLLUTION_DIR / f"{raw}_1.fastq"
        r2 = POLLUTION_DIR / f"{raw}_2.fastq"
    else:
        raise ValueError(f"Unknown sample prefix: {sample}")
    return r1, r2


def load_assignment(assignment_file: Path):
    orf2family = {}
    orf2status = {}
    with open(assignment_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            orf_id = row["orf_id"]
            orf2family[orf_id] = row["assigned_family"]
            orf2status[orf_id] = row["status"]
    return orf2family, orf2status


def load_coverage_and_sum(coverage_file: Path, orf2family, orf2status, total_reads: int):
    fam_ab = defaultdict(float)
    fam_numreads = defaultdict(float)
    fam_orf_count = defaultdict(int)

    seen_orfs = set()

    with open(coverage_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            orf_id = row["#rname"]

            if orf_id not in orf2family:
                continue
            if orf2status.get(orf_id) != "confident":
                continue
            if orf2family.get(orf_id) == "ambiguous":
                continue

            family = orf2family[orf_id]

            start = int(row["startpos"])
            end = int(row["endpos"])
            numreads = float(row["numreads"])

            orf_len_bp = end - start + 1
            if orf_len_bp <= 0:
                continue

            # RPKM-like abundance
            abundance = numreads * 1e9 / (orf_len_bp * total_reads)

            fam_ab[family] += abundance
            fam_numreads[family] += numreads
            if orf_id not in seen_orfs:
                fam_orf_count[family] += 1
                seen_orfs.add(orf_id)

    return fam_ab, fam_numreads, fam_orf_count


def main():
    assignment_files = sorted(TARGET_DIR.glob("*.target_orfs.assignment.tsv"))

    long_rows = []
    all_families = set()
    total_reads_rows = []

    for assignment_file in assignment_files:
        sample = assignment_file.name.replace(".target_orfs.assignment.tsv", "")
        coverage_file = COV_DIR / f"{sample}.coverage.tsv"

        if not coverage_file.exists():
            print(f"[WARN] coverage missing for {sample}", file=sys.stderr)
            continue

        try:
            r1, r2 = get_raw_fastqs(sample)
        except ValueError as e:
            print(f"[WARN] {e}", file=sys.stderr)
            continue

        if not r1.exists() or not r2.exists():
            print(f"[WARN] raw fastq missing for {sample}", file=sys.stderr)
            continue

        print(f"[INFO] processing {sample}", file=sys.stderr)

        total_reads = count_reads_from_fastq(r1) + count_reads_from_fastq(r2)
        total_reads_rows.append([sample, total_reads, str(r1), str(r2)])

        orf2family, orf2status = load_assignment(assignment_file)
        fam_ab, fam_numreads, fam_orf_count = load_coverage_and_sum(
            coverage_file, orf2family, orf2status, total_reads
        )

        for family in sorted(fam_ab):
            long_rows.append([
                sample,
                family,
                f"{fam_ab[family]:.6f}",
                f"{fam_numreads[family]:.0f}",
                fam_orf_count[family],
                total_reads
            ])
            all_families.add(family)

    # 1. long table
    long_out = OUT_DIR / "family_abundance.long.tsv"
    with open(long_out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample", "family", "abundance_rpkm_like", "family_numreads", "family_orf_count", "total_reads"])
        writer.writerows(long_rows)

    # 2. matrix
    sample2fam = defaultdict(dict)
    for row in long_rows:
        sample, family, abundance = row[0], row[1], row[2]
        sample2fam[sample][family] = abundance

    matrix_out = OUT_DIR / "family_abundance.matrix.tsv"
    fam_list = sorted(all_families)
    with open(matrix_out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample"] + fam_list)
        for sample in sorted(sample2fam):
            writer.writerow([sample] + [sample2fam[sample].get(fam, "0") for fam in fam_list])

    # 3. total reads table
    total_reads_out = OUT_DIR / "sample_total_reads.tsv"
    with open(total_reads_out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample", "total_reads", "r1", "r2"])
        writer.writerows(total_reads_rows)

    print(f"[OK] written: {long_out}", file=sys.stderr)
    print(f"[OK] written: {matrix_out}", file=sys.stderr)
    print(f"[OK] written: {total_reads_out}", file=sys.stderr)


if __name__ == "__main__":
