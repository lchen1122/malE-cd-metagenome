#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# ============================================================
# Script: extract_target_orfs_with_assignment.py
# Purpose: Extract target ORFs from Prodigal genes.fna based on
#          DIAMOND hits and assign them to curated target families
# Stage: Annotation
# ============================================================
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract target ORFs from Prodigal genes.fna using DIAMOND hits, keep all qualified hits, and assign family."
    )
    p.add_argument("-d", "--diamond", required=True,
                   help="DIAMOND result TSV (outfmt 6: qseqid sseqid pident length evalue bitscore ...)")
    p.add_argument("-g", "--genes", required=True,
                   help="Prodigal genes.fna")
    p.add_argument("--out-all-fna", required=True,
                   help="Output FASTA of all unique target ORFs")
    p.add_argument("--out-confident-fna", required=True,
                   help="Output FASTA of confident target ORFs only")
    p.add_argument("--out-all-hits", required=True,
                   help="Output TSV of all filtered hits")
    p.add_argument("--out-assignment", required=True,
                   help="Output TSV of final ORF assignment")
    p.add_argument("--max-evalue", type=float, default=1e-5,
                   help="Maximum e-value to keep hit (default: 1e-5)")
    p.add_argument("--min-bitscore", type=float, default=50.0,
                   help="Minimum bitscore to keep hit (default: 50)")
    p.add_argument("--min-aln-len", type=int, default=60,
                   help="Minimum amino-acid alignment length to keep hit (default: 60)")
    p.add_argument("--min-pident", type=float, default=0.0,
                   help="Minimum percent identity to keep hit (default: 0)")
    p.add_argument("--bitscore-gap", type=float, default=10.0,
                   help="Minimum bitscore gap between top family and second family to call confident (default: 10)")
    return p.parse_args()


def family_from_target(sseqid: str) -> str:
    # 例如 MalE__1 -> MalE
    return sseqid.split("__")[0]


def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(re.sub(r"\s+", "", line))
        if header is not None:
            yield header, "".join(seq_chunks)


def prodigal_orf_id(header: str) -> str:
    # Prodigal header 常见格式：
    # >k141_62860_1 # 1 # 819 # 1 # ID=...
    # 取第一个空白前的 token 作为 ORF ID
    return header.split()[0]


def load_orf_sequences(genes_fna):
    seqs = {}
    for header, seq in fasta_iter(genes_fna):
        orf_id = prodigal_orf_id(header)
        seqs[orf_id] = seq
    return seqs


def load_filtered_hits(diamond_file, max_evalue, min_bitscore, min_aln_len, min_pident):
    """
    保留所有通过过滤条件的 hits
    返回: hits_by_orf[orf_id] = [hit1, hit2, ...]
    """
    hits_by_orf = defaultdict(list)

    with open(diamond_file) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 6:
                continue

            qseqid = row[0].strip()
            sseqid = row[1].strip()
            pident = float(row[2])
            aln_len = int(float(row[3]))
            evalue = float(row[4])
            bitscore = float(row[5])

            if evalue > max_evalue:
                continue
            if bitscore < min_bitscore:
                continue
            if aln_len < min_aln_len:
                continue
            if pident < min_pident:
                continue

            hit = {
                "orf_id": qseqid,
                "target": sseqid,
                "family": family_from_target(sseqid),
                "pident": pident,
                "aln_len": aln_len,
                "evalue": evalue,
                "bitscore": bitscore
            }
            hits_by_orf[qseqid].append(hit)

    # 每个 ORF 的 hits 排序：bitscore 高优先，其次 evalue 小优先
    for orf_id in hits_by_orf:
        hits_by_orf[orf_id].sort(key=lambda x: (-x["bitscore"], x["evalue"]))

    return hits_by_orf


def assign_family(hits, bitscore_gap):
    """
    对一个 ORF 的所有 hits 做 family assignment：
    1. 先按 family 聚合，每个 family 只保留最佳 hit
    2. 如果只有 1 个 family -> confident
    3. 如果 top family 比 second family 的 bitscore 至少高 bitscore_gap -> confident
    4. 否则 -> ambiguous
    """
    best_by_family = {}

    for h in hits:
        fam = h["family"]
        if fam not in best_by_family:
            best_by_family[fam] = h
        else:
            old = best_by_family[fam]
            if (h["bitscore"] > old["bitscore"]) or (
                h["bitscore"] == old["bitscore"] and h["evalue"] < old["evalue"]
            ):
                best_by_family[fam] = h

    fam_hits = list(best_by_family.values())
    fam_hits.sort(key=lambda x: (-x["bitscore"], x["evalue"]))

    top = fam_hits[0]

    if len(fam_hits) == 1:
        return {
            "assigned_family": top["family"],
            "status": "confident",
            "top_target": top["target"],
            "top_bitscore": top["bitscore"],
            "top_evalue": top["evalue"],
            "second_family": "",
            "second_bitscore": "",
            "bitscore_gap": "",
            "n_families": 1,
            "n_hits": len(hits)
        }

    second = fam_hits[1]
    gap = top["bitscore"] - second["bitscore"]

    if gap >= bitscore_gap:
        return {
            "assigned_family": top["family"],
            "status": "confident",
            "top_target": top["target"],
            "top_bitscore": top["bitscore"],
            "top_evalue": top["evalue"],
            "second_family": second["family"],
            "second_bitscore": second["bitscore"],
            "bitscore_gap": gap,
            "n_families": len(fam_hits),
            "n_hits": len(hits)
        }
    else:
        return {
            "assigned_family": "ambiguous",
            "status": "ambiguous",
            "top_target": top["target"],
            "top_bitscore": top["bitscore"],
            "top_evalue": top["evalue"],
            "second_family": second["family"],
            "second_bitscore": second["bitscore"],
            "bitscore_gap": gap,
            "n_families": len(fam_hits),
            "n_hits": len(hits)
        }


def write_fasta(orf_ids, seqs, out_fna):
    Path(out_fna).parent.mkdir(parents=True, exist_ok=True)
    with open(out_fna, "w") as out:
        for orf_id in sorted(orf_ids):
            if orf_id not in seqs:
                continue
            seq = seqs[orf_id]
            out.write(f">{orf_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def main():
    args = parse_args()

    seqs = load_orf_sequences(args.genes)
    hits_by_orf = load_filtered_hits(
        args.diamond,
        args.max_evalue,
        args.min_bitscore,
        args.min_aln_len,
        args.min_pident
    )

    # 只保留 genes.fna 里确实存在的 ORF
    valid_orfs = sorted([orf for orf in hits_by_orf if orf in seqs])

    missing = sorted(set(hits_by_orf) - set(valid_orfs))
    if missing:
        print(f"Warning: {len(missing)} ORFs found in DIAMOND hits but not found in genes.fna")
        for x in missing[:10]:
            print("  missing:", x)

    Path(args.out_all_hits).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_assignment).parent.mkdir(parents=True, exist_ok=True)

    # 1. all hits
    with open(args.out_all_hits, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "orf_id", "family", "target", "pident", "aln_len", "evalue", "bitscore", "orf_len_bp"
        ])
        for orf_id in valid_orfs:
            orf_len = len(seqs[orf_id])
            for h in hits_by_orf[orf_id]:
                writer.writerow([
                    orf_id,
                    h["family"],
                    h["target"],
                    f"{h['pident']:.2f}",
                    h["aln_len"],
                    f"{h['evalue']:.3g}",
                    f"{h['bitscore']:.2f}",
                    orf_len
                ])

    # 2. assignment
    assignment = {}
    with open(args.out_assignment, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "orf_id", "assigned_family", "status",
            "top_target", "top_bitscore", "top_evalue",
            "second_family", "second_bitscore", "bitscore_gap",
            "orf_len_bp", "n_hits", "n_families"
        ])

        for orf_id in valid_orfs:
            decision = assign_family(hits_by_orf[orf_id], args.bitscore_gap)
            decision["orf_len_bp"] = len(seqs[orf_id])
            assignment[orf_id] = decision

            writer.writerow([
                orf_id,
                decision["assigned_family"],
                decision["status"],
                decision["top_target"],
                f"{decision['top_bitscore']:.2f}",
                f"{decision['top_evalue']:.3g}",
                decision["second_family"],
                decision["second_bitscore"] if decision["second_bitscore"] != "" else "",
                f"{decision['bitscore_gap']:.2f}" if decision["bitscore_gap"] != "" else "",
                decision["orf_len_bp"],
                decision["n_hits"],
                decision["n_families"]
            ])

    # 3. fasta: all valid target ORFs
    write_fasta(valid_orfs, seqs, args.out_all_fna)

    # 4. fasta: confident only
    confident_orfs = [orf for orf in valid_orfs if assignment[orf]["status"] == "confident"]
    write_fasta(confident_orfs, seqs, args.out_confident_fna)

    print(f"Total ORFs passing filters: {len(valid_orfs)}")
    print(f"Confident ORFs: {len(confident_orfs)}")
    print(f"All target ORFs FASTA: {args.out_all_fna}")
    print(f"Confident ORFs FASTA:  {args.out_confident_fna}")
    print(f"All hits table:        {args.out_all_hits}")
    print(f"Assignment table:      {args.out_assignment}")


if __name__ == "__main__":
    main()
