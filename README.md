# malE-cd-metagenome

Computational workflows for metagenomic assembly, target gene annotation, abundance quantification, and association analysis of **malE** and curated cadmium-resistance genes in environmental metagenomic datasets.

## Overview

This repository contains the analysis pipeline used to evaluate the environmental relevance of **malE** and a curated cadmium-resistance gene set in publicly available metagenomic datasets from cadmium-contaminated and control environments.

The workflow includes:

1. de novo assembly of metagenomic reads  
2. ORF prediction from assembled contigs  
3. target gene annotation against a curated cadmium-resistance reference set  
4. extraction of high-confidence target ORFs  
5. read mapping back to sample-specific target references  
6. family-level abundance quantification  
7. downstream association analysis in R

## Repository structure

```text
malE-cd-metagenome/
├── README.md
├── sample_metadata.tsv
└── scripts/
    ├── assembly/
    │   └── assembly_diamond_array.sh
    ├── annotation/
    │   ├── run_extract_target_orfs_with_assignment.sh
    │   └── extract_target_orfs_with_assignment.py
    ├── mapping/
    │   ├── make_mapping_manifest.sh
    │   └── map_back_array.sh
    ├── quantification/
    │   └── summarize_all_family_abundance.py
    └── statistics/
        └── malE_cdset_association_analysis.R
```
## Sample metadata

Sample-level metadata are provided in `sample_metadata.tsv`, including the following fields:

- `Sample`
- `Group`
- `Reads`

`Group` indicates whether the sample belongs to the control or cadmium-contaminated dataset.

## Workflow

### 1. Assembly and target annotation

Metagenomic paired-end reads were assembled on a per-sample basis using **MEGAHIT**. ORFs were predicted from assembled contigs using **Prodigal** in metagenomic mode. Predicted proteins were searched against a curated reference set of **malE** and cadmium-resistance-related genes using **DIAMOND blastp**.

**Main script**

- `scripts/assembly/assembly_diamond_array.sh`

### 2. Extraction of high-confidence target ORFs

Hits passing predefined thresholds for e-value, bitscore, and alignment length were retained. High-confidence target ORFs were extracted from the corresponding nucleotide gene sequences to construct sample-specific target reference sets.

**Main scripts**

- `scripts/annotation/run_extract_target_orfs_with_assignment.sh`
- `scripts/annotation/extract_target_orfs_with_assignment.py`

### 3. Mapping back to sample-specific target references

Raw reads from each sample were mapped back to their own target ORF reference set using **Bowtie2**. Alignments were processed and coverage was calculated with **SAMtools**.

**Main scripts**

- `scripts/mapping/make_mapping_manifest.sh`
- `scripts/mapping/map_back_array.sh`

### 4. Family-level abundance quantification

ORF-level read counts were normalized by gene length and total sample read number to generate family-level relative abundance estimates.

**Main script**

- `scripts/quantification/summarize_all_family_abundance.py`

### 5. Downstream association analysis

Downstream statistical analysis and visualization were performed in **R**. The main analysis focused on the relationship between **malE** abundance and a curated cadmium-resistance gene set.

**Main script**

- `scripts/statistics/malE_cdset_association_analysis.R`

## Software

This workflow uses the following software:

- **MEGAHIT**
- **Prodigal**
- **DIAMOND**
- **Bowtie2**
- **SAMtools**
- **R**

## Inputs

The main inputs include:

- paired-end metagenomic FASTQ files
- curated target gene reference database
- sample metadata table

## Outputs

Typical outputs include:

- assembled contigs
- predicted protein and gene sequences
- DIAMOND annotation tables
- extracted target ORF FASTA files
- mapping coverage tables
- family-level abundance tables
- downstream statistical results and figures

## Notes

- The curated cadmium-resistance gene set used in this repository was defined prior to downstream analysis.
- Sample-specific target ORF references were used for read mapping to reduce bias associated with direct mapping to short external reference sequences.
- For visualization, abundance values were log10(x + 1)-transformed where appropriate.
- Correlation analyses were primarily based on **Spearman correlation**, while interaction models were used as supporting analyses to evaluate group-dependent trends.

## Citation

If you use this workflow, please cite the associated study and acknowledge the original public metagenomic datasets analyzed in this project.
