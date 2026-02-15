# Metagenomic Analysis Pipeline

This repository contains a complete workflow for metagenomic data analysis, including quality control, host removal, taxonomic classification, assembly, functional annotation, and MAG reconstruction.

---

## Requirements

- SeqKit
- Trimmomatic v0.38
- FastQC / MultiQC
- Bowtie2
- Samtools
- Bedtools
- Kraken2 / Bracken
- MEGAHIT
- QUAST
- Prokka
- eggNOG-mapper
- MetaBAT2
- CheckM
- KOfam

---

## Workflow

### 1. Quality Control
- Trimmomatic trimming
- FastQC / MultiQC quality assessment

### 2. Host Removal
- Bowtie2 mapping against host genome
- Extract unmapped reads

### 3. Taxonomic Classification
- Kraken2 classification
- Bracken abundance estimation

### 4. Assembly
- MEGAHIT metagenome assembly
- QUAST assembly quality assessment

### 5. Functional Annotation
- Prokka gene prediction
- eggNOG-mapper functional annotation (KEGG / GO)

### 6. MAG Reconstruction
- Bowtie2 read mapping
- MetaBAT2 binning
- CheckM quality assessment

### 7. MAG Functional Annotation
- Prokka annotation of MAGs
- KEGG annotation using KOfam

---

## Output

- Quality-controlled reads
- Host-removed reads
- Taxonomic profiles
- Assembled contigs
- Functional annotations (KEGG / GO)
- High-quality MAGs
- KEGG gene profiles
