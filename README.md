# **Fusarium Evolution Pipeline**  
A reproducible phylogenomics workflow for extracting single-copy orthogroups from OrthoFinder results and constructing species phylogenies.

## Overview  
This repository provides an end-to-end pipeline for:  
1.Extracting single-copy orthogroups (SC-OGs)  
2.Performing multiple sequence alignment and trimming  
3.Constructing concatenated supermatrices  
4.Inferring phylogenetic trees  
The workflow is implemented using Nextflow, ensuring scalability and reproducibility across computing environments.  

## Workflow Description  
The pipeline consists of the following steps:  
### 1.Identification of single-copy orthogroups  
Based on OrthoFinder outputs:  
Orthogroups.tsv  
Orthogroups_SingleCopyOrthologues.txt    
### 2.Sequence extraction  
Protein sequences for each orthogroup are retrieved from species proteomes.
### 3.Quality control of orthogroups  
Detection of missing taxa  
Detection of duplicated genes  
### 4.Multiple sequence alignment  
Using MAFFT  
### 5.Alignment trimming
Using trimAl to remove poorly aligned regions
### 6.Alignment validation
Ensuring equal sequence lengths within each alignment
### 7.Concatenation
All alignments are concatenated into a supermatrix
→ Partition file is generated simultaneously
### 8.Phylogenetic inference (optional)
Using IQ-TREE2 with model selection and bootstrap support

## Requirements  
### Core dependencies  
Python ≥ 3.8  
Nextflow  
MAFFT  
trimAl  
### Python packages  
pandas  
biopython  
### Optional
IQ-TREE2  

## Usage
### Basic Run
    nextflow run main.nf \
      --orthogroups_tsv /path/to/Orthogroups.tsv \
      --single_copy_txt /path/to/Orthogroups_SingleCopyOrthologues.txt \
      --proteomes_dir /path/to/proteomes \
      --species_ref /path/to/species.txt \
      --expect_species_n N \
      --total_species N \
      --min_coverage 0.5 \
      --outdir results \
      --do_iqtree false    
### Run with phylogenetic inference
    nextflow run main.nf \
      --orthogroups_tsv /path/to/Orthogroups.tsv \
      --single_copy_txt /path/to/Orthogroups_SingleCopyOrthologues.txt \
      --proteomes_dir /path/to/proteomes \
      --species_ref /path/to/species.txt \
      --do_iqtree true \
      --iqtree_extra "-T AUTO"

## Output Structure
    results/
    ├── aln/                      # MAFFT alignments
    ├── trimmed/                 # trimAl outputs
    ├── concat/
    │   ├── concatenated_alignment.fasta
    │   └── partition.tsv
    ├── iqtree/                  # phylogenetic results (optional)

## Reproducibility  
This workflow is designed for fully reproducible phylogenomics analysis:  
All steps are scripted and version-controlled   
Nextflow ensures:  
1.parallel execution  
2.resume capability  
3.environment consistency  

## Notes & Best Practices
Protein FASTA files must follow naming convention:
`SpeciesName.xxx`  
→ prefix before . is used as species ID  
Species names must match exactly with:  
`species.txt`  
Current pipeline does not filter short OGs before concatenation  
Default IQ-TREE parameters:  
`-m MFP -bb 1000 -alrt 1000`
