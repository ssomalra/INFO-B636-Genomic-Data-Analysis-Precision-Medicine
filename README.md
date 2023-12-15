# Genomic-Data-Analysis-Precision-Medicine
This repository contains the code for assignments from the INFO-B 636: Genomic Data Analysis & Precision Medicine course.

## Assignment 1: De Novo Assembly
A comparative analysis of the de novo assembly tools, Velvet and Oases, in their ability to assemble a bacterial genome given short-read DNA sequencing data for Escherichia coli.

Programming Language: Python and Linux

Required Files:
- de_novo_assembly.sh
- assembly_metrics.ipynb

Execution:
1. Fastq files were downloaded from SRA (SRR21904868) 
2. Reads were checked for quality using FASTQC
3. Adapters and poor quality reads were removed using Trimgalore
4. De novo assembly was ran using Velvet and Oases. 
5. The output file from Oases, transcript.fa, was downloaded and used as input in assembly_metrics.ipynb to calculate assembly metrics. 

## Assignment 2: RNA Sequencing
RNA-seq analysis to evaluate changes in transcriptomic profiles between the two conditions, SARS-CoV-2 infected vs non-infected samples, and the two time points, 24H and 72H. 

Programming Language: R and Linux

Required Files:
- assignment2.sh
- DESeq2.R
- assignment_2_info.csv : metadata

Softwares Used:
- RStudio
- Cytoscape with ClueGo add-on

Output Files:
- featureCounts.txt : count table of expression levels between samples
- deseq_treatment_results.csv : DEGs between the two conditions
- deseq_timepoint_results.csv : DEGs between the two time points
- DEG_timepoint_GOterms.xlsx : GO Term enrichment analysis for DEGs between the two time points

Execution:
1. Fastq files were downloaded from SRA (PRJNA901149)
2. Reads were checked for quality using FASTQC
3. Adapters and poor quality reads were removed using Trimgalore
4. Genome was assembled using STAR. The reference genome and gene annotation files were downloaded from Gencode release 44: GRCh38.primary_assembly.genome.fa and gencode.v44.annotation.gtf. The assembly represents the primary assembly of the Homo sapiens genome (GRCh38).
5. Trimmed reads were aligned using STAR
6. Reads were quantified using featureCounts.
7. The count table generated from the previous step along with the assignment_2_info.csv file, which contains the metadata, served as input for the DESeq2.R script for conducting a differential gene expression analysis through DESeq2.
9. Differentially expressed genes (DEGs) were used as input for functional enrichment analysis using ClueGo. NOTE: ClueGo did not find GO terms for the DEGs from the two conditions. 
