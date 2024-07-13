# Assignment 2
This assignment introduces students to the RNAseq computational pipeline.

## Assignment Description
SARS-CoV-2 is the respiratory virus responsible for the COVID-19 pandemic. Viral disease pathology can be studied by evaluating changes in transcriptomic profiles between infected and non-infected samples. Assessing changes at different time points can also illuminate how the immune system coordinates its response to viral infection. The data for this assignment focuses on miRNA expression levels in cultivated human respiratory cells: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA901149

In this assignment, you will perform transcriptome analysis for samples across:
- two condition sets : Mock Control and SARS-CoV-2 Infected
- two time points : 24H and 72H

## Assignment Goals
Metadata for the assignment is located in assignment_2_info.csv.
- Use this data to align, quantify and construct transcriptome profiles for control and SARS-CoV-2 infected respiratory cells for each time point.  
- Identify differentially-expressed genes (DEGS) between both (condition at 24 vs 72 H) and (control vs SARS-CoV-2). 
- Perform GO Term enrichment analysis on each set of DEGs. 

Provide every step involved in your analysis in README.txt file along with the software/s used and the version of the human genome build you are using for this assignment. Your output should include a text or excel file with the expression levels of every gene in the human genome across each condition/sample in each of these two sets.

## README.txt Submission
_**Description:**_ The steps involved in running RNA-seq analysis to evaluate changes in transcriptomic profiles between the two conditions, SARS-CoV-2 infected vs non-infected samples, and the two time points, 24H and 72H.

_**Required Files:**_
- assignment2.sh
- DESeq2.R

_**Softwares Used:**_
- RStudio
- Cytoscape with ClueGo plugin

_**Output Files:**_
- featureCounts.txt : count table of expression levels between samples
- deseq_treatment_results.csv : DEGs between the two conditions
- deseq_timepoint_results.csv : DEGs between the two time points
- DEG_timepoint_GOterms.xlsx : GO Term enrichment analysis for DEGs between the two time points

_**Execution:**_
1. Fastq files were downloaded from SRA (PRJNA901149)
2. Reads were checked for quality using FASTQC
3. Adapters and poor quality reads were removed using Trimgalore
4. Genome was assembled using STAR. The reference genome and gene annotation files were downloaded from Gencode release 44: GRCh38.primary_assembly.genome.fa and gencode.v44.annotation.gtf. The assembly represents the primary assembly of the Homo sapiens genome (GRCh38).
5. Trimmed reads were aligned using STAR
6. Reads were quantified using featureCounts.
7. The count table outputted from the previous step was downloaded and used as input for the DESeq2.R script for differential gene expression analysis using DESeq2.
8. Differentially expressed genes (DEGs) were used as input for functional enrichment analysis using ClueGo. NOTE: ClueGo did not find GO terms for the DEGs from the two conditions.
