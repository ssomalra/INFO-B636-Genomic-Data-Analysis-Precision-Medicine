#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --time=1-6:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=RNAseq
#SBATCH -o RNAseq.out
#SBATCH -A c00473

# fastq file information
# treatment - 2 replicates; paired_end
ENCFF597AES.fastq > RNAseq_rep1_1.fastq
ENCFF604CLV.fastq > RNAseq_rep1_2.fastq

ENCFF905PJJ.fastq > RNAseq_rep2_1.fastq
ENCFF996ORV.fastq > RNAseq_rep2_2.fastq

# control - paired_end
ENCFF226MZL.fastq > RNAseq_control1_1.fastq
ENCFF108KML.fastq > RNAseq_control1_2.fastq

ENCFF512QIR.fastq -> RNAseq_control2_1.fastq
ENCFF816KYC.fastq -> RNAseq_control2_2.fastq


# quality control - fastqc
module load fastqc
fastqc RNAseq_rep1_1.fastq RNAseq_rep1_2.fastq RNAseq_rep2_1.fastq RNAseq_rep2_2.fastq RNAseq_control1_1.fastq RNAseq_control1_2.fastq RNAseq_control2_1.fastq RNAseq_control2_2.fastq


# trim adaptors
pip install cutadapt
module load fastqc/0.11.9
module load python/3.10.5
module load trimgalore

trim_galore --phred33 --fastqc --paired RNAseq_rep1_1.fastq RNAseq_rep1_2.fastq
trim_galore --phred33 --fastqc --paired RNAseq_rep2_1.fastq RNAseq_rep2_2.fastq
trim_galore --phred33 --fastqc --paired RNAseq_control1_1.fastq RNAseq_control1_2.fastq
trim_galore --phred33 --fastqc --paired RNAseq_control2_1.fastq RNAseq_control2_2.fastq


# perform alignment
module load hisat2

hisat2 -x /N/slate/ssomalra/INFO_636/final_project/assembly/genome_110 -1 RNAseq_control1_1_val_1.fq -2 RNAseq_control1_2_val_2.fq -S RNAseq_control1.sam --summary-file RNAseq_control1_alignment.txt
hisat2 -x /N/slate/ssomalra/INFO_636/final_project/assembly/genome_110 -1 RNAseq_control2_1_val_1.fq -2 RNAseq_control2_2_val_2.fq -S RNAseq_control2.sam --summary-file RNAseq_control2_alignment.txt
hisat2 -x /N/slate/ssomalra/INFO_636/final_project/assembly/genome_110 -1 RNAseq_rep1_1_val_1.fq -2 RNAseq_rep1_2_val_2.fq -S RNAseq_rep1.sam --summary-file RNAseq_rep1_alignment.txt
hisat2 -x /N/slate/ssomalra/INFO_636/final_project/assembly/genome_110 -1 RNAseq_rep2_1_val_1.fq -2 RNAseq_rep2_2_val_2.fq -S RNAseq_rep2.sam --summary-file RNAseq_rep2_alignment.txt


# binary index SAM
module unload gcc
module load gcc/9.3.0
module unload python
module load python/3.10.5
module load samtools


samtools view -S -b RNAseq_control1.sam > RNAseq_control1.bam
samtools view -S -b RNAseq_control2.sam > RNAseq_control2.bam
samtools view -S -b RNAseq_rep1.sam > RNAseq_rep1.bam
samtools view -S -b RNAseq_rep2.sam > RNAseq_rep2.bam


# sort BAM file
samtools sort RNAseq_control1.bam -o RNAseq_control1.sorted.bam
samtools sort RNAseq_control2.bam -o RNAseq_control2.sorted.bam
samtools sort RNAseq_rep1.bam -o RNAseq_rep1.sorted.bam
samtools sort RNAseq_rep2.bam -o RNAseq_rep2.sorted.bam


# download gtf file
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip gunzip Homo_sapiens.GRCh38.110.gtf.gz


# feature counts
module load subread
featureCounts -M -T 12 -a Homo_sapiens.GRCh38.110.gtf -o featureCounts.txt alignment/*.sorted.bam
