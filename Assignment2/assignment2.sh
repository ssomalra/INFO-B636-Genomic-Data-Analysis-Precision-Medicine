#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-6:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=assignment2
#SBATCH -o assignment2.out
#SBATCH -A c00473

# download files
module load sra-toolkit/3.0.0
for ((i=72;i<=83;i++));do fasterq-dump SRR222698$i;done


# quality control - fastqc
module load fastqc
for ((i=72;i<=83;i++));do fastqc SRR222698$i.fastq;done


# trim adaptors
pip install cutadapt
module load fastqc/0.11.9
module load python/3.10.5
module load trimgalore
for ((i=72;i<=83;i++));do trim_galore --phred33 --illumina -q 30 SRR222698$i_trimmed.fq;done


# download human reference genome and gene annotation file from gencode 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz 
gunzip GRCh38.primary_assembly.genome.fa.gz


wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz


# assemble genome
module load star
STAR --runThreadN 12 --runMode genomeGenerate --genomeSAsparseD 2 --genomeDir GRCh38_index/release_44 --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v44.annotation.gtf


# alignment
for ((i=72;i<=83;i++));do STAR --runThreadN 12 --runMode alignReads --genomeDir GRCh38_release_44_assembly --readFilesIn SRR222698${i}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix SRR222698${i}.sam; done


# feature counts
module load subread
featureCounts -M -T 12 -a gencode.v44.annotation.gtf -o featureCounts.txt *.bam

