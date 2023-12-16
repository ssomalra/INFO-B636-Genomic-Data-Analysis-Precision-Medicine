#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-6:59:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=de_novo_assembly
#SBATCH -o de_novo_assembly.out
#SBATCH -A c00473

# load fastq file
module load sra-toolkit/3.0.0
prefetch SRR21904868


# quality control - fastqc
module load fastqc
fastqc SRR21904868_1.fastq SRR21904868_2.fastq


# trimming - trimgalore
module load trimgalore
trim_galore --phred33 --fastqc --paired SRR21904868_1.fastq SRR21904868_2.fastq


# de novo assembly - velvet
# build the hash index
velveth <kmer_size>_assembly_results <kmer size> -short -separate -fastq SRR21904868_1_val_1.fq SRR21904868_2_val_2.fq

# run assembly
velvetg <kmer_size>_assembly_results/ -read_trkg yes


# de novo assembly - oases
oases <kmer_size>_assembly_results/

