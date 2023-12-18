#!/bin/bash
#SBATCH --mail-user=ssomalra@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --time=1-6:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=eCLIP
#SBATCH -o eCLIP.out
#SBATCH -A c00473

# fastq file information
# treatment - 2 replicates; paired_end
ENCFF365UUW.fastq > LARP4_rep1_1.fastq
ENCFF276CTM.fastq > LARP4_rep1_2.fastq

ENCFF362EMN.fastq > LARP4_rep2_1.fastq
ENCFF163HLQ.fastq > LARP4_rep2_2.fastq

# control - paired_end
ENCFF512AHQ.fastq > LARP4_control_1.fastq
ENCFF649CHF.fastq > LARP4_control_2.fastq


# quality control
module load fastqc
fastqc LARP4_rep1_1.fastq LARP4_rep1_2.fastq LARP4_rep2_1.fastq LARP4_rep2_2.fastq LARP4_control_1.fastq LARP4_control_2.fastq


# trimming
module load python/3.9.8
​module load trimgalore
trim_galore --phred33 --fastqc --paired LARP4_rep1_1.fastq LARP4_rep1_2.fastq
trim_galore --phred33 --fastqc --paired LARP4_rep2_1.fastq LARP4_rep2_2.fastq
trim_galore --phred33 --fastqc --paired LARP4_control_1.fastq LARP4_control_2.fastq


# download reference genome
# download genome from ensembl: https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


#Index genome
module load hisat2
mkdir assembly
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa assembly/genome_110


# perform alignment
module load hisat2

hisat2 -x assembly/genome_110 -1 LARP4_control_1_val_1.fq -2 LARP4_control_2_val_2.fq -S LARP4_control.sam --summary-file LARP4_control_alignment.txt
hisat2 -x assembly/genome_110 -1 LARP4_rep1_1_val_1.fq -2 LARP4_rep1_2_val_2.fq -S LARP4_rep1.sam --summary-file LARP4_rep1_alignment.txt
hisat2 -x assembly/genome_110 -1 LARP4_rep2_1_val_1.fq -2 LARP4_rep2_2_val_2.fq -S LARP4_rep2.sam --summary-file LARP4_rep2_alignment.txt


# binary index SAM
module unload gcc
module load gcc/9.3.0
module unload python
module load python/3.10.5
module load samtools

samtools view -S -b LARP4_control.sam > LARP4_control.bam
samtools view -S -b LARP4_rep1.sam > LARP4_rep1.bam
samtools view -S -b LARP4_rep2.sam > LARP4_rep2.bam


# sort BAM file
samtools sort LARP4_control.bam -o LARP4_control.sorted.bam
samtools sort LARP4_rep1.bam -o LARP4_rep1.sorted.bam
samtools sort LARP4_rep2.bam -o LARP4_rep2.sorted.bam


# index BAM file
samtools index -b LARP4_control.sorted.bam LARP4_control.sorted.bam.bai
samtools index -b LARP4_rep1.sorted.bam LARP4_rep1.sorted.bam.bai
samtools index -b LARP4_rep2.sorted.bam LARP4_rep2.sorted.bam.bai


# peakcalling
module unload python
module load macs
macs2 callpeak -t alignment/bam_files/LARP4_rep1.sorted.bam -c alignment/bam_files/LARP4_control.sorted.bam -f BAM -n LARP4_rep1 --nomodel --broad --outdir LARP4_rep1_peaks

macs2 callpeak -t alignment/bam_files/LARP4_rep2.sorted.bam -c alignment/bam_files/LARP4_control.sorted.bam -f BAM -n LARP4_rep2 --nomodel --broad --outdir LARP4_rep2_peaks


# consensus peaks
module load bedtools
bedtools intersect -a LARP4_rep1_peaks/LARP4_rep1_peaks.broadPeak -b LARP4_rep2_peaks/LARP4_rep1_peaks.broadPeak > LARP4_consensus_peaks.bed


#Annotate Peaks
#Before running ~ modify bed file so that the first column (chromosome name) begins with 'chr'
awk -F'\t' '{$1="chr"$1; print}' OFS='\t' LARP4_consensus_peaks.bed > LARP4_consensus_peaks2.bed 

module load homer/4.11.1 
annotatePeaks.pl differential_peaks/diff_c1_vs_c2_c3.0_cond2.bed hg38 -annStats annotation_stats.txt > homer_ann.txt


# functional enrichment
# from homer_ann.txt, I filtered the annotations for 3'UTR, then took all corresponding genes for functional enrichment analysis using the R package, ClusterProfiler


# filter 3'UTR peaks taken from homer_ann.txt
grep -Fwf 3_utr_peaks.txt LARP4_consensus_peaks.bed > 3utr_LARP4_peaks.bed


# extracting the sequences for the peak coordinates
module load bedtools
bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed peakcalling/3utr_LARP4_peaks.bed -s -fo 3utr_LARP4_sequences.fa


# motif discovery
module load perl/5.30.1
module load python/3.9.8
module load meme

meme 3utr_LARP4_sequences.fa -dna -maxw <n> -nmotifs 5 -o 3utr_motif_<n>
