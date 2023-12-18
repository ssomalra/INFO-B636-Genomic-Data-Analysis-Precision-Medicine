# Unraveling LARP4 Interactions with PABPs and Exploring Functional Significance in HepG2 Cells

Poly(A) binding proteins (PABPs) represent a crucial class of RNA binding proteins (RBPs) that play a pivotal role in mRNA regulation. Specifically, these proteins bind to the 3’ poly(A) tail region in transcripts, exerting control over mRNA decay and facilitating translation initiation. Among these, poly(A)-binding protein cytoplasmic 1 (PABPC1) is a specialized PABP known for its interaction with poly(A) tails in the cytoplasm. This interaction serves to safeguard mRNA from degradation, thereby promoting efficient mRNA translation.

In the intricate landscape of mRNA regulation, La‐related RBP 4 (LARP4) emerges as a key player. LARP4 functions by regulating mRNA translation through its interactions with the poly(A) tail and other PABPs, such as PABPC1. By doing so, LARP4 not only shields mRNA from cytoplasmic deadenylation but also contributes to poly(A) tail lengthening. Notably, disruptions in LARP4 have been associated with various cancers, emphasizing its critical role in cellular processes. Depletion of LARP4 has been shown to enhance cell migration and invasion, underscoring its significance in cancer biology.

Given the pivotal role of LARP4 in the regulation of mRNA and its implications in cancer, this research aims to utilize enhanced cross-linking and immunoprecipitation (eCLIP) data obtained from HepG2 cells targeting LARP4. The primary objective is to validate the binding interactions between LARP4, poly(A) tails, and potentially PABPs for advancing our understanding of the intricate molecular mechanisms governing mRNA regulation, particularly in the context of cancer cells.

## Data Retrieval
<b>eCLIP</b>
- [ENCSR805SRN](https://www.encodeproject.org/experiments/ENCSR805SRN/): 2 replicates of HepG2 Cells treated with LARP4
- [ENCSR192PUP](https://www.encodeproject.org/experiments/ENCSR192PUP/): Control dataset with no specific target in HepG2 cells

<b>RNAseq</b>
- [ENCSR744PAQ](https://www.encodeproject.org/experiments/ENCSR744PAQ/): 2 replicates of HepG2 cells subjected to shRNA knockdown of LARP4
- [ENCSR135LXL](https://www.encodeproject.org/experiments/ENCSR135LXL/): Control shRNA against no target in HepG2 cells

## eCLIP Analysis Execution
1. FASTQ files were downloaded from ENCODE.
2. Reads were checked for quality using FastQC.
3. Adapters and poor quality reads were removed using TrimGalore.
4. Genome was assembled using HISAT2. The reference genome was downloaded from Ensembl release 110 (Homo_sapiens.GRCh38.dna.primary_assembly.fa).
5. Preprocessed reads were aligned to the human GRCh38 genome with HISAT2.
6. SAMtools was employed to generate sorted BAM files.
7. The sorted BAM files served as input for peakcalling using mACS2 with the --broadpeaks parameter.
8. Consensus peaks, totaling 9,120 peaks, were by intersecting peaks identified in the two replicates using BEDtools.
9. Consensus peaks were used for peak annotations with Homer annotatePeaks.pl.
10. Peaks corresponding to the 3' UTR region were used for motif discovery using MEME and functional enrichment analysis with the R package, ClusterProfiler.

## RNAseq Analysis
1. FASTQ files were downloaded from ENCODE.
2. Reads were checked for quality using FastQC.
3. Adapters and poor quality reads were removed using TrimGalore.
4. Preprocessed reads were aligned to hte human GRCh38 genome with HISAT2.
5. The resulting SAM files were converted to sorted BAM files using SAMtools.
6. The sorted BAM files served as input for quantification of read expression using FeatureCounts.
7. The count table was downloaded and was used as input for differentially gene expression with the R package, DESeq2. Genes exhibiting differential expression were identified by a log2FoldChange less than -2 or greater than 2, coupled with a false discovery rate less than 0.01.
