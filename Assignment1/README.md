# Assignment 1 Report
When a reference genome is not present, reads must be aligned using de novo assembly. A comparative analysis was conducted on two de novo assembly tools, Velvet and Oases, both of which are based on the de Bruijn graph algorithm. This analysis was performed using approximately 2 million short reads obtained from Illumina DNA sequencing data of Escherichia coli. 

Fastq files were downloaded from the Sequence Read Archive (SRA) with the accession number SRR21904868. These files underwent preprocessing, which involved adapter trimming and the removal of low-quality reads before assembly. Roughly about 5,000 reads were removed during this process. Both Velvet and Oases were executed with the recommended parameters, with the exception that Velvetg employed the read tracking parameter to offer a more detailed description of the assembly. This parameter is crucial for Oases' algorithm as Oases relies on Velvet for its operations. The two tools were run using seven different kmer lengths: 75, 95, 97, 99, 101, 105, 115.

Velvet assembly process involves two steps: the creation of a hash index using the Fastq files, followed by the actual assembly. Oases is reliant on the preliminary assembly produced by the Velvet assembler. 

For the seven different kmer sizes, a comparison was made between the two tools using several assembly metrics, including the N50 score, number of contigs, length of the longest contig, and the total length of assembled contigs. Velvet provided these assembly metrics upon the completion of its run. In contrast, Oases did not directly output these metrics. Instead, it generated a file named ‘transcripts.fa’, which contained the sequences of all contigs in FASTA format. The header of each contig included its length, which was extracted and used to calculate the assembly metrics through a Python script.

<p align="center">
  <img width="793" alt="velvet_oases_n50" src="https://github.com/ssomalra/INFO-B636-Genomic-Data-Analysis-Precision-Medicine/assets/116883221/ed2cdd9d-de8e-4881-92e2-7cb8b91052ce">
</p>

Figure 1 visually depicts the impact of kmer length on the N50 score for both Velvet and Oases. The data indicates that, for both tools, the N50 score reached its peak at a kmer length of 97, which seems to be the optimal choice considering the specific parameters used. Notably, Oases consistently achieved a N50 score higher than Velvet, with values of 130,000 and 70,000, respectively, at a kmer size of 97. Moreover, Oases produced assemblies with fewer, longer contigs, resulting in an overall longer assembled genome compared to Velvet, as detailed in Table 1.

<p align="center">
  <b>Table 1. Comparison of Velvet and Oases Assembly Metrics</b>
  <img width="629" alt="assembler_metric_results" src="https://github.com/ssomalra/INFO-B636-Genomic-Data-Analysis-Precision-Medicine/assets/116883221/9a4b5a84-544a-45e5-bc7c-6d5c8ad8a374">
</p>

The justification for the much better performance of Oases can be attributed to the corrective measures it applies after using Velvet’s assembly as an input. For instance, Oases uses algorithms similar to TourBus searches to help identify and correct potentially erroneous chimeric contigs. Chimeric contigs can introduce errors in genome assembly, making their correction essential. This distinction could explain why Oases consistently produces superior assembly metrics compared to Velvet, establishing it as the better tools for de novo assembly. 
