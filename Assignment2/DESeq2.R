# 0. INSTALL PACKAGES ==========================================================
install.packages('DESeq2')
install.packages('AnnotationDbi')
install.packages('org.Hs.eg.db')

# 1. PREPEARE DATA =============================================================
# 1. read in count table
count_matrix <- read.csv('~/featureCounts.txt', sep='\t')

count_matrix <- count_matrix[, -c(2,3,4,5,6)] # remove columns from count_matrix
rownames(count_matrix) <- count_matrix$Geneid # set rownames
count_matrix <- count_matrix[, -1] # remove first column

# 2. read in sample info
coldata <- read.csv('~/assignment_2_info.csv')

# add rownames to coldata
rownames(coldata) <- c('alignment.SRR22269872.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269873.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269874.samAligned.sortedByCoord.out.bam',
                       'alignment.SRR22269875.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269876.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269877.samAligned.sortedByCoord.out.bam',
                       'alignment.SRR22269878.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269879.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269880.samAligned.sortedByCoord.out.bam',
                       'alignment.SRR22269881.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269882.samAligned.sortedByCoord.out.bam', 'alignment.SRR22269883.samAligned.sortedByCoord.out.bam')


# 2. Construct DESeq2DataSetObject =============================================
# Differentially expressed genes between control vs SARS-CoV-2
dds_treatment <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ Condition)

# Differentially expressed genes between condition at 24 vs 72H
dds_timepoint <- DESeqDataSetFromMatrix(countData = count_matrix,
                                        colData = coldata,
                                        design = ~ Time.Point)

# pre-filtering: removing rows with low gene counts
# keeping rows with at least 10 reads total
# Condition
keep <- rowSums(counts(dds_treatment)) >= 10
dds_treatment <- dds_treatment[keep, ]

# Time Point
keep <- rowSums(counts(dds_timepoint)) >= 10
dds_timepoint <- dds_timepoint[keep, ]

# 3. RUN DESEQ2 ================================================================
# condition
dds_treatment <- DESeq(dds_treatment)
res <- results(dds_treatment, alpha=0.01)

filtered_res <- na.omit(res)
filtered_res <- filtered_res[filtered_res$padj < 0.01 & abs(filtered_res$log2FoldChange) > 2, ]

filtered_res <- as.data.frame(filtered_res)

# time point
dds_timepoint <- DESeq(dds_timepoint)
res_time <- results(dds_timepoint, alpha=0.01)

filtered_res_time <- na.omit(res_time)
filtered_res_time <- filtered_res_time[filtered_res_time$padj < 0.01 & abs(filtered_res_time$log2FoldChange) > 0.02, ]

filtered_res_time <- as.data.frame(filtered_res_time)

# 4. GET GENE NAMES ============================================================
# condition
ensembl_id <- rownames(filtered_res) # extract ensemble IDs
ensembl_id <- sub("\\.\\d+", "", ensembl_id) # remove decimal points

# get gene names corresponding to ensemble IDs
geneSymbols <- select(org.Hs.eg.db, keys=ensembl_id, column='SYMBOL', keytype='ENSEMBL')

filtered_res <- cbind(filtered_res, geneSymbols[, 2]) # add gene names to dataframe
colnames(filtered_res)[ncol(filtered_res)] <- "geneName" # rename column

write.csv(filtered_res, file='~/Assignment_2/deseq_condition_results.csv')

# time point
ensembl_id_time <- rownames(filtered_res_time)
ensembl_id_time <- sub("\\.\\d+", "", ensembl_id_time)

geneSymbols_time <- select(org.Hs.eg.db, keys=ensembl_id_time, column='SYMBOL', keytype='ENSEMBL')

filtered_res_time <- cbind(filtered_res_time, geneSymbols_time[, 2])
colnames(filtered_res_time)[ncol(filtered_res_time)] <- "geneName"

write.csv(filtered_res_time, file='~/Assignment_2/deseq_timepoint_results.csv')

