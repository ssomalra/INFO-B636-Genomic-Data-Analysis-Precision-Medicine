# 0. INSTALL PACKAGES ==========================================================
install.packages('DESeq2')
install.packages('AnnotationDbi')
install.packages('org.Hs.eg.db')
install.packages('ggrepel')

library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggrepel)

# 1. PREPEARE DATA =============================================================
# 1. read in count table
count_matrix <- read.csv('~/Documents/IUPUI/23-24 Courses/636_Precision_Medicine/final_project/featureCounts.txt', sep='\t')

count_matrix <- count_matrix[, -c(2,3,4,5,6)] # remove columns from count_matrix
rownames(count_matrix) <- count_matrix$Geneid # set rownames
count_matrix <- count_matrix[, -1] # remove first column

# 2. read in sample info
metadata <- read.csv('~/Documents/IUPUI/23-24 Courses/636_Precision_Medicine/final_project/metadata.csv')

# add rownames to coldata
rownames(metadata) <- c('RNAseq_control1.sorted.bam', 'RNAseq_control2.sorted.bam',
                        'RNAseq_rep1.sorted.bam', 'RNAseq_rep2.sorted.bam')

# 2. Construct DESeq2DataSetObject =============================================
# Differentially expressed genes between control vs SARS-CoV-2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                        colData = metadata,
                                        design = ~ condition)

# pre-filtering: removing rows with low gene counts
# keeping rows with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 3. RUN DESEQ2 ================================================================
# condition
dds <- DESeq(dds)
res <- results(dds, alpha=0.01)

filtered_res <- na.omit(res)
filtered_res <- filtered_res[filtered_res$padj < 0.01 & abs(filtered_res$log2FoldChange) > 2, ]

res <- as.data.frame(res)
filtered_res <- as.data.frame(filtered_res)

# 4. GET GENE NAMES ============================================================
ensembl_id <- rownames(filtered_res) # extract ensemble IDs

# get gene names corresponding to ensemble IDs
geneSymbols <- select(org.Hs.eg.db, keys=ensembl_id, column='SYMBOL', keytype='ENSEMBL')
geneSymbols$ENSEMBL <- rownames(geneSymbols)

geneSymbols_unique <- geneSymbols %>% distinct(ENSEMBL, .keep_all = TRUE) # remove duplicates

filtered_res <- cbind(filtered_res, geneSymbols_unique[, 2]) # add gene names to dataframe
colnames(filtered_res)[ncol(filtered_res)] <- "geneName" # rename column

write.csv(filtered_res, file='~/Documents/IUPUI/23-24 Courses/636_Precision_Medicine/final_project/LARP4_DEG.csv')

# 5. VIOLIN PLOT ===============================================================
# prepare res dataset
ensemble_id_res <- rownames(res)
geneSymbols <- mapIds(org.Hs.eg.db, keys=ensemble_id_res, column='SYMBOL', keytype='ENSEMBL')
geneSymbols <- as.data.frame(geneSymbols)

res2 <- cbind(res, geneSymbols[, 1]) # add gene names to dataframe
colnames(res2)[ncol(res2)] <- "geneName" # rename column

res2 <- na.omit(res2) # remove NAs

# label genes
res3 <- res2 %>%
  mutate(
    Expression = case_when(log2FoldChange >= 2 & padj <= 0.01 ~ "Up-regulated",
                           log2FoldChange <= -2 & padj <= 0.01 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

#write.csv(res3, file='~/Documents/IUPUI/23-24 Courses/636_Precision_Medicine/final_project/LARP4_DEG_labelled.csv')

# extract top 10 upregulated and top 10 downregulated genes
top_upregulated <- res3 %>%
  filter(Expression == "Up-regulated") %>%
  top_n(10, wt = log2FoldChange * -log(padj, 10))

top_downregulated <- res3 %>%
  filter(Expression == "Down-regulated") %>%
  top_n(10, wt = -log2FoldChange * -log(padj, 10))

top_genes <- rbind(top_upregulated, top_downregulated)

# plot volcano plot
p <- ggplot(res3, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color=Expression), size=2/5) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c('dodgerblue3', 'gray50', 'firebrick3')) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_vline(xintercept=c(-2, 2), col="grey10", linetype = "dotted") +
  geom_hline(yintercept=-log10(0.01), col="grey10", linetype = 'dotted') +
  geom_text(data = top_genes, aes(label = geneName), size=2, vjust = -0.50, hjust = 0.50) +
  theme_minimal()