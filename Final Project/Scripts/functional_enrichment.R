# 1. INSTALL PACKAGES ==========================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.17")
BiocManager::install("GO.db")
BiocManager::install("HDO.db", force=TRUE)
BiocManager::install("digest")
BiocManager::install("pathview")
BiocManager::install("enrichplot", force=TRUE)
library(clusterProfiler)
library(enrichplot)
library(ggplot2) # we use ggplot2 to add x axis labels (ex: ridgeplot)

# GET ANNOTATIONS
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# 2. FUNCTIONAL ENRICHMENT =====================================================
# open file
UTR_genes <- read.xlsx('Documents/IUPUI/23-24 Courses/636_Precision_Medicine/final_project/homer_ann.xlsx', sheet="3UTR_genes")
genes <- UTR_genes$Gene.Name # extract genes
genes <- na.omit(genes)
sorted_genes <- sort(genes, decreasing = TRUE)

# biological processes
functional_enrichment <- enrichGO(gene = sorted_genes,
                                  OrgDb = 'org.Hs.eg.db',
                                  keyType = 'SYMBOL',
                                  pAdjustMethod = 'none',
                                  ont = 'BP',
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)

plot(dotplot(functional_enrichment, showCategory=10))
