# wczytanie paczek

library(DESeq2)
library(pheatmap)
library(ggplot2)


# read in data

args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
instrument = sub("\\_counts.txt$", "", counts_file)

count_matrix <- read.table(counts_file, header=TRUE, row.names=1)

count_data <- count_matrix[,6:length(count_matrix[1,])]
colnames(count_data) <-  sub("\\.bam$", "", colnames(count_data))


# metadane

sample_names <- names(count_data)

condition <- as.factor(
  c("mock", "mock",
    "zika", "zika"))

meta_data <- data.frame(samples=sample_names, 
                        condition=condition)


# objekt DESeq 
dds <- DESeqDataSetFromMatrix(countData=count_data, 
                              colData=meta_data, 
                              design = ~ condition)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)


# Normalizacja
log_data <- rlog(dds)
norm_data_matrix <- assay(log_data)
norm_data <- as.data.frame(norm_data.matrix)



# Heatmapa

my_heatmap <- pheatmap(as.matrix(norm_data),
        name = instrument,
        cluster_cols=FALSE,
        show_rownames=FALSE,
        column_labels = colnames(count_data),
        row_title = "Genes")

ggsave(paste("Projekt/DE_output/", instrument, "_heatmap.png", sep = ""), 
       my_heatmap, 
       device = "png",
       width=5,
       height=5)


# Analiza głównych składowych (PCA)

my_PCA <- plotPCA("Projekt/DE_output/",log_data, intgroup="condition")
ggsave(paste(instrument, "_pca.png", sep = ""), 
       my_PCA, 
       device = "png",
       width=5,
       height=5)

# istotne różnice w eksprsji genów

res <- results(dds)
res_no_0 <- res[res$baseMean!=0,]

there_is_padj = !is.na(res_no_0$padj)
res_no_0_with_padj = res_no_0[there_is_padj,]

res_signif = res_no_0_with_padj[res_no_0_with_padj$padj<0.05,]
write.csv(res_signif, 
          file = paste("Projekt/DE_output/", instrument, "_signif_results.csv", sep=""), 
          row.names = TRUE)
