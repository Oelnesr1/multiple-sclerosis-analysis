library(ggplot2)
library(umap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

count_file_name = args[1]
new_file_name = args[2]

count_file_name = "/Users/omarelnesr/MS-Project/Datasets/ENCODE/Count-Matrices/all-Cell.Count-Matrix.txt"
new_file_name = "/Users/omarelnesr/MS-Project/Datasets/ENCODE/cell_types/T_Cells/CD8_T_cells"

count_matrix <- read.table(count_file_name, sep="\t", row.names = "Gene", header = T)
count_matrix <- count_matrix[!grepl("PAR", row.names(count_matrix)),]
count_matrix <- count_matrix[!grepl("__", row.names(count_matrix)),]
rownames(count_matrix) = gsub("\\..*", "", row.names(count_matrix))
count_matrix

condition = rep(c(rep("Healthy",3), rep("MS",3)),9)
cell_types <- c(rep("CD14+ Monocyte", 6), rep("CD4+ Memory T cell", 6), rep("CD4+ Regulatory T cell",6), rep("CD8+ Memory T cell",6), rep("IgD- Memory B cell",6), 
                rep("Immature NK cell",6), rep("Naive B cell",6), rep("Naive Thymus CD4+ T cell",6), rep("Naive Thymus CD8+ T cell",6))

broader_cell_types <- c(rep("CD14+ Monocyte", 6), rep("T cell", 18), rep("B cell", 6), rep("NK cell", 6), rep("B cell", 6), rep("T cell", 12))
T_cells <- c(rep("CD14+ Monocyte", 6), rep("CD4+ T cell", 12), rep("CD8+ T cell", 6), rep("IgD- Memory B cell",6), rep("NK cell", 6), rep("Naive B cell", 6), rep("CD4+ T cell", 6), rep("CD8+ T cell", 6))
donor = names(count_matrix)
metadata = data.frame(donor, condition, cell_types, broader_cell_types)
row.names(metadata) = names(count_matrix)
metadata$condition1 <- paste(condition, cell_types)
metadata$condition2 <- paste(condition, broader_cell_types)
metadata$condition3 <- paste(condition, T_cells)

count_matrix <- count_matrix[,row.names(metadata)]

dds_before = DESeqDataSetFromMatrix(countData = round(count_matrix), colData = metadata, design = ~ condition3)
# dds_middle <- collapseReplicates(dds_before, metadata$Sample.Name)
dds_after = DESeq(dds_before, betaPrior = TRUE)
res = results(dds_after)
res <- results(dds_after, contrast=c("condition1","MS CD4+ Regulatory T cell","Healthy CD4+ Regulatory T cell"))
res <- results(dds_after, contrast=c("condition1","MS CD8+ Memory T cell","Healthy CD8+ Memory T cell"))
res <- results(dds_after, contrast=c("condition1","MS IgD- Memory B cell","Healthy IgD- Memory B cell"))
res <- results(dds_after, contrast=c("condition1","MS Immature NK cell","Healthy Immature NK cell"))
res <- results(dds_after, contrast=c("condition1","MS Naive B cell","Healthy Naive B cell"))
res <- results(dds_after, contrast=c("condition1","MS Naive Thymus CD4+ T cell","Healthy Naive Thymus CD4+ T cell"))
res <- results(dds_after, contrast=c("condition1","MS Naive Thymus CD8+ T cell","Healthy Naive Thymus CD8+ T cell"))
res <- results(dds_after, contrast=c("condition3","MS CD4+ T cell","Healthy CD4+ T cell"))
res <- results(dds_after, contrast=c("condition3","MS CD8+ T cell","Healthy CD8+ T cell"))

orderedResultspadj = res[order(res$padj),]
orderedResultsStat = res[order(-res$stat),]
stat_df = data.frame(orderedResultsStat)

GSEA <- as.data.frame(stat_df[4][!is.na(stat_df$stat),])
row.names(GSEA) <- row.names(stat_df)[!is.na(stat_df$stat)]
colnames(GSEA) <- c("")

write.table(res, paste(new_file_name, ".DESeq2_Results.tsv", sep=""), sep = "\t", quote=F)
write.table(orderedResultspadj, paste(new_file_name, ".DESeq2_Results.padj.tsv", sep=""), sep = "\t", quote=F)
write.table(GSEA, paste(new_file_name, ".DESeq2_Results.gsea.tsv", sep=""), sep = "\t", quote=F)

results = data.frame(res)

png(paste(new_file_name, ".DESeq2.svg", sep = ""), width = 10, height = 10)
ggplot(data = results, aes(x = log2FoldChange, y = -log10(padj))) + geom_point() +
        xlab("Log 2 Fold Change") + ylab("Log 10 Statistical Significance") +
        ggtitle("Volcano Plot of Differential Gene Expression", subtitle = "CD8+ T cells in MS") +
theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20))
dev.off()

upReg = rownames(results)[results$log2FoldChange >= 1 & results$padj < 0.05]
upReg = upReg[!is.na(upReg)]
write.table(upReg, paste(new_file_name, ".DESeq2_Up.tsv", sep=""), sep="\t", quote=F, row.names=F, col.names=F)

downReg = rownames(results)[results$log2FoldChange <= -1 & results$padj < 0.05]
downReg = downReg[!is.na(downReg)]
write.table(downReg, paste(new_file_name, ".DESeq2_Down.tsv", sep=""), sep="\t", quote=F, row.names=F, col.names=F)

