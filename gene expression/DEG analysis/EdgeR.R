library("edgeR")
library("ggplot2")

args = commandArgs(trailingOnly = TRUE)

count_file_name <- args[1]
new_file_name <- args[2]

count_file_name <- "/Users/omarelnesr/MS-Project/B_cell_hypomethylation_count_matrix.tsv"
new_file_name <- "B_cell_hypomethylation"

count_file_name <- "/Users/omarelnesr/MS-Project/Count_Matrices/B_cell_hypomethylation_hg38_gencode.Count_Matrix.tsv"
new_file_name <- "B_cell_hypomethylation_hg38_gencode"

count_matrix <- read.table(count_file_name, sep="\t", row.names = "Genes", header = T)
count_matrix <- count_matrix[!grepl("PAR", row.names(count_matrix)),]
count_matrix <- count_matrix[!grepl("__", row.names(count_matrix)),]
row.names(count_matrix) <- gsub("\\..*", "", row.names(count_matrix))

disease_key = c(rep("multiple sclerosis", 23), rep("healthy", 14))
group = factor(disease_key)
dge <- DGEList(counts = count_matrix, group=group)
keep <- filterByExpr(y = dge)
y = dge[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y)
fit <- glmQLFit(y, design)
tr <- glmTreat(fit, lfc=0)
qlf <- glmQLFTest(fit)
et <- exactTest(object = y)

top_degs = topTags(object = qlf, n = "Inf")
results = data.frame(top_degs)

png(paste(new_file_name, ".EdgeR_PVal.png", sep = ""), res = 1000, units = "in", width = 10, height = 10)
ggplot(data = results, aes(x = logFC, y = -log10(PValue))) + geom_point() +
  xlab("Log 2 Fold Change") + ylab("Log 10 Statistical Significance") +
  ggtitle("Volcano Plot of Differential Gene Expression", subtitle = "CD19+ B-Cells in MS") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20))
dev.off()

png(paste(new_file_name, ".EdgeR_FDR.png", sep = ""), res = 1000, units = "in", width = 10, height = 10)
ggplot(data = results, aes(x = logFC, y = -log10(FDR))) + geom_point() +
  xlab("Log 2 Fold Change") + ylab("Log 10 False Discovery Rate") +
  ggtitle("Volcano Plot of Differential Gene Expression", subtitle = "CD19+ B-Cells in MS") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20))
dev.off()


write.table(results, sep = "\t", file = paste(new_file_name, ".EdgeR.tsv", sep = ""), quote = FALSE)

upReg = Genes[results$logFC >= 1 & results$FDR <= 0.05]
write.table(upReg, paste(new_file_name, ".EdgeR_Up.tsv", sep=""), sep="\t", quote=F, row.names=F, col.names=F)

downReg = Genes[results$logFC <= -1 & results$FDR <= 0.05]
write.table(downReg, paste(new_file_name, ".EdgeR_Down.tsv", sep=""), sep="\t", quote=F, row.names=F, col.names=F)

