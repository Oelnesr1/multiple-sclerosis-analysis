library(ggplot2)
library(umap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)

count_file_name <- "/Users/omarelnesr/MS-Project/Datasets/Linsley/Linsley.Count_Matrix.tsv"
new_file_name <- "/Users/omarelnesr/MS-Project/Datasets/Linsley/Linsley"

count_matrix <- read.table(count_file_name, sep="\t", row.names = "Genes", header = T)
count_matrix <- count_matrix[!grepl("PAR", row.names(count_matrix)),]
count_matrix <- count_matrix[!grepl("__", row.names(count_matrix)),]
row.names(count_matrix) = gsub("\\..*", "", row.names(count_matrix))

metadata_file <- read.csv("/Users/omarelnesr/MS-Project/Datasets/Linsley/Linsley_Metadata.csv", header = T)
metadata <- metadata_file[order(metadata_file$Run),c(1, 10,17,18,20,29,33)]
metadata <- rbind(metadata[(metadata$diseasestatus == "MS pretreatment"),],metadata[(metadata$diseasestatus == "Healthy Control"),])
metadata <- metadata[order(metadata$Run),]
metadata$condition <- paste(metadata$celltype, metadata$diseasestatus)
row.names(metadata) <- metadata$Run
metadata <- metadata[1:(nrow(metadata)-1),2:ncol(metadata)]
metadata

count_matrix <- count_matrix[,row.names(metadata)]

interest <- c(1, 10, 18, 25, 27, 29, 31, 33)
metada <- metadata_file[(metadata_file$Sequencing_method != "3' single cell RNAseq 10x Genomics v2"),interest]
colnames(metadata)

diseaseKey <- rep(c(rep("PPMS", 4), rep("SPMS", 4), rep("RRMS", 3), rep("Healthy", 7)))

pca <- prcomp(transposed_tempData)

pcs <- data.frame(pca$x)

png(paste(new_file_name, "_PCA.png"), res = 500, units = "in", width = 10, height = 10)
ggplot(data = pcs, aes(x = PC1, y = PC2, color=diseaseKey )) + geom_point()
dev.off()


# tempData=read.table(args[1], header=TRUE, sep="\t", row.names=1)
transposed_tempData=t(count_matrix)

umapOutput=umap(transposed_tempData, n_components = 2, n_neighbors = 20)
umapOutput=data.frame(umapOutput$layout)

# png(paste(new_file_name, "_UMAP_NoColors.png"), res=500, units="in", width=10, height=10)
# ggplot(data=umapOutput, aes(x=X1, y=X2))+geom_point()
# dev.off()

png(paste(new_file_name, "_UMAP_GroupColors.png"), res=500, units="in", width=10, height=10)
ggplot(data=umapOutput, aes(x=X1, y=X2, color=metadata$condition))+geom_point()
dev.off()

png(paste(new_file_name, "_UMAP_SampleColors.png"), res=500, units="in", width=10, height=10)
ggplot(data=umapOutput, aes(x=X1, y=X2, color=metadata$Sample.Name))+geom_point()
dev.off()

png(paste(new_file_name, "_UMAP_GenderColors.png"), res=500, units="in", width=10, height=10)
ggplot(data=umapOutput, aes(x=X1, y=X2, color=metadata$gender))+geom_point()
dev.off()
