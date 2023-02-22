library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape2)



CD14_positive_monocyte_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/CD14_positive_monocyte.tsv", sep = "\t", header = TRUE, row.names = 1)
CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
CD4_positive_alpha_beta_memory_T_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/CD4_positive_alpha_beta_memory_T_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
CD8_positive_alpha_beta_memory_T_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/CD8_positive_alpha_beta_memory_T_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
IgD_negative_memory_B_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/IgD_negative_memory_B_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
immature_natural_killer_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/immature_natural_killer_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
naive_B_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/naive_B_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
naive_thymus_derived_CD4_positive_alpha_beta_T_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/naive_thymus_derived_CD4_positive_alpha_beta_T_cell.tsv", sep = "\t", header = TRUE, row.names = 1)
naive_thymus_derived_CD8_positive_alpha_beta_T_cell_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_hallmark_sets/naive_thymus_derived_CD8_positive_alpha_beta_T_cell.tsv", sep = "\t", header = TRUE, row.names = 1)

cell_type_array <- list(CD14_positive_monocyte_df,
CD4_positive_alpha_beta_memory_T_cell_df,
CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_df,
CD8_positive_alpha_beta_memory_T_cell_df,
IgD_negative_memory_B_cell_df, immature_natural_killer_cell_df,
naive_B_cell_df, naive_thymus_derived_CD4_positive_alpha_beta_T_cell_df,
naive_thymus_derived_CD8_positive_alpha_beta_T_cell_df)

CD14_positive_monocyte_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/CD14_positive_monocyte_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
CD4_positive_alpha_beta_memory_T_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/CD4_positive_alpha_beta_memory_T_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
CD8_positive_alpha_beta_memory_T_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/CD8_positive_alpha_beta_memory_T_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
IgD_negative_memory_B_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/IgD_negative_memory_B_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
immature_natural_killer_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/immature_natural_killer_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
naive_B_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/naive_B_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)
naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df <- read.table("/Users/omarelnesr/MS-Project/GSEA_c2_sets/naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2.tsv", sep = "\t", header = TRUE, row.names = 1)

for (cell in cell_type_array) {
  cell[order(row.names(cell)), ]
}

for (x in 1:10) {
  needed <- which(rownames(CD14_positive_monocyte_c2_df) %in% rownames(CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df))
  CD14_positive_monocyte_c2_df <- CD14_positive_monocyte_c2_df[needed,]
  CD4_positive_alpha_beta_memory_T_cell_c2_df <- CD4_positive_alpha_beta_memory_T_cell_c2_df[needed,]
  CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df <- CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df[needed,]
  CD8_positive_alpha_beta_memory_T_cell_c2_df <- CD8_positive_alpha_beta_memory_T_cell_c2_df[needed,]
  IgD_negative_memory_B_cell_c2_df <- IgD_negative_memory_B_cell_c2_df[needed,]
  immature_natural_killer_cell_c2_df <- immature_natural_killer_cell_c2_df[needed,]
  naive_B_cell_c2_df <- naive_B_cell_c2_df[needed,]
  naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df <- naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df[needed,]
  naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df <- naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df[needed,]


  CD14_positive_monocyte_c2_df[(rownames(CD14_positive_monocyte_c2_df) %in% rownames(CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df)),]
  CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df[(rownames(CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df) %in% rownames(CD4_positive_alpha_beta_memory_T_cell_c2_df)),]
  CD4_positive_alpha_beta_memory_T_cell_c2_df[(rownames(CD4_positive_alpha_beta_memory_T_cell_c2_df) %in% rownames(CD8_positive_alpha_beta_memory_T_cell_c2_df)),]
  CD8_positive_alpha_beta_memory_T_cell_c2_df[(rownames(CD8_positive_alpha_beta_memory_T_cell_c2_df) %in% rownames(IgD_negative_memory_B_cell_c2_df)),]
  IgD_negative_memory_B_cell_c2_df[(rownames(IgD_negative_memory_B_cell_c2_df) %in% rownames(immature_natural_killer_cell_c2_df)),]
  immature_natural_killer_cell_c2_df[(rownames(immature_natural_killer_cell_c2_df) %in% rownames(naive_B_cell_c2_df)),]
  naive_B_cell_c2_df[(rownames(naive_B_cell_c2_df) %in% rownames(naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df)),]
  naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df[(rownames(naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df) %in% rownames(naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df)),]
  naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df[(rownames(naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df) %in% rownames(CD14_positive_monocyte_c2_df)),]
}

CD14_positive_monocyte_df <- CD14_positive_monocyte_df[order(IgD_negative_memory_B_cell_df$NES), ]
CD4_positive_alpha_beta_memory_T_cell_df <- CD4_positive_alpha_beta_memory_T_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_df <- CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
CD8_positive_alpha_beta_memory_T_cell_df <- CD8_positive_alpha_beta_memory_T_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
naive_thymus_derived_CD4_positive_alpha_beta_T_cell_df <- naive_thymus_derived_CD4_positive_alpha_beta_T_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
IgD_negative_memory_B_cell_df <- IgD_negative_memory_B_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
naive_B_cell_df <- naive_B_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
naive_thymus_derived_CD8_positive_alpha_beta_T_cell_df <- naive_thymus_derived_CD8_positive_alpha_beta_T_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]
immature_natural_killer_cell_df <- immature_natural_killer_cell_df[order(IgD_negative_memory_B_cell_df$NES), ]

CD14_positive_monocyte_c2_df <- CD14_positive_monocyte_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
CD4_positive_alpha_beta_memory_T_cell_c2_df <- CD4_positive_alpha_beta_memory_T_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df <- CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
CD8_positive_alpha_beta_memory_T_cell_c2_df <- CD8_positive_alpha_beta_memory_T_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df <- naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
IgD_negative_memory_B_cell_c2_df <- IgD_negative_memory_B_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
naive_B_cell_c2_df <- naive_B_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df <- naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]
immature_natural_killer_cell_c2_df <- immature_natural_killer_cell_c2_df[order(IgD_negative_memory_B_cell_c2_df$NES), ]

head(CD14_positive_monocyte_df)
print(CD14_positive_monocyte_df[3])

NES_h_df <- cbind(CD14_positive_monocyte_df[3], CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_df[3], CD4_positive_alpha_beta_memory_T_cell_df[3], CD8_positive_alpha_beta_memory_T_cell_df[3], IgD_negative_memory_B_cell_df[3], immature_natural_killer_cell_df[3], naive_B_cell_df[3], naive_thymus_derived_CD4_positive_alpha_beta_T_cell_df[3], naive_thymus_derived_CD8_positive_alpha_beta_T_cell_df[3])
names(NES_h_df) <- c("CD14_positive_monocyte", "CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell", "CD4_positive_alpha_beta_memory_T_cell", "CD8_positive_alpha_beta_memory_T_cell_df", "IgD_negative_memory_B_cell", "immature_natural_killer_cell", "naive_B_cell", "naive_thymus_derived_CD4_positive_alpha_beta_T_cell", "naive_thymus_derived_CD8_positive_alpha_beta_T_cell")
head(NES_h_df)



# head(CD14_positive_monocyte_c2_df)

NES_c2_df <- cbind(CD14_positive_monocyte_c2_df[3], CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell_c2_df[3], CD4_positive_alpha_beta_memory_T_cell_c2_df[3], CD8_positive_alpha_beta_memory_T_cell_c2_df[3], IgD_negative_memory_B_cell_c2_df[3], immature_natural_killer_cell_c2_df[3], naive_B_cell_c2_df[3], naive_thymus_derived_CD4_positive_alpha_beta_T_cell_c2_df[3], naive_thymus_derived_CD8_positive_alpha_beta_T_cell_c2_df[3])
names(NES_c2_df) <- c("CD14_positive_monocyte", "CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell", "CD4_positive_alpha_beta_memory_T_cell", "CD8_positive_alpha_beta_memory_T_cell", "IgD_negative_memory_B_cell", "immature_natural_killer_cell", "naive_B_cell", "naive_thymus_derived_CD4_positive_alpha_beta_T_cell", "naive_thymus_derived_CD8_positive_alpha_beta_T_cell")
# print(NES_c2_df)

# heatmap.2(data.matrix(NES_h_df), trace = "none", dendrogram = "none", Rowv = "none", Colv = "none")

head(NES_h_df)

NES_h_data <- cor(NES_h_df[sapply(NES_h_df, is.numeric)])
NES_h_data1 <- melt(data.matrix(NES_h_df), value.name = "NES")
# , value.name = "NES")

# head(NES_h_data)
head(NES_h_data1)

# head(data.matrix(NES_h_df))



NES_c2_data <- cor(NES_c2_df[sapply(NES_c2_df, is.numeric)])

# NES_c2_df <- NES_c2_df[order("CD14_positive_monocyte"),]
# head(NES_c2_df)

NES_c2_data1 <- melt(data.matrix(NES_c2_df), value.name = "NES")


head(NES_c2_data1)

png("Hallmark_Heatmap.png", res = 1000, units = "in", width = 10, height = 10)
ggplot(NES_h_data1, aes(x = Var2, y = Var1, fill = NES)) +
geom_tile() +
scale_fill_distiller(palette = "Spectral") +
theme_gray(base_size =  8) +
# scale_x_discrete(expand = c(0,0), breaks = c("Monocyte", "reg T", "CD4 mem T", "CD8 mem T", "IgD- mem B", "ink", "naive B", "naive CD4 T", "naive CD8 T")) +
guides(fill = guide_colourbar(barwidth = 3, barheight = 25), ) +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
    # axis.text.x = c("Monocyte", "reg T", "CD4 mem T", "CD8 mem T", "IgD- mem B", "ink", "naive B", "naive CD4 T", "naive CD8 T"),
    axis.title = element_text(size = (20))) +
labs(y = "Gene Sets", x = "Cell Type")
# ggplot(NES_h_data1, aes(x = Var2, y = Var1, fill = NES)) + geom_tile() + scale_fill_distiller(palette = "Spectral") + theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = (20))) + guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) + labs(y = "Gene Sets", x = "Cell Types")
dev.off()

png("C2_Heatmap.png", res = 1000, units = "in", width = 10, height = 10)
ggplot(NES_c2_data1, aes(x = Var2, y = Var1, fill = NES)) +
geom_tile() +
scale_fill_distiller(palette = "Spectral") +
theme_gray(base_size =  8) +
# scale_x_discrete(expand = c(0,0), breaks = c("Monocyte", "reg T", "CD4 mem T", "CD8 mem T", "IgD- mem B", "ink", "naive B", "naive CD4 T", "naive CD8 T")) +
guides(fill = guide_colourbar(barwidth = 3, barheight = 25), ) +
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
    # axis.text.x = c("Monocyte", "reg T", "CD4 mem T", "CD8 mem T", "IgD- mem B", "ink", "naive B", "naive CD4 T", "naive CD8 T"),
    axis.title = element_text(size = (20))) +
labs(y = "Gene Sets", x = "Cell Type")
dev.off()

# coul <- colorRampPalette(brewer.pal(8, "Spectral"))(25)
# heatmap.2(data.matrix(NES_c2_df), trace = "none", dendrogram = "none", Rowv = "none", Colv = "none", col = coul)
