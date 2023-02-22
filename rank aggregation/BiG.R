library("RobustRankAggreg")
library("data.table")
source("/Users/omarelnesr/MS-Project/R-Scripts/BiG_code_platform_changed.R")
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

args = commandArgs(trailingOnly = TRUE)
DupsAndNorm <- function(df, type="rGEO") {
  df <- df[!grepl("NULL", df[,1]),]
  
  duplicates <- df[duplicated(df[,1]) | duplicated(df[,1], fromLast=TRUE),]
  df <- df[!(duplicated(df[,1]) | duplicated(df[,1], fromLast=TRUE)),]
  
  duplicates <- duplicates[order(abs(duplicates[,2])),]
  duplicates <- duplicates[!duplicated(duplicates[,1], fromLast=TRUE),]
  
  df <- rbind(df, duplicates)
  rank <- df[order(-df[,2]),]
  
  
  if(type == "rGEO") {rank[,2] <- -log10(rank[,2])}
  
  rank[,2] <- -(((rank[,2] - rank[nrow(rank),2]) / (rank[1,2] - rank[nrow(rank),2])) * 20 - 10)
  rank <- rank[order(rank[,2], decreasing = FALSE),]
  
  duplicate_ranks <- rank[duplicated(rank[,2]) | duplicated(rank[,2], fromLast=TRUE),]
  while(nrow(duplicate_ranks) != 0) {
    duplicate_ranks <- rank[duplicated(rank[,2]) | duplicated(rank[,2], fromLast=TRUE),]
    rank <- rank[!(duplicated(rank[,2]) | duplicated(rank[,2], fromLast=TRUE)),]
    duplicates_to_change <- duplicate_ranks[duplicated(duplicate_ranks[,2]),]
    duplicate_ranks <- subset(duplicate_ranks, !(symbol %in% duplicates_to_change$symbol))
    duplicates_to_change[,2] <- duplicates_to_change[,2] - 0.000000001
    rank <- rbind(rank, duplicate_ranks, duplicates_to_change)
    duplicate_ranks <- rank[duplicated(rank[,2]) | duplicated(rank[,2], fromLast=TRUE),]
    print(nrow(duplicate_ranks))
  }
  rank <- as.data.frame(rank[order(rank[,2], decreasing = TRUE),])
  row.names(rank) <- rank[,1]
  return(rank)
}
args[1] <- 3
args[2] <- "blobby"
args[3] <- "/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation_hg38_gencode.DESeq2-Results.gsea.tsv"
args[4] <- "/Users/omarelnesr/MS-Project/Datasets/ENCODE/cell_types/B_Cells/B_cells.DESeq2_Results.gsea.tsv"
args[5] <- "/Users/omarelnesr/MS-Project/Datasets/Linsley/Linsley_B_cells.DESeq2_Results.gsea.tsv"

my_file_list <- args[3:(2+as.numeric(args[1]))]
y <- lapply(my_file_list, read.table)
y <- lapply(y, function(df){df[order(df$V2, decreasing=TRUE),1]})
r <- matrix_transfer(y, full=FALSE, mean = TRUE)

NTlength <- rep(1000,length(y))
for (i in 1:length(y)) {
  NTlength[i] <- length(y[[i]])
}

print("Beginning rGEO Rank Aggregation...")
res <- aggregateRanks(y, method="geom.mean")

names(res) <- c("symbol", "score")
res <- DupsAndNorm(res)

y <- append(y, list(res[,1]))
r <- matrix_transfer(y, full=FALSE, mean = TRUE)
r <- r[order(r[,ncol(r)]),]

heat <- melt(r)
png("B_cell_Heatmap.png", res = 1000, units = "in", width = 10, height = 10)
ggplot(heat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  theme_gray(base_size =  8) +
  guides(fill = guide_colourbar(barwidth = 3, barheight = 25), ) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.title = element_text(size = (20))) +
  labs(y = "Placement in List", x = "Gene Ranks")
dev.off()


fwrite(res, paste(args[2], ".rGEO.tsv", sep=""), sep="\t", quote=F, row.names=F, col.names=F)
print("rGEO Rank Aggregation Complete")

print("Beginning BiG Rank Aggregation...")
result <- BiG_diffuse(r=r, n_T=NTlength, n_p1=0, M=2000, burnin=1000, prior="IG")
print("BiG Rank Aggregation Complete")

rankedEntities <- cbind(row.names(r), result)
rankedEntities <- rankedEntities[order(result),]
fwrite(rankedEntities, paste(args[2], ".BiG.tsv", sep=""), sep="\t", quote=F, row.names=F, col.names=F)


