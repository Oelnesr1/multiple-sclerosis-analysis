library(r2r)
library(data.table)
library(fgsea)
library(ggplot2)

DupsAndNorm <- function(df, type=" ") {
  df <- df[!grepl("NULL", df[,1]),]
  
  duplicates <- df[duplicated(df[,1]) | duplicated(df[,1], fromLast=TRUE),]
  df <- df[!(duplicated(df[,1]) | duplicated(df[,1], fromLast=TRUE)),]
  
  duplicates <- duplicates[order(abs(duplicates[,2])),]
  duplicates <- duplicates[!duplicated(duplicates[,1], fromLast=TRUE),]
  
  df <- rbind(df, duplicates)
  rank <- df[order(-df[,2]),]
  
  
  if(type == "rGEO") {rank[,2] <- -log10(rank[,2])}
  
  rank[,2] <- -((rank[,2] - rank[nrow(rank),2]) / (rank[1,2] - rank[nrow(rank),2])) * 10 - 5
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
  }
  
  rank <- as.data.frame(rank[order(rank[,2], decreasing = TRUE),])
  row.names(rank) <- rank[,1]
  return(rank)
}


args = commandArgs(trailingOnly = TRUE)
args[1] = "rGEO"
args[2] = "/Users/omarelnesr/MS-Project/RA_FGSEA/CD8_T_cells"
file_names <- paste(args[2], ".", args[1], sep="")
ensembl_ids <- read.table(paste(file_names, ".tsv", sep=""), sep = "\t", header=F)

chip <- read.delim("/Users/omarelnesr/MS-Project/Human_Ensembl_Gene_ID_MSigDB.v2022.1.Hs.chip", sep = "\t", header=T)

m <- hashmap()
m[chip$Probe.Set.ID] <- chip$Gene.Symbol
symbol <- as.matrix(m[ensembl_ids[,1]])
symbol_df <- data.frame(symbol, ensembl_ids[,2])
rank <- DupsAndNorm(symbol_df, args[1])


openTargets <- read.delim("/Users/omarelnesr/MS-Project/RA_FGSEA/Open_Targets_MS.tsv", sep="\t", header = TRUE)
openTargets <- as.data.frame(openTargets[order(openTargets$overallAssociationScore, decreasing = TRUE),1:2])
openTargetsStatic <- DupsAndNorm(openTargets)

rank_rows <- row.names(rank)
open_rows <- row.names(openTargetsStatic)

open_same <- openTargetsStatic[open_rows %in% rank_rows,]
rank_same <- rank[rank_rows %in% open_rows,]
open_same <- open_same[order(row.names(open_same)),]
rank_same <- rank_same[order(row.names(rank_same)),]

mix <- merge(open_same, rank_same)
colnames(mix) <- c("symbol", "openTargets", "rank")

res <- cor.test(mix$rank, mix$openTargets, method = "pearson")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
res <- cor.test(mix$rank, mix$openTargets, method = "kendall")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
res <- cor.test(mix$rank, mix$openTargets, method = "spearman")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}


mix <- mix[order(-mix$openTargets),]
mix1 <- mix[1:100,]

res <- cor.test(mix1$rank, mix1$openTargets, method = "pearson")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
res <- cor.test(mix1$rank, mix1$openTargets, method = "kendall")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
res <- cor.test(mix1$rank, mix1$openTargets, method = "spearman")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}

mix <- mix[order(-mix$rank),]
mix1 <- mix[1:100,]


res <- cor.test(mix1$rank, mix1$openTargets, method = "pearson")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
res <- cor.test(mix1$rank, mix1$openTargets, method = "kendall")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
res <- cor.test(mix1$rank, mix1$openTargets, method = "spearman")
if (res$estimate >= 0.1 || res$estimate <= -0.1) {print(res$estimate)}
nrow(mix)
