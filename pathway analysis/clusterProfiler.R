library(clusterProfiler)
library(DOSE)
library(ChIPseeker)
library(enrichplot)
library(meshes)
library(GOSemSim)
library(ReactomePA)
library(dplyr)
library(msigdbr)
library(data.table)
library(r2r)
library(fgsea)
rm(list = ls())

DupsAndNorm <- function(df, type=" ") {
  df <- df[!grepl("NULL", df[,1]),]
  
  duplicates <- df[duplicated(df[,1]) | duplicated(df[,1], fromLast=TRUE),]
  df <- df[!(duplicated(df[,1]) | duplicated(df[,1], fromLast=TRUE)),]
  
  duplicates <- duplicates[order(abs(duplicates[,2])),]
  duplicates <- duplicates[!duplicated(duplicates[,1], fromLast=TRUE),]
  
  df <- rbind(df, duplicates)
  rank <- df[order(-df[,2]),]
  
  
  if(type == "rGEO") {rank[,2] <- -log10(rank[,2])}
  
  rank[,2] <- ((rank[,2] - rank[nrow(rank),2]) / (rank[1,2] - rank[nrow(rank),2])) * 20 - 10
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
fix_list <- function(symbol_df) {
  res <- list(NULL)
  for (i in symbol_df) {
    symbol <- as.vector(as.character(i[,1]))
    score <- as.vector(as.numeric(i[,2]))
    df <- list(data.frame(symbol, score))
    res <- append(res, df)
  }
  return(res[-1])
}
GSEA_cp <- function(gl, collapsedPathways) {
  res <- list(NULL)
  for (i in 1:length(gl)) {
    print(i)
    gsea_res <- GSEA(gl[[i]], TERM2GENE = collapsedPathways[[i]], nPermSimple = 1000, eps = 0, minGSSize = 1, maxGSSize = 10000)
    res <- append(res, gsea_res)
  }
  return(res[-1])
}
convert_hash <- function(hashmap, object) {
  return(hashmap[object[,1]])
}
make_cp <- function(rank) {
  gl <- rank[,2]
  names(gl) <- as.character(rank[,1])
  gl <- sort(gl, decreasing = TRUE)
  return(gl)
}
make_fg <- function(rank, gene_sets_fg, gene_sets_cp) {
  print("making ranks")
  ranks <- setNames(rank[,2], rank[,1])
  
  print("running fgsea")
  fgseaRes <- fgsea(gene_sets_fg, ranks, nPermSimple = 5000)
  
  print("collapsing pathways")
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], gene_sets_fg, ranks)
  
  print("making pathways for clusterProfiler")
  collapsedPathways <- gene_sets_cp[gene_sets_cp$term %in% collapsedPathways$mainPathways,]
  return(collapsedPathways)
}

args = commandArgs(trailingOnly = TRUE)

RA <- "rGEO"
cell_types <- list("B_cells", "CD4_T_cells", "CD8_T_cells", "Monocytes", "NAWM", "NK_cells")
base_dir <- "/Users/omarelnesr/MS-Project/RA_FGSEA"

files <- lapply(function(X) paste(base_dir, "/", X, ".", RA, ".tsv", sep=""), X = cell_types)
ensembl_ids <- lapply(files, read.delim, sep = "\t", header = F)

chip <- read.delim("/Users/omarelnesr/MS-Project/Human_Ensembl_Gene_ID_MSigDB.v2022.1.Hs.chip", sep = "\t", header=T)
gmt.file <- "/Users/omarelnesr/MS-Project/msigdb.v2022.1.Hs.symbols.gmt"
gene_sets_cp <- read.gmt(gmt.file)
gene_sets_fg <- gmtPathways(gmt.file)

m <- hashmap()
m[chip$Probe.Set.ID] <- chip$Gene.Symbol

symbol <- lapply(lapply(ensembl_ids, convert_hash, hashmap = m), as.matrix)

symbol_df <- lapply(mapply(function(X,Y) cbind(X, as.numeric(Y[,2])), X = symbol, Y = ensembl_ids), data.frame)

symbol_df <- fix_list(symbol_df)

rank_symbol <- lapply(symbol_df, DupsAndNorm, type=RA)
gl <- lapply(rank_symbol, make_cp)

symbol_entrez <- lapply(X = rank_symbol, function(X) mapIds(org.Hs.eg.db, keys = X[,1], keytype = "GENENAME", 
                                                     column = "ENTREZID", multiVals = "first"))

collapsedPathways <- lapply(rank, make_fg, gene_sets_fg = gene_sets_fg, gene_sets_cp = gene_sets_cp)

hmm <- GSEA_cp(gl, collapsedPathways)

for (i in 1:length(hmm)) {
  msig <- arrange(hmm[[1]], abs(p.adjust))
  png(paste(cell_types[[1]],".dotplot.png", sep = ""), res=1000, units="in", width=6, height=10)
  enrichplot::dotplot(msig, showCategory = 25 ,x = "NES", size = "GeneRatio", font.size = 6.5 )
  dev.off()
}
