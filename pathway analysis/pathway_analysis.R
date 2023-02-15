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
library(org.Hs.eg.db)
library(AnnotationDbi)
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
convert_ensembl_entrez <- function(ensembl_ids, gl_entrez) {
  res = list(NULL)
  for (i in 1:length(gl_entrez)) {
    score <- as.vector(as.numeric(ensembl_ids[[i]][,2]))
    symbol <- as.vector(as.character(gl_entrez[[i]]))
    df <- data.frame(symbol, score)
    df <- list(df[!is.na(df$symbol),])
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

symbol_df <- fix_list(ensembl_ids)

rank <- lapply(symbol_df, DupsAndNorm, type=RA)
gl_ensembl <- lapply(rank, make_cp)

entrez_convert <- lapply(X = ensembl_ids, function(X) mapIds(org.Hs.eg.db, keys = X[,1], keytype = "ENSEMBL",
                                                             column = "ENTREZID", multiVals = "first"))

entrez_df <- convert_ensembl_entrez(ensembl_ids, entrez_convert)
rank_entrez <- lapply(entrez_df, DupsAndNorm, type=RA)
gl_entrez <- lapply(rank_entrez, make_cp)

# GO_res <- lapply(X = gl_ensembl, function(X) gseGO(geneList = X, OrgDb = org.Hs.eg.db,
#                                            keyType = 'ENSEMBL', minGSSize = 1, maxGSSize = 10000,
#                                            eps = 0.0, ont = "ALL", pAdjustMethod = "BH",
#                                            pvalueCutoff  = 0.01))

# GO_res <- lapply(GO_res, simplify)

# DO_res <- lapply(X = gl_entrez, function(X) gseDO(geneList = X, minGSSize = 1, maxGSSize = 10000,
#                                                    eps = 0.0,
#                                                    pvalueCutoff  = 0.01))


# DGN_res <- lapply(X = gl_entrez, function(X) gseDGN(geneList = X, minGSSize = 1, maxGSSize = 10000,
#                                                    eps = 0.0,
#                                                    pvalueCutoff  = 0.01))

# KEGG_res <- lapply(X = gl_entrez, function(X) gseKEGG(geneList = X, organism = "hsa", minGSSize = 1,
#                                                       maxGSSize = 10000, eps = 0.0, pvalueCutoff  = 0.01))
# MKEGG_res <- lapply(X = gl_entrez, function(X) gseMKEGG(geneList = X, organism = "hsa", minGSSize = 1,
#                                                       maxGSSize = 10000, eps = 0.0, pvalueCutoff  = 0.01))
# REACTOME_res <- lapply(X = gl_entrez, function(X) gsePathway(geneList = X, minGSSize = 1,
#                                                              maxGSSize = 10000, eps = 0.0, pvalueCutoff  = 0.01))

# WP_res <- lapply(X = gl_entrez, function(X) gseWP(geneList = X, minGSSize = 1, organism = "Homo sapiens",
#                                                        maxGSSize = 10000, eps = 0.0, pvalueCutoff  = 0.01))

WPPrint <- as.data.frame(WP_res[[6]])[order(-WP_res[[6]]$NES),1:10]
fwrite(WPPrint, file = paste(cell_types[[6]], ".WPRes.txt", sep = ""), sep="\t", sep2=c("", " ", ""))


msig <- arrange(WP_res[[1]], abs(p.adjust))
png(paste(cell_types[[1]],".Reactome.dotplot.png", sep = ""), res=1000, units="in", width=5, height=8)
enrichplot::dotplot(msig, showCategory = 15 ,x = "NES", size = "GeneRatio", font.size = 9)
dev.off()


msig <- arrange(REACTOME_res[[5]], abs(p.adjust))
svg(paste(cell_types[[5]],".REACTOME.dotplot.svg", sep = ""), width=5, height=8)
enrichplot::dotplot(msig, showCategory = 15 ,x = "NES", size = "GeneRatio", font.size = 9)
dev.off()
