
Libraries <- function() {
  library(GEOquery)
  library(BiocGenerics)
  library(limma)
  library(umap)
  library(ggplot2)
  library(RColorBrewer)
  library(oligo)
  library(arrayQualityMetrics)
  library(pd.hg.u133a.2)
  library(hgu133plus2.db)
  library(hgu133a2.db)
  library(AnnotationDbi)
  library(gplots)
  library(geneplotter)
  library(pheatmap)
  library(enrichplot)
  library(topGO)
  library(ReactomePA)
  library(clusterProfiler)
  library(BiocStyle)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(matrixStats)
  library(genefilter)
  library(openxlsx)
  library(ArrayExpress)
  library(pd.hugene.1.0.st.v1)
  library(hugene10sttranscriptcluster.db)
  library(Biobase)
  library(oligoClasses)
  library(maEndToEnd)
  library(rmarkdown)
  library(devtools)
  library(r2r)
  library(data.table)
  library(fgsea)
  library(RCurl)
  library(ACME)
  library(genefilter)
  library(maqcExpression4plex)
  library(knitr)
  library(RUnit)
}

Libraries()

celFiles <- oligoClasses::list.celfiles(path = "~/MS-Project/GSE62584_RAW", full.name = TRUE, listGzipped = TRUE)
affyRaw <- oligo::read.celfiles(celFiles)
eset <- oligo::rma(affyRaw)

gset <- getGEO("GSE62584", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
row.names(pData(gset)) <- colnames(exprs(eset))
pData(eset) <- pData(gset)
fData(eset) <- fData(gset)

pd <- pData(eset)[,c(1,36:39)]
names(pd) <- c("individual", "age", "cell_type", "disease_state", "gender")
row.names(pd) <- row.names(pData(eset))
pData(eset) <- pd




df <- AnnotationDbi::select(hgu133a2.db, keys=row.names(eset), columns=c("ENSEMBL", "SYMBOL"), keytype="PROBEID")
df <- df[!duplicated(df$PROBEID),]

individual <- row.names(pData(eset))

cell_type <- c(rep("CD4", 15), rep("CD19", 15), rep("CD14", 15), rep("CD8", 15))

disease_state <- c(rep(c(rep("HC", 9), rep("MS", 6)), 4))

design <- model.matrix(~ 0 + disease_state[1:15])
row.names(design) <- individual[cell_type == "CD4"]
colnames(design) <- c("HC", "MS")

cont.mat <- makeContrasts(MS-HC, levels = design)

CD4_fit <- eBayes(contrasts.fit(lmFit(eset[,cell_type == "CD4"], design = design), cont.mat))
row.names(design) <- individual[cell_type == "CD8"]
CD8_fit <- eBayes(contrasts.fit(lmFit(eset[,cell_type == "CD8"], design = design), cont.mat))
row.names(design) <- individual[cell_type == "CD14"]
CD14_fit <- eBayes(contrasts.fit(lmFit(eset[,cell_type == "CD14"], design = design), cont.mat))
row.names(design) <- individual[cell_type == "CD19"]
CD19_fit <- eBayes(contrasts.fit(lmFit(eset[,cell_type == "CD19"], design = design), cont.mat))

CD4 <- topTable(CD4_fit,  adjust.method = "BH",
                sort.by="t", number=Inf)
CD4 <- subset(CD4, select=c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
CD4 <- CD4[order(CD4$adj.P.Val),]



CD8 <- topTable(CD8_fit,  adjust.method = "BH",
                sort.by="t", number=Inf)
CD8 <- subset(CD8, select=c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
CD8[order(CD8$adj.P.Val),]
  

CD14 <- topTable(CD14_fit,  adjust.method = "BH",
                sort.by="t", number=Inf)
CD14 <- subset(CD14, select=c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
CD14[order(CD14$t),]

CD19 <- topTable(CD19_fit,  adjust.method = "BH",
                 sort.by="t", number=Inf)
CD19 <- subset(CD19, select=c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
CD19[order(CD19$t),]

hist(CD4$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
hist(CD8$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
hist(CD14$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
hist(CD19$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

CD4_fit$coefficients
CD4_fit$genes$Gene.symbol
volcano_names <- ifelse(abs(CD4_fit$coefficients)>=1, 
                        CD4_fit$genes$Gene.symbol, NA)
volcanoplot(CD4_fit, coef = 1L, style = "p-value", highlight = 100, 
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

volcanoplot(palmieri_fit_CD, coef = 1L, style = "p-value", highlight = 100, 
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

df <- AnnotationDbi::select(hgu133a2.db, keys=row.names(CD4), columns=c("ENSEMBL", "SYMBOL"), keytype="PROBEID")
nrow(df)
df <- df[!duplicated(df$PROBEID),]
nrow(df)
df1 <- as.data.frame(CD4$t)
row.names(df1) <- row.names(CD4)
head(df1)
head(CD4[,3])
df <- merge(df, CD4$t)
df
df <- df[!duplicated(df$ENSEMBL),]
nrow(df)
df
x


CD4_fit
?eBayes

?model.matrix

disease_state
?select

medians <- rowMedians(Biobase::exprs(eset))

man_threshold <- 4

hist_res <- hist(medians, 100, col = "cornsilk", freq = FALSE, 
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)


pData(eset) <- pd
eset




row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(eset)))

RLE_data <- sweep(Biobase::exprs(eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4", 
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))
