# I did not make this file! I copied it from this paper and made minor edits: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6063319/
library("maEndToEnd")
library(Biobase)
library(oligoClasses)
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(oligo)
library(arrayQualityMetrics)
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)
library(rmarkdown)
library(BiocStyle)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)

raw_data_dir <- "/Users/omarelnesr/MS-Project/E-MTAB-2967"
if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}

sdrf_location <- file.path(raw_data_dir, "E-MTAB-2967.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Source.Name",
                                                         "Characteristics.individual.",
                                                         "Factor.Value.disease.",
                                                         "Factor.Value.phenotype.")]
exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = pData(raw_data)$Factor.Value.disease.,
                     Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                     Individual = pData(raw_data)$Characteristics.individual.)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
oligo::boxplot(raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

palmieri_eset <- oligo::rma(raw_data, target = "core", normalize = FALSE)
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))

RLE_data <- sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

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
palmieri_eset_norm <- oligo::rma(raw_data, target = "core")
exp_palmieri <- Biobase::exprs(palmieri_eset_norm)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease = 
                       Biobase::pData(palmieri_eset_norm)$Factor.Value.disease.,
                     Phenotype = 
                       Biobase::pData(palmieri_eset_norm)$Factor.Value.phenotype.)


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
phenotype_names <- ifelse(str_detect(pData
                                     (palmieri_eset_norm)$Factor.Value.phenotype.,
                                     "non"), "non_infl.", "infl.")

disease_names <- ifelse(str_detect(pData
                                   (palmieri_eset_norm)$Factor.Value.disease.,
                                   "Crohn"), "CD", "UC")

annotation_for_heatmap <- 
  data.frame(Phenotype = phenotype_names,  Disease = disease_names)

row.names(annotation_for_heatmap) <- row.names(pData(palmieri_eset_norm))
dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))

rownames(dists) <- row.names(pData(palmieri_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  Phenotype = c(non_infl. = "chartreuse4", infl. = "burlywood3"),
  Disease = c(CD = "blue4", UC = "cadetblue2")
)
pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")
palmieri_medians <- rowMedians(Biobase::exprs(palmieri_eset_norm))

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
man_threshold <- 4

hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE, 
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)
no_of_samples <- 
  table(paste0(pData(palmieri_eset_norm)$Factor.Value.disease., "_", 
               pData(palmieri_eset_norm)$Factor.Value.phenotype.))
no_of_samples 
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(palmieri_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
palmieri_manfiltered <- subset(palmieri_eset_norm, idx_man_threshold)
anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys = (featureNames(palmieri_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))
anno_grouped <- group_by(anno_palmieri, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)
anno_filtered <- filter(anno_summarized, no_of_matches > 1)

head(anno_filtered)
probe_stats <- anno_filtered 

nrow(probe_stats)
ids_to_exlude <- (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)

table(ids_to_exlude)
palmieri_final <- subset(palmieri_manfiltered, !ids_to_exlude)

validObject(palmieri_final)
fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))
fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID 

validObject(palmieri_final)

individual <- 
  as.character(Biobase::pData(palmieri_final)$Characteristics.individual.)

tissue <- str_replace_all(Biobase::pData(palmieri_final)$Factor.Value.phenotype.,
                          " ", "_")

tissue <- ifelse(tissue == "non-inflamed_colonic_mucosa",
                 "nI", "I")

disease <- 
  str_replace_all(Biobase::pData(palmieri_final)$Factor.Value.disease.,
                  " ", "_")
disease <- 
  ifelse(str_detect(Biobase::pData(palmieri_final)$Factor.Value.disease., 
                    "Crohn"), "CD", "UC")
i_CD

i_CD <- individual[disease == "CD"]
tissue
disease
design_palmieri_CD <- model.matrix(~ 0 + tissue[disease == "CD"] + i_CD)
colnames(design_palmieri_CD)[1:2] <- c("I", "nI")
rownames(design_palmieri_CD) <- i_CD 

i_UC <- individual[disease == "UC"]
design_palmieri_UC <- model.matrix(~ 0 + tissue[disease == "UC"] + i_UC )
colnames(design_palmieri_UC)[1:2] <- c("I", "nI")
rownames(design_palmieri_UC) <- i_UC 

design_palmieri_CD
design_palmieri_UC
