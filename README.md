# multiple-sclerosis-analysis

This is a broad overview of how the computational pipeline works:

# RNA-sequencing:

For each dataset, the SRR numbers for all of the samples in the dataset need to be collected. The range of the SRR numbers is an acceptable input (i.e. SRR100 - SRR200 will download and parse all of the sample records for that range).

./gene expression/RNA-seq/download_srr.sh will use the NCBI SRA Toolkit to download each of the SRR sample data and open it for further processing.

./gene expression/RNA-seq/fastp.sh will take each of these files and perform read trimming and quality control. It also creates interactive plots for examining the quality before and after trimming and adapter removal.

./gene expression/RNA-seq/HISAT2.sh will then align each of these files to the GENCODE hg38 genome.

./gene expression/RNA-seq/htseq-count.sh will take the output .bam files and create a file with the number of each expressed gene for each sample.

./other/gene_counter.py will combine all of the files into one matrix of sample vs. gene, with each cell signifying the amount of expression of that gene (row) in that sample (column).

./gene expression/DEG analysis/DESeq2.R will take all the gene count matrix and perform differential gene expression analysis. As of now, it needs to be modified for each experiment design and other factors, but I am planning on improving this later. This will output several files: a volcano plot for the change in gene expression, the list of significantly up and down regulated genes, the total result of the analysis, and a ranked list of each of the genes with the 'stat' value that DESeq2 outputs.

# Microarray:

Right now, I have only created a script for analyzing microarray data from this type of array: Affymetrix Human Genome U133A. I currently plan on expanding it further. However, the ./gene expression/microarray/MyMicroarray.R file is able to analyze microarray and created a ranked list of differentially expressed genes from the analysis.

# Rank Aggregation:

This is the step where the ranked lists are combined. ./rank aggregation/BiG.R takes the ranked lists output from the differential gene expression analysis and performs rank aggregation on it to produce a single list.

# Pathway Analysis:

./pathway analysis/clusterProfiler.R performs gene set enrichment analysis on the Molecular Signatures Database for each of the input ranked lists.
./pathway analysis/pathway_analysis.R (I'm so creative at coming up with names) performs the same gene set enrichment analysis algorithm on the following databases: Gene Ontology, Reactome, Kyoto Encyclopedia of Genes and Genomes, Kyoto Encylopedia of Genes and Genomes Modules, Disease Ontology, Disease Gene Network, and WikiPathways.

For each analysis a dot plot and file for the enriched / depleted gene sets and pathways are created. 

# Protein Analysis:

Coming soon.
