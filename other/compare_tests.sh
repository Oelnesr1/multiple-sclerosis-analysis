#!/bin/bash


hg38_DESeq2_up=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation_hg38_gencode.DESeq2_Up.tsv
hg38_DESeq2_down=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation_hg38_gencode.DESeq2_Down.tsv

DESeq2_up=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation.DESeq2_Up.tsv
DESeq2_down=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation.DESeq2_Down.tsv

hg38_edger_up=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation_hg38_gencode.EdgeR_Up.tsv
hg38_edger_down=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation_hg38_gencode.EdgeR_Down.tsv

edger_up=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation.EdgeR_Up.tsv
edger_down=/Users/omarelnesr/MS-Project/Datasets/B_cell_hypomethylation/B_cell_hypomethylation.EdgeR_Down.tsv

their_upregulated=/Users/omarelnesr/Documents/PNAS_B_Cell_Study_Upregulated.tsv
their_downregulated=/Users/omarelnesr/Documents/PNAS_B_Cell_Study_Downregulated.tsv

gene_converter=/Users/omarelnesr/gsea_home/output/c2/B_cell_hypomethylation_c2.GseaPreranked.1670296392015/Symbol_to_probe_set_mapping_details.tsv

study_down=/Users/omarelnesr/Documents/PNAS_B_Cell_Study_Downregulated_Converted.tsv
study_up=/Users/omarelnesr/Documents/PNAS_B_Cell_Study_Upregulated_Converted.tsv

my_converted=/Users/omarelnesr/B_cell_hypomethylation_converted.tsv


awk  'FNR==NR{a[$1]; next} ($1) in a{print $3}' $their_downregulated $gene_converter | sort -u > $study_down
awk  'FNR==NR{a[$1]; next} ($1) in a{print $3}' $their_upregulated $gene_converter | sort -u > $study_up

awk  'FNR==NR{a[$1]; next} ($3) in a{print $1}' $my_upregulated $gene_converter > $my_converted
awk  'FNR==NR{a[$1]; next} ($1) in a{print $3}' $my_downregulated $gene_converter >> $my_converted


num_deseq2_up=( $(wc -l $DESeq2_up) )
num_deseq2_down=( $(wc -l $DESeq2_down) )


num_hg38_deseq2_up=( $(wc -l $hg38_DESeq2_up) )
num_hg38_deseq2_down=( $(wc -l $hg38_DESeq2_down) )


num_hg38_edger_up=( $(wc -l $hg38_edger_up) )
num_hg38_edger_down=( $(wc -l $hg38_edger_down) )


num_edger_up=( $(wc -l $edger_up) )
num_edger_down=( $(wc -l $edger_down) )


num_study_up=( $(wc -l $study_up) )
num_study_down=( $(wc -l $study_down) )

common_hg38_deseq2_deseq2_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_up $DESeq2_up | wc -l) )
common_hg38_deseq2_hg38_edger_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_up $hg38_edger_up | wc -l) )
common_hg38_deseq2_edger_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_up $edger_up | wc -l) )
common_hg38_deseq2_study_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_up $study_up | wc -l) )
common_deseq2_hg38_edger_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $DESeq2_up $hg38_edger_up | wc -l) )
common_deseq2_edger_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $DESeq2_up $edger_up | wc -l) )
common_deseq2_study_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $DESeq2_up $study_up | wc -l) )
common_hg38_edger_edger_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $edger_up $hg38_edger_up | wc -l) )
common_hg38_edger_study_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $study_up $hg38_edger_up | wc -l) )
common_edger_study_up=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $study_up $edger_up | wc -l) )

common_hg38_deseq2_deseq2_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_down $DESeq2_down | wc -l) )
common_hg38_deseq2_hg38_edger_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_down $hg38_edger_down | wc -l) )
common_hg38_deseq2_edger_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_down $edger_down | wc -l) )
common_hg38_deseq2_study_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $hg38_DESeq2_down $study_down | wc -l) )
common_deseq2_hg38_edger_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $DESeq2_down $hg38_edger_down | wc -l) )
common_deseq2_edger_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $DESeq2_down $edger_down | wc -l) )
common_deseq2_study_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $DESeq2_down $study_down | wc -l) )
common_hg38_edger_edger_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $edger_down $hg38_edger_down | wc -l) )
common_hg38_edger_study_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $study_down $hg38_edger_down | wc -l) )
common_edger_study_down=( $(awk 'FNR==NR {x[$1];next} ($1 in x)' $study_down $edger_down | wc -l) )

printf "Analysis\tUp\tDown\n"
printf "hg38 DESeq2\t$num_hg38_deseq2_up\t$num_hg38_deseq2_down\n"
printf "DESeq2\t\t$num_deseq2_up\t$num_deseq2_down\n"
printf "EdgeR\t\t$num_edger_up\t$num_edger_down\n"
printf "hg38 EdgeR\t$num_hg38_edger_up\t$num_hg38_edger_down\n"
printf "Study\t\t$num_study_up\t$num_study_down\n\n"

printf "Common Genes\n\n"

printf "Analysis\t\t\tUp\tDown\n"

printf "hg38 DESeq2 + DESeq2\t\t$common_hg38_deseq2_deseq2_up\t$common_hg38_deseq2_deseq2_down\n"
printf "hg38 DESeq2 + hg38 EdgeR\t$common_hg38_deseq2_hg38_edger_up\t$common_hg38_deseq2_hg38_edger_down\n"
printf "hg38 DESeq2 + EdgeR\t\t$common_hg38_deseq2_edger_up\t$common_hg38_deseq2_edger_down\n"
printf "DESeq2 + hg38 EdgeR\t\t$common_deseq2_hg38_edger_up\t$common_deseq2_hg38_edger_down\n"
printf "DESeq2 + EdgeR\t\t\t$common_deseq2_edger_up\t$common_deseq2_edger_down\n"
printf "hg38 EdgeR + EdgeR\t\t$common_hg38_edger_edger_up\t$common_hg38_edger_edger_down\n"
printf "hg38 DESeq2 + Study\t\t$common_hg38_deseq2_study_up\t$common_hg38_deseq2_study_down\n"
printf "DESeq2 + Study\t\t\t$common_deseq2_study_up\t$common_deseq2_study_down\n"
printf "hg38 EdgeR + Study\t\t$common_hg38_edger_study_up\t$common_hg38_edger_study_down\n"
printf "EdgeR + Study\t\t\t$common_edger_study_up\t$common_edger_study_down\n"
