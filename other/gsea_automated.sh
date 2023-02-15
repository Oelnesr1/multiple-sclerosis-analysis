#!/bin/bash

cell_types=(CD14_positive_monocyte CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell CD4_positive_alpha_beta_memory_T_cellCD8_positive_alpha_beta_memory_T_cell IgD_negative_memory_B_cell immature_natural_killer_cell naive_B_cell naive_thymus_derived_CD4_positive_alpha_beta_T_cell naive_thymus_derived_CD8_positive_alpha_beta_T_cell)

gene_collections=(h c1 c2 c3 c4 c5 c6 c7 c8)
cmd="/Users/omarelnesr/Downloads/GSEA_4.3.2/gsea-cli.sh"
operation="GSEAPreranked"
collapse="-collapse Collapse"
mode="-mode Abs_max_of_probes"
nperm="-nperm 1000"
rnd_seed="-rnd_seed timestamp"
scoring_scheme="-scoring_scheme weighted"
chip="-chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2022.1.Hs.chip"
create_svgs="-create_svgs false"
include_only_symbols="-include_only_symbols true"
make_sets="-make_sets true"
plot_top_x="-plot_top_x 20"
set_max="-set_max 5000"
set_min="-set_min 1"
zip_report="-zip_report false"


for cell_type in ${cell_type[@]}; do
    for gene_collection in ${gene_collections[@]}; do
        gmx="-gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/$gene_collection.all.v2022.1.Hs.symbols.gmt"
        rnk="-rnk /Users/omarelnesr/MS-Project/DESeq2-Results/GSEA-Compatible/Public_Datasets/$cell_type.rnk"
        rpt_label="-rpt_label $cell_type"_"$gene_collection"
        out_folder = "-out /Users/omarelnesr/gsea_home/output/$gene_collection"
        $cmd $operation $gmx $collapse $mode $nperm $rnd_seed $rnk $rpt_label $chip $create_svgs $include_only_symbols $make_sets $plot_top_x $set_max $set_min $zip_report $out_folder
        printf "GSEA analysis for $cell_type in $gene_collection was successful!"
    done
done

