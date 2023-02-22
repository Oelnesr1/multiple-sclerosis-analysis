import sys
import shlex
import subprocess
cmd = '/Users/omarelnesr/Downloads/GSEA_4.3.2/gsea-cli.sh'
operation_name = 'GSEAPreranked'
gene_set = ['h', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8']
# rank_file = ['CD14_positive_monocyte', 'CD4_positive_CD25_positive_alpha_beta_regulatory_T_cell', 
#              'CD4_positive_alpha_beta_memory_T_cell', 'CD8_positive_alpha_beta_memory_T_cell', 
#              'IgD_negative_memory_B_cell', 'immature_natural_killer_cell', 
#              'naive_B_cell', 'naive_thymus_derived_CD4_positive_alpha_beta_T_cell', 'naive_thymus_derived_CD8_positive_alpha_beta_T_cell']

rank_file = ["B_cell_hypomethylation"]

collapse = '-collapse Collapse'
norm = '-norm meandiv'
nperm = '-nperm 1000'
mode = '-mode Abs_max_of_probes'
rnd_seed = '-rnd_seed timestamp'
scoring_scheme = '-scoring_scheme weighted'
create_svgs = '-create_svgs false'
include_only_symbols = '-include_only_symbols true'
make_sets = '-make_sets true'
plot_top_x = '-plot_top_x 20'   
set_max = '-set_max 5000'
set_min = '-set_min 1'
zip_report = '-zip_report false'
ensemble_chip = '-chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2022.1.Hs.chip'
chip = ensemble_chip
end_command_list = []
for temp_gene_set in gene_set:
    gmx = '-gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/'+temp_gene_set+'.all.v2022.1.Hs.symbols.gmt'
    for temp_rank_file in rank_file:
        temp_output_file = temp_rank_file.replace("-", "_")
        rnk = '-rnk /home/elnesro/GSEA_Compatible/'+temp_rank_file+'.rnk'
        rnk = ' -rnk /Users/omarelnesr/MS-Project/DESeq2-Results/GSEA-Compatible/Public_Datasets/' + temp_rank_file + '.rnk'
        rpt_label = '-rpt_label ' + temp_output_file + '_' + temp_gene_set
        out_folder = '-out /Users/omarelnesr/gsea_home/output/' + temp_gene_set
        end_command = cmd + ' ' + operation_name + ' ' + gmx + ' ' + collapse + ' ' + mode + ' ' + nperm + ' ' + rnd_seed + ' ' + rnk + ' ' + scoring_scheme + ' ' + rpt_label + ' ' + chip + ' ' + create_svgs + ' ' + include_only_symbols + ' ' + make_sets + ' ' + plot_top_x + ' ' + set_max + ' ' + set_min + ' ' + zip_report + ' ' + out_folder
        end_command_list.append(end_command)
        print(temp_gene_set + ' ' + temp_rank_file)
print(end_command)
print(end_command_list[0])
for command in end_command_list:
    i = 1
    subprocess.run(shlex.split(command))
    print('command ' + str(i) + ' has run successfully!')
    i += 1
