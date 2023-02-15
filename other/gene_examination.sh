#!/bin/bash

gsea_dir=/Users/omarelnesr/MS-Project/DESeq2-Results/GSEA-Compatible/GSEA_Dropbox
gsea_dir2=/Users/omarelnesr/MS-Project/DESeq2-Results/Full-Results
ensemble_genes=(ENSG00000117020 ENSG00000006071 ENSG00000100604 ENSG00000077279 ENSG00000197635 ENSG00000109911 ENSG00000125798 ENSG00000150907 ENSG00000152254 ENSG00000115263 ENSG00000106633 ENSG00000135100 ENSG00000121351 ENSG00000254647 ENSG00000173404 ENSG00000016082 ENSG00000135363 ENSG00000204103 ENSG00000162992 ENSG00000122859 ENSG00000125820 ENSG00000163623 ENSG00000077264 ENSG00000106331 ENSG00000007372 ENSG00000175426 ENSG00000125851 ENSG00000139515 ENSG00000143627 ENSG00000079689 ENSG00000140612 ENSG00000163581 ENSG00000114902 ENSG00000140319 ENSG00000143742 ENSG00000144867 ENSG00000157005 ENSG00000136854 ENSG00000019505 ENSG00000111424)
cd $gsea_dir

for file in *.rnk; do
i=0
for gene_name in ${ensemble_genes[@]}; do
val=$(grep $gene_name $file | sort -u | wc -l) 
if (( $val == 1 )); then
i=$(( $i + 1 ))
fi
done
echo $file $i
done

for gene_name in ${ensemble_genes[@]}; do
i=0
for file in *.rnk; do
val=$(grep $gene_name $file | sort -u | wc -l) 
if (( $val > 0 )); then
i=$(( $i + 1 ))
fi
done
printf "$gene_name $i\t"
done