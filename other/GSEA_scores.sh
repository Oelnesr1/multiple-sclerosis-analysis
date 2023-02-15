#!/bin/bash


BaseDir='/home/elnesro/gsea_home/output'
big_sets=(h c1 c2 c3 c4 c5 c6 c7 c8)
big_sets=(h)
array=()
gsea_report="gsea_report_for_na_"

cd ~
header=('Gene Collections' 'Cell Type' 'Gene Set Name' 'Size' 'Enrichment Score (ES)' 'Normalized Enrichment Score (NES)' 'Nominal p-value' 'FDR q-value' 'FWER p-value')
for ((i=0;i<9;i++));
do
    printf "${header[$i]}\t"
done > ENCODE_GSEA_Rank.txt
cd "$BaseDir"
for big_set in ${big_sets[@]} ;
do
    cd "$BaseDir"/$big_set
    for cell_type in *;
    do
        cd "$BaseDir"/$big_set/$cell_type
        cell_type21=($(echo $cell_type | ggrep -Po '.*?(?=\.)'))
        cell_type2=${cell_type21[0]}
        for file in $gsea_report*;
        do
            html=$(cat "./$file")
            extracted_beginning=($(echo $html | ggrep -Po '(?<=(</a></td><td>)).*?(?=</td><td>t)'))
            extracted=($(echo $html | ggrep -Po '(?<=(<tr>)).*?(?=</td><td>t)'))
            if [[ ! -z "${extracted[@]}" ]]
            then               
                # printf "\n${extracted[13]}" >> /Users/omarelnesr/gsea_output1.txt
                # extracted_new=("${extracted[13]}" "${extracted[18]}" "${extracted[21]}" "${extracted[24]}" "${extracted[27]}")
                for ((i=1;i<200000;i++));
                do
                if [[ ! -z "${extracted[$i]}" ]]
                then
                    # printf "\n${extracted[$i]} $i"
                    if (( $i < 100 ))
                    then
                        if (( $i % 5 == 2 ))
                        then
                            gene_set_name_20=$(echo ${extracted[$i]} | ggrep -Po '(?<=(>)).*?(?=</a>)')
                            printf "\n$big_set\t$cell_type2\t$gene_set_name_20"
                        fi
                        if (( $i % 5 == 4 ))
                        then
                            statistics=($(echo ${extracted[$i]} | ggrep -Po '(?<=(<td>)).*?(?=</td>)'))
                            for stat in ${statistics[@]}; do
                            printf "\t$stat"
                            done
                        fi
                    fi
                    if (( $i > 100 || $i == 100 )) &&  (( $i % 2 == 1 ))
                    then
                        all_values=($(echo ${extracted[$i]} | ggrep -Po '(?<=(<td>)).*?(?=</td>)'))
                        for value in ${all_values[@]}; do
                        if [[ $value == ${all_values[0]} ]]
                        then
                            printf "\n$big_set\t$cell_type2\t$value"
                        else
                            printf "\t$value"
                        fi
                        done
                    fi

                    # if (( $i % 3 == 1 ))
                    # then   
                    #     gene_set_name=$(echo "${extracted[$i]}" | ggrep -Po '(?<=(href=)).*?(?=.html)')
                    #     echo ${gene_set_name:1}
                    # fi

                    # printf "\n${extracted[$i]} $i"
                fi
                done
            fi
        done
    done
done >> /home/elnesro/ENCODE_GSEA_Rank.txt
