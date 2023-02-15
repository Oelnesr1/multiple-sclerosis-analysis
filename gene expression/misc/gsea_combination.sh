#!/bin/bash

dir=$1
file_name=$2
cd $dir


columns_less_than_20=("0" "3" "4" "5" "6" "7" "8" "9")
columns_greater_than_20=("0" "2" "3" "4" "5" "6" "7" "8")

j=0
for file in gsea_report_for_na*.tsv; do
i=1
while IFS=$'\t' read -r -a my_array; do 

    if (($i < 22)); then
        if (($i == 1)); then
            if (($j == 0)); then
                for column in ${columns_less_than_20[@]}; do
                printf -- "${my_array[$column]} \t"
                done
            fi
            else
                for column in ${columns_less_than_20[@]}; do
                printf -- "${my_array[$column]} \t"
                done
            fi
        else
            for column in ${columns_greater_than_20[@]}; do
            printf -- "${my_array[$column]} \t"
            done
    fi
    if (( !($i == 1 & $j != 0) )); then
    printf "\n"
    fi
    i=$(($i + 1))
    
done < "$file" >> $file_name
j=$(($j + 1))
done 
