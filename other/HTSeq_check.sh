#!/bin/bash

dir=/home/elnesro/HTSeq_counts
cd $dir
for file in *.counts; do
data=$(cat $file)
gene_names=($(awk '{print $1}' $file))
counts=($(awk '{print $2}' $file))
num_lines=$(awk '{print $2}' $file | wc -l)
printf "$file $num_lines ${#gene_names[@]} ${#counts[@]} \n"
# printf "${gene_names[@]} ${counts[@]} \n"

for ((line=0;line<$num_lines;i++)); do
if [[ ${counts[$line]} != '0' ]]; then
printf "$file ${gene_names[$lines]} ${counts[$line]} \n"
fi
done
done

