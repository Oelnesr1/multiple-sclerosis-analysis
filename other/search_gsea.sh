# !/bin/bash
cd ~
str=(c1 c2 c3 c4 c5 c6 c7 c8)

wc -l ENCODE_GSEA_total_ranks.txt
printf "Number of Unique Sets:\tNumber of Cells:\tNumber of Actual Sets:\tDifference:\tPercent Similarity:\n"
val=$(grep "HALLMARK" ENCODE_GSEA_total_ranks.txt | awk '{print $3}' | sort -u | wc -l)

val1=$(grep "HALLMARK" ENCODE_GSEA_total_ranks.txt | awk '{print $2}' | sort -u | wc -l)

val2=$(grep "HALLMARK" ENCODE_GSEA_total_ranks.txt | sort -u | wc -l)
val3=$(( $val * $val1 ))
val4=$(( $val3 - $val2 ))
val5=$(echo "$val2 * 100 / $val3" | bc -l)
printf "HALLMARK $val * $val1 \t\t\t - $val2 \t\t = $val4 \t --> $val5\n"

for string in ${str[@]}; do
val=$(grep "$string.*$string" ENCODE_GSEA_total_ranks.txt | awk '{print $3}' | sort -u | wc -l)
val1=$(grep "$string.*$string" ENCODE_GSEA_total_ranks.txt | awk '{print $2}' | sort -u | wc -l)
val2=$(grep "$string.*$string" ENCODE_GSEA_total_ranks.txt | sort -u | wc -l)
val3=$(( $val * $val1 ))
val4=$(( $val3 - $val2 ))
val5=$(echo "$val2 * 100 / $val3" | bc -l)
printf "$string \t $val * $val1 \t\t\t - $val2 \t\t = $val4 \t --> $val5\n"
done

for string in ${str[@]}; do
val=$(grep "$string.*$string" gsea_output1.txt | awk '{print $2}' | sort -u | wc -l)
printf "$string $val\n"
done

for string in ${str[@]}; do
val=$(grep "$string.*$string" gsea_output1.txt | sort -u | wc -l)
printf "$string $val\n"
done


array=$(grep HALLMARK ENCODE_GSEA_total_ranks.txt | awk '{print $3}' | sort -u)
for set_name in ${array[@]}; do
val=$(grep $set_name ENCODE_GSEA_total_ranks.txt | sort -u | wc -l)
echo "$set_name $val" 
done
