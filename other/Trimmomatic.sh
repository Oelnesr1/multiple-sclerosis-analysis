#!/bin/bash

dir=$1
number=$2

FASTQDir=$1"/FASTQ"
mkdir $1"/Trimmomatic"
workingDir=$1"/Trimmomatic"
cd $FASTQDir
array=()
if [$number == "paired"]; 
then 
    for file in *_1.fastq;
    do
        array+=("${file%_*}")
    done 

    cd $workingDir

    for file in ${array[@]};
    do
    java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 \
     $FASTQDir"/"$file"_1.fastq" $FASTQDir"/"$file"_2.fastq" $file"_1.P.fq" $file"_1.U.fq" $file"_2.P.fq" $file"_2.U.fq" \
     SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
    done

else 
    for file in *;
    do
        array+=("${file%.*}")
    done

    cd $workingDir

    for file in ${array[@]};
    do
     java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 16 \
     "$FASTQDir/$file.fastq" $file".fq" \
     SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
    done



cd $workingDir

for file in ${array[@]};
do
    java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 \
     $FASTQDir"/"$file"_1.fastq" $FASTQDir"/"$file"_2.fastq" $file"_1.P.fq" $file"_1.U.fq" $file"_2.P.fq" $file"_2.U.fq" \
     SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
done
