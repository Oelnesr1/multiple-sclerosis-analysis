#!/bin/bash

dir="/data/musers/elnesro/MS-data/Datasets/$1"
genome=~/hg38/genome

cd $dir/HISAT2
export PATH=$PATH:~/samtools-1.16.1:~/hisat2-2.2.1

N=$3

if [ "$2" == "single" ]; then
    for file in $dir/FASTP/*_fp.html; do
    (
    file="${file%_*}"
    file="${file#*/*/*/*/*/*/*/*/}"

    if [ ! -f $file.bam ]; then

    echo "starting task $file.."

    echo "Mapping..."
    time hisat2 -p 8 -x $genome -U $dir/FASTP/$file"_fp.fq" -S $file.sam |
    echo "Coverting from sam to bam..."
    time samtools view -bSu $file.sam -@ 8 -o $file.unsorted.bam
    echo "Sorting..."
    time samtools sort -@ 8 $file.unsorted.bam -o $file.bam
    echo "Removing unnecessary files..."
    time rm $file.sam
    time rm $file.unsorted.bam

    fi

    echo "task $file is done"
        
        sleep $(( (RANDOM % 3) + 1))
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
    done
elif [ "$2" == "double" ]; then
   for file in $dir/FASTP/*_fp.html; do
    (
    file="${file%_*}"
    file="${file#*/*/*/*/*/*/*/*/}"

    if [ ! -f $file.bam ]; then

    echo "starting task $file.."

    echo "Mapping..."
    time hisat2 -p 8 -x $genome -1 $dir/FASTP/$file"_1_fp.fq" -2 $dir/FASTP/$file"_2_fp.fq" -S $file.sam
    echo "Coverting from sam to bam..."
    time samtools view -bSu $file.sam -@ 8 -o $file.unsorted.bam
    echo "Sorting..."
    time samtools sort -@ 8 $file.unsorted.bam -o $file.bam
    echo "Removing unnecessary files..."
    time rm $file.sam
    time rm $file.unsorted.bam

    fi

    echo "task $file is done"
    
        sleep $(( (RANDOM % 3) + 1))
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
    done 
else
  echo "Error: the second argument must be one of 'single' and 'double'" >&2
  exit 1
fi
