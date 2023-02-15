#!/bin/bash


dataset_name=$1
dir="/data/musers/elnesro/MS-data/Datasets/$dataset_name"

cd $dir/FASTP

N=$3

if [ "$2" == "single" ]; then
    for file in $dir/FASTERQ/*.fastq; do
    (
    file="${file%.*}"
    file="${file#*/*/*/*/*/*/*/*/}"
    time ~/fastp -i $dir/FASTERQ/$file".fastq" -o $file"_fp.fq" -h $file"_fp.html" -j $file"_fp.json" -V -x -l 36 -r -w 4
    rm $file"_fp.json"
        
        echo "starting task $file.."
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
   for file in $dir/FASTERQ/*_1.fastq; do
    (
    file="${file%_*}"
    file="${file#*/*/*/*/*/*/*/*/}"
    time ~/fastp -i $dir/FASTERQ/$file"_1.fastq" -I $dir/FASTERQ/$file"_2.fastq" -o $file"_1_fp.fq" -O $file"_2_fp.fq" -h $file"_fp.html" -j $file"_fp.json" -V -x -l 36 -r -w 4
    rm $file"_fp.json"
    
        echo "starting task $file.."
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

