#!/bin/bash

dir=/data/musers/elnesro/MS-data/Datasets/$1
gtf=/home/elnesro/GENCODE/gencode.v42.annotation.gtf

cd $dir/Counts

N=$3

for file in $dir/HISAT2/*.bam; do
    (
    file="${file%.*}"
    file="${file#*/*/*/*/*/*/*/*/}"

        echo "starting task $file"
        python -m HTSeq.scripts.count -f bam -s $2 -r pos $dir/HISAT2/$file.bam $gtf > $file.counts
        echo "finished task $file!"

        sleep $(( (RANDOM % 3) + 1))
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi

done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

echo "all done"

