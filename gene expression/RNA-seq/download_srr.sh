#!/bin/bash

dataset_name=$1
num_ids=$2
first_id_num=$3

cd ~
export PATH=$PATH:$PWD/sratoolkit.3.0.0-ubuntu64/bin
mkdir /data/musers/elnesro/MS-data/Datasets/$dataset_name
mkdir /data/musers/elnesro/MS-data/Datasets/$dataset_name/FASTERQ

N=$4

for ((i=0;i<$num_ids;i++)); do
    (
        current_id=$(($first_id_num + $i))
        current_full_id="SRR${current_id}"
        echo "starting task $i.."
        prefetch "$current_full_id"
        fasterq-dump /data/musers/elnesro/MS-data/SRA_toolkit/sra/${current_full_id}.sra --outdir /data/musers/elnesro/MS-data/Datasets/$dataset_name/FASTERQ -t /tmp --mem 2G --threads 8 -v
        echo "finished task $i"
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
