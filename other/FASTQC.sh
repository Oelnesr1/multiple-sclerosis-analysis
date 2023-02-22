#!/bin/bash

dir=$1
cd /data/musers/elnesro/MS-data/Datasets/$dir

mkdir /data/musers/elnesro/MS-data/Datasets/$dir/FASTQC

cd FASTERQ

~/FastQC/fastqc ${@:2} --outdir /data/musers/elnesro/MS-data/Datasets/$dir/FASTQC/ --threads 16
