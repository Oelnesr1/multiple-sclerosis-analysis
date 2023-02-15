#!/bin/bash

num_files=$1
# the file name must have the directory to which files must be outputted, i.e. /data/musers/b_cells --> b_cells is the file name, but the directory is within it
file_name=$2

Rscript ~/MS-Project/R-Scripts/BiG.R $num_files $file_name ${@:3}
Rscript ~/MS-Project/R-Scripts/fgsea.R "rGEO" $file_name
# Rscript ~/MS-Project/R-Scripts/fgsea.R "BiG" $file_name
