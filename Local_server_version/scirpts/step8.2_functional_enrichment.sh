#!/bin/bash
#./step8_functional_enrichment.sh genome_version input_directory output_directory scripts_path
#
#2022-04-22
#Haojie Chen

genome_version=$1
input_dir=$2
output_dir=$3
scripts_path=$4


Rscript "$scripts_path/Differential_motifscan_enrichment_analysis.R" --input=$input_dir --species=$genome_version --outdir=$output_dir --adjusted_p_value_cutoff=0.1
