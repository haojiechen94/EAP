#!/bin/bash
#./step7_differential_analysis.sh metadata.csv input_directory output_directory variable_of_interest
#
#2022-04-22
#Haojie Chen

metadata=$1
input_dir=$2
output_dir=$3
variable_of_interest=$4

script_dir="/picb/rsgeno/chenhj/test_docker_images/scripts"
Rscript "$script_dir/Pairwise_differential_analysis.R" --input="$input_dir/step6_reads_counting/NA_profile_bins.xls" --metadata=$metadata --interested_variable=$variable_of_interest  --outdir=$output_dir --adjusted_p_value_cutoff=0.1

