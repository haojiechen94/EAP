#!/bin/bash
#./step7_differential_analysis.sh input_directory output_directory metacsv script_dir reference_directory ref_index variable_of_interest
#
#2022-05-06
#Haojie Chen/Zhijie Guo

input_dir=${1}
output_dir=${2}
metacsv=${3}
script_dir=${4}
ref_dir=${5}
ref_idx=${6}
variable_of_interest=${7}

Rscript "$script_dir/Pairwise_differential_analysis.R" --input="$input_dir/step6_reads_counting/NA_profile_bins.xls" --metadata=$metacsv --interested_variable=$variable_of_interest  --outdir="$output_dir/step7_differential_analysis" --adjusted_p_value_cutoff=0.1
Rscript "$script_dir/Differential_motifscan_enrichment_analysis.R" --input="$input_dir/step7_differential_analysis" --species=$ref_idx --outdir="$output_dir/step8_functional_enrichment" --adjusted_p_value_cutoff=0.1
