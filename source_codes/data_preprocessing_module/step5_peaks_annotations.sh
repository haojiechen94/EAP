#!/bin/bash
#./step5_peaks_annotations.sh input_directory output_directory script_dir reference_directory ref_index
#
#2022-05-05
#Haojie Chen/Zhijie Guo


input_dir=${1}
output_dir=${2}
script_dir=${3}
ref_dir=${4}
ref_idx=${5}

annotation_file="$ref_dir/genomic_annotations/$ref_idx/$ref_idx.refGene.gtf"

python "$script_dir/link_peaks_to_genes.py" --pathname="$input_dir/step3_peaks_calling/*/*_peaks.bed" --ref=$annotation_file --outdir="$output_dir/step5_peaks_annotations/link_peaks_to_genes/"
python "$script_dir/peaks_annotation.py" --pathname="$input_dir/step3_peaks_calling/*/*_peaks.bed" --ref=$annotation_file --outdir="$output_dir/step5_peaks_annotations/peaks_annotation/"
