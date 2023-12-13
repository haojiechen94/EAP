#!/bin/bash
#./step8_functional_enrichment.sh input_directory output_directory script_dir reference_directory ref_index
#
#2022-05-05
#Haojie Chen/Zhijie Guo


input_dir=${1}
output_dir=${2}
script_dir=${3}
ref_dir=${4}
ref_idx=${5}

annotation_file="$ref_dir/genomic_annotations/$ref_idx/$ref_idx.refGene.gtf"

if [ "$ref_idx" = "hg19" ] || [ "$ref_idx" = "hg38" ];
    then
        gene_sets="$ref_dir/functional_annotation_libraries/KEGG_2021_Human.gmt"
elif [ "$ref_idx" = "mm9" ] || [ "$ref_idx" = "mm10" ];
    then
        gene_sets="$ref_dir/functional_annotation_libraries/KEGG_2019_Mouse.gmt"
else
    echo "Unknown genome"
fi

python "$script_dir/functional_annotation_of_differential_peaks.py" --input_dir="$input_dir/step7_differential_analysis/" --annotation_file=$annotation_file --gene_sets=$gene_sets --outdir="$output_dir/step8_functional_enrichment"
