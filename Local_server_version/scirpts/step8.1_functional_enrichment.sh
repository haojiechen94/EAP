#!/bin/bash
#./step8_functional_enrichment.sh genome_version input_directory output_directory scripts_path
#
#2022-04-22
#Haojie Chen

genome_version=$1
input_dir=$2
output_dir=$3
scripts_path=$4

annotation_file="$scripts_path/genomic_annotations/$genome_version/$genome_version.refGene.gtf"

if [ "$genome_version" = "hg19" ] || [ "$genome_version" = "hg38" ];
    then
        gene_sets="$scripts_path/functional_annotation_libraries/KEGG_2021_Human.gmt"
elif [ "$genome_version" = "mm9" ] || [ "$genome_version" = "mm10" ];
    then
        gene_sets="$scripts_path/functional_annotation_libraries/KEGG_2019_Mouse.gmt"
else
    echo "Unknown genome"
fi


python "$scripts_path/functional_annotation_of_differential_peaks.py" --input_dir=$input_dir --annotation_file=$annotation_file --gene_sets=$gene_sets --outdir=$output_dir

