#!/bin/bash
#./step4_peaks_annotations.sh genome_version input_directory output_directory
#
#2022-04-22
#Haojie Chen

genome_version=$1
input_dir=$2
output_dir=$3

echo "Genome version: $genome_version"

annotation_file="/picb/rsgeno/chenhj/test_docker_images/genomic_annotations/$genome_version/$genome_version.refGene.gtf"
script_dir="/picb/rsgeno/chenhj/test_docker_images/scripts"

python "$script_dir/link_peaks_to_genes.py" --pathname="$input_dir/*/*_peaks.bed" --ref=$annotation_file --outdir="$output_dir/link_peaks_to_genes/"
python "$script_dir/peaks_annotation.py" --pathname="$input_dir/*/*_peaks.bed" --ref=$annotation_file --outdir="$output_dir/peaks_annotation/"
