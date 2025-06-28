#!/bin/bash
#./step2.2_reads_mapping_and_removing_duplicates.sh  genome_version name shift_size input_file output_directory
#
#2022-04-22s
#Haojie Chen

genone_version=$1
name=$2
shift_size=$3
drm_bed=$4
output_dir=$5

script_dir="/picb/rsgeno/chenhj/test_docker_images/scripts"
annotation_file="/picb/rsgeno/chenhj/test_docker_images/genomic_annotations/$genone_version/$genone_version.refGene.gtf"
python "$script_dir/bed_coverage_heatmap.py" --input=$drm_bed --ref=$annotation_file --shift_size=$shift_size --outdir=$output_dir --name=$name

