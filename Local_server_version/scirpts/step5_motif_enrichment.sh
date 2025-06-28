#!/bin/bash
#./step4_peaks_annotations.sh genome_version name input_directory output_directory
#
#2022-04-22
#Haojie Chen

genome_version=$1
name=$2
input_dir=$3
output_dir="$4/$name"

mkdir $output_dir

peaks_bed="$input_dir/$name/$name""_peaks.bed"
homer_reference="/picb/rsgeno/chenhj/homer/data/genomes/$genome_version"
findMotifsGenome.pl $peaks_bed $homer_reference $output_dir
