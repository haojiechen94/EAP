#!/bin/bash
#./step6_reads_counting.sh ATAC|ChIPPE|ChIPSE metadata.csv genome_version typical_bin_size input_directory output_directory
#
#2022-04-22
#Haojie Chen

sequencing_type=$1
metadata=$2
genome_version=$3
typical_bin_size=$4
input_dir=$5
output_dir=$6

script_dir="/picb/rsgeno/chenhj/test_docker_images/scripts"

if [ "$sequencing_type" = "ATAC" ];
    then
        python "$script_dir/parameters.py" --peaks="$input_dir/step3_peaks_calling/*/*_peaks.bed" --summits="$input_dir/step3_peaks_calling/*/*_summits.bed" --reads="$input_dir/step2_mapping/*/*_drm.bed" --metadata=$metadata --genome=$genome_version --outdir=$output_dir --typical_bin_size=$typical_bin_size --sequencing_type=ATAC
elif [ "$sequencing_type" = "ChIPPE" ] || [ "$sequencing_type" = "ChIPSE" ];
    then
        python "$script_dir/parameters.py" --peaks="$input_dir/step3_peaks_calling/*/*_peaks.bed" --summits="$input_dir/step3_peaks_calling/*/*_summits.bed" --reads="$input_dir/step2_mapping/*/treatment/*_drm.bed" --metadata=$metadata --genome=$genome_version --outdir=$output_dir --typical_bin_size=$typical_bin_size --sequencing_type=ChIP
else
    echo "Unknown sequencing type"
    exit 1
fi


cd $output_dir && profile_bins --parameters="$output_dir/parameters.txt"

python "$script_dir/separate_proximal_and_distal_peak_regions.py" --input="$output_dir/NA_profile_bins.xls" --ref="/picb/rsgeno/chenhj/test_docker_images/genomic_annotations/$genome_version/$genome_version.refGene.gtf" --outdir=$output_dir




