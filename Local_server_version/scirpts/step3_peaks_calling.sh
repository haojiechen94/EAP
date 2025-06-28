#!/bin/bash
#./step3_peaks_calling.sh ATAC|ChIPPE|ChIPSE genome_version name input_directory output_directory
#
#2022-04-22
#Haojie Chen

sequencing_type=$1
genone_version=$2
name=$3
input_dir=$4
output_dir=$5

if [ "$genone_version" = "hg19" ] || [ "$genone_version" = "hg38" ];
    then
        genome="hs"
elif [ "$genone_version" = "mm9" ] || [ "$genone_version" = "mm10" ];
    then
        genome="mm"
else
    echo "Unknown genome"
    exit 1
fi

echo "Sequencing type: $sequencing_type"
echo "Sample name: $name"
echo "Genome version: $genome"

if [ "$sequencing_type" = "ATAC" ];
    then
        temp_dir="$output_dir/$name";
        mkdir $temp_dir;
        drm_bed="$input_dir/$name/$name""_ss1_drm.bed"
        macs_output="$temp_dir/$name.macs_output"

        cd $temp_dir && macs -t $drm_bed -n $name -f BED -g $genome --nomodel --shiftsize=1 --keep-dup=all > $macs_output 2>&1
elif [ "$sequencing_type" = "ChIPPE" ] || [ "$sequencing_type" = "ChIPSE" ];
    then
        temp_dir="$output_dir/$name";
        mkdir $temp_dir;
        treatment_drm_bed="$input_dir/$name/treatment/$name""_treatment_ss100_drm.bed"
        control_drm_bed="$input_dir/$name/control/$name""_control_ss100_drm.bed"
        macs_output="$temp_dir/$name.macs_output"
        cd $temp_dir && macs -t $treatment_drm_bed -c $control_drm_bed -n $name -f BED -g $genome --nomodel --shiftsize=100 --keep-dup=all > $macs_output 2>&1
else
    echo "Unknown sequencing type"
    exit 1
fi        