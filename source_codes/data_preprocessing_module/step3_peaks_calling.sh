#!/bin/bash
#./step3_peaks_calling.sh metainfo input_directory output_directory ref_index ATAC|ChIPPE|ChIPSE
#
#2022-05-05
#Haojie Chen/Zhijie Guo

meta_info=${1}
input_dir=${2}
output_dir=${3}
ref_idx=${4}
seq_type=${5}

file_info=(`echo $meta_info |cut -d ',' -f 1,2,3 |tr ',' ' '`)
name=${file_info[0]}

if [ "$ref_idx" = "hg19" ] || [ "$ref_idx" = "hg38" ];
    then
        genome="hs"
elif [ "$ref_idx" = "mm9" ] || [ "$ref_idx" = "mm10" ];
    then
        genome="mm"
else
    echo "Unknown genome"
    exit 1
fi

if [ "$seq_type" = "ATAC" ];
    then
        temp_dir="$output_dir/step3_peaks_calling/$name";
        mkdir -m 775 -p $temp_dir;
        drm_bed="$input_dir/step2_mapping/$name/$name""_ss1_drm.bed"
        macs_output="$temp_dir/$name.macs_output"

        cd $temp_dir && macs -t $drm_bed -n $name -f BED -g $genome --nomodel --shiftsize=1 --keep-dup=all > $macs_output 2>&1
elif [ "$seq_type" = "ChIPPE" ] || [ "$seq_type" = "ChIPSE" ];
    then
        temp_dir="$output_dir/step3_peaks_calling/$name";
        mkdir -m 775 -p $temp_dir;
        treatment_drm_bed="$input_dir/step2_mapping/$name/treatment/$name""_treatment_ss100_drm.bed"
        control_drm_bed="$input_dir/step2_mapping/$name/control/$name""_control_ss100_drm.bed"
        macs_output="$temp_dir/$name.macs_output"
        cd $temp_dir && macs -t $treatment_drm_bed -c $control_drm_bed -n $name -f BED -g $genome --nomodel --shiftsize=100 --keep-dup=all > $macs_output 2>&1
else
    echo "Unknown sequencing type"
    exit 1
fi
