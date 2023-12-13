#!/bin/bash
#./step6_reads_counting.sh input_directory output_directory metacsv script_dir reference_directory ref_index ATAC|ChIPPE|ChIPSE typical_bin_size
#
#2022-05-05
#Haojie Chen/Zhijie Guo

input_dir=${1}
output_dir=${2}
metacsv=${3}
script_dir=${4}
ref_dir=${5}
ref_idx=${6}
seq_type=${7}
typical_bin_size=${8}

step6_dir="$output_dir/step6_reads_counting"
black_list="$ref_dir/black_lists/$ref_idx""_blacklist.bed"


if [ "$seq_type" = "ATAC" ];
    then
        python "$script_dir/parameters.py" --peaks="$input_dir/step3_peaks_calling/*/*_peaks.bed" --summits="$input_dir/step3_peaks_calling/*/*_summits.bed" --reads="$input_dir/step2_mapping/*/*_drm.bed" --metadata=$metacsv --black_list=$black_list --outdir=$step6_dir --typical_bin_size=$typical_bin_size --sequencing_type=ATAC
elif [ "$seq_type" = "ChIPPE" ] || [ "$seq_type" = "ChIPSE" ];
    then
        python "$script_dir/parameters.py" --peaks="$input_dir/step3_peaks_calling/*/*_peaks.bed" --summits="$input_dir/step3_peaks_calling/*/*_summits.bed" --reads="$input_dir/step2_mapping/*/treatment/*_drm.bed" --metadata=$metacsv --black_list=$black_list --outdir=$step6_dir --typical_bin_size=$typical_bin_size --sequencing_type=ChIP
else
    echo "Unknown sequencing type"
    exit 1
fi

cd $step6_dir
profile_bins --parameters="$step6_dir/parameters.txt"

python "$script_dir/separate_proximal_and_distal_peak_regions.py" --input="$step6_dir/NA_profile_bins.xls" --ref="$ref_dir/genomic_annotations/$ref_idx/$ref_idx.refGene.gtf" --outdir=$step6_dir --bed

python "$script_dir/link_peaks_to_genes.py" --pathname="$step6_dir/proximal_peak_regions*.bed" --ref="$ref_dir/genomic_annotations/$ref_idx/$ref_idx.refGene.gtf" --outdir="$step6_dir/"


