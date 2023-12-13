#!/bin/bash
#./step4_motif_enrichment.sh metainfo input_directory output_directory reference_directory ref_index ATAC|ChIPPE|ChIPSE
#
#2022-05-05
#Haojie Chen/Zhijie Guo

meta_info=${1}
input_dir=${2}
output_dir=${3}
ref_dir=${4}
ref_idx=${5}
seq_type=${6}

step4_motif_enrichment_dir="$output_dir/step4_motif_enrichment"
homer_reference="$ref_dir/homer/data/genomes/$ref_idx"

if [ "$seq_type" = "ATAC" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,3 |tr ',' ' '`)
        name=${file_info[0]};
        mkdir -m 775 -p "$step4_motif_enrichment_dir/$name"
        peaks_bed="$input_dir/step3_peaks_calling/$name/$name""_peaks.bed"
        findMotifsGenome.pl $peaks_bed $homer_reference "$step4_motif_enrichment_dir/$name"

elif [ "$seq_type" = "ChIPPE" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,3,4,5 |tr ',' ' '`)
        name=${file_info[0]};
        mkdir -m 775 -p "$step4_motif_enrichment_dir/$name"
        peaks_bed="$input_dir/step3_peaks_calling/$name/$name""_peaks.bed"
        findMotifsGenome.pl $peaks_bed $homer_reference "$step4_motif_enrichment_dir/$name"

elif [ "$seq_type" = "ChIPSE" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,4 |tr ',' ' '`)
        name=${file_info[0]};
        mkdir -m 775 -p "$step4_motif_enrichment_dir/$name"
        peaks_bed="$input_dir/step3_peaks_calling/$name/$name""_peaks.bed"
        findMotifsGenome.pl $peaks_bed $homer_reference "$step4_motif_enrichment_dir/$name"

else
    echo "Unknown sequencing type"
    exit 1
fi

