#!/bin/bash
#./step0_mkdir.sh output_directory
#
#2022-05-02
#Haojie Chen/Zhijie Guo


output_dir=$1
mkdir -m 775 -p $output_dir

# Step1 FastQC and cutting adapters
step1_dir="$output_dir/step1_fastqc_and_trim_galore/"
mkdir -m 775 -p $step1_dir
step1_fastqc_dir="$step1_dir/fastqc"
mkdir -m 775 -p $step1_fastqc_dir
step1_trim_galore_dir="$step1_dir/trim_galore"
mkdir -m 775 -p $step1_trim_galore_dir

# Step2 Reads mapping and removing duplicates
step2_mapping_dir="$output_dir/step2_mapping"
mkdir -m 775 -p $step2_mapping_dir

# Step3 Peaks calling
step3_peaks_calling_dir="$output_dir/step3_peaks_calling"
mkdir -m 775 -p $step3_peaks_calling_dir

# Step4 Motif enrichment
step4_motif_enrichment_dir="$output_dir/step4_motif_enrichment"
mkdir -m 775 -p $step4_motif_enrichment_dir

# Step5 Peaks annotation
step5_peaks_annotation_dir="$output_dir/step5_peaks_annotations"
mkdir -m 775 -p $step5_peaks_annotation_dir
step5_link_peaks_to_genes_dir="$step5_peaks_annotation_dir/link_peaks_to_genes"
mkdir -m 775 -p $step5_link_peaks_to_genes_dir
step5_peaks_annotation_dir1="$step5_peaks_annotation_dir/peaks_annotation"
mkdir -m 775 -p $step5_peaks_annotation_dir1

#Step6 Reads counting
step6_reads_counting_dir="$output_dir/step6_reads_counting"
mkdir -m 775 -p $step6_reads_counting_dir

# Step7 Differential analysis
step7_differential_analysis_dir="$output_dir/step7_differential_analysis"
mkdir -m 775 -p $step7_differential_analysis_dir

# Step8 Functional enrichment
step8_functional_enrichment_dir="$output_dir/step8_functional_enrichment"
mkdir -m 775 -p $step8_functional_enrichment_dir

