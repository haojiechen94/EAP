#!/bin/bash
#./eap_pipeline.sh ATAC|ChIPPE|ChIPSE genome_version typical_bin_size metadata.csv input_directory output_diectory variable_of_interest analysis_name scripts_path
#
# author: Haojie Chen & Zhijie Guo
# date: 2025-06-27
# email: chenhaojie2017@sinh.ac.cn & guozhijie2018@sinh.ac.cn

sequencing_type=$1
genome_version=$2
typical_bin_size=$3
metadata=$4
input_dir=$5
output_dir=$6
variable_of_interest=$7
analysis_name=$8
scripts_path=$9

echo "Step1 FastQC and cutting adapters"
step1_dir="$output_dir/step1_fastqc_and_trim_galore/"
mkdir $step1_dir
step1_fastqc_dir="$step1_dir/fastqc"
mkdir $step1_fastqc_dir
step1_trim_galore_dir="$step1_dir/trim_galore"
mkdir $step1_trim_galore_dir

if [ "$sequencing_type" = "ATAC" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1,2,3`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               treatment_read1_fq_gz=${arr[1]};
               treatment_read2_fq_gz=${arr[2]};
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path -v $scripts_path:$scripts_path  eap_image1 "$scripts_path/step1_fastqc_and_trim_galore.sh" $sequencing_type $name $treatment_read1_fq_gz $treatment_read2_fq_gz $input_dir $step1_dir
        done

elif [ "$sequencing_type" = "ChIPPE" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1,2,3,4,5`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               treatment_read1_fq_gz=${arr[1]}
               treatment_read2_fq_gz=${arr[2]}
               control_read1_fq_gz=${arr[3]}
               control_read2_fq_gz=${arr[4]}
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path  eap_image1 "$scripts_path/step1_fastqc_and_trim_galore.sh" $sequencing_type $name $treatment_read1_fq_gz $treatment_read2_fq_gz $control_read1_fq_gz $control_read2_fq_gz $input_dir $step1_dir
        done

elif [ "$sequencing_type" = "ChIPSE" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1,2,4`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               treatment_fq_gz=${arr[1]}
               control_fq_gz=${arr[2]}
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path  eap_image1 "$scripts_path/step1_fastqc_and_trim_galore.sh" $sequencing_type $name $treatment_fq_gz $control_fq_gz $input_dir $step1_dir
        done
else
    echo "Unknown sequencing type"
    exit 1
fi

echo "Step2 Reads mapping and removing duplicates"
step2_mapping_dir="$output_dir/step2_mapping"
mkdir $step2_mapping_dir

if [ "$sequencing_type" = "ATAC" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1,2,3`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               treatment_read1_fq_gz=${arr[1]};
               treatment_read2_fq_gz=${arr[2]};
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image1 "$scripts_path/step2.1_reads_mapping_and_removing_duplicates.sh" $sequencing_type $genome_version $name $treatment_read1_fq_gz $treatment_read2_fq_gz "$output_dir/step1_fastqc_and_trim_galore/trim_galore" $step2_mapping_dir
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image2 "$scripts_path/step2.2_reads_mapping_and_removing_duplicates.sh" $genome_version $name 1 "$step2_mapping_dir/$name/$name""_ss1_drm.bed" "$step2_mapping_dir/$name/"
        done

elif [ "$sequencing_type" = "ChIPPE" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1,2,3,4,5`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               treatment_read1_fq_gz=${arr[1]}
               treatment_read2_fq_gz=${arr[2]}
               control_read1_fq_gz=${arr[3]}
               control_read2_fq_gz=${arr[4]}
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image1 "$scripts_path/step2.1_reads_mapping_and_removing_duplicates.sh" $sequencing_type $genome_version $name $treatment_read1_fq_gz $treatment_read2_fq_gz $control_read1_fq_gz $control_read2_fq_gz "$output_dir/step1_fastqc_and_trim_galore/trim_galore" $step2_mapping_dir
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image2 "$scripts_path/step2.2_reads_mapping_and_removing_duplicates.sh" $genome_version $name 100 "$step2_mapping_dir/$name/treatment/$name""_treatment_ss100_drm.bed" "$step2_mapping_dir/$name/treatment/" 
        done

elif [ "$sequencing_type" = "ChIPSE" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1,2,4`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               treatment_fq_gz=${arr[1]}
               control_fq_gz=${arr[2]}
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image1 "$scripts_path/step2.1_reads_mapping_and_removing_duplicates.sh" $sequencing_type $genome_version $name $treatment_fq_gz $control_fq_gz "$output_dir/step1_fastqc_and_trim_galore/trim_galore" $step2_mapping_dir
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image2 "$scripts_path/step2.2_reads_mapping_and_removing_duplicates.sh" $genome_version $name 100 "$step2_mapping_dir/$name/treatment/$name""_treatment_ss100_drm.bed" "$step2_mapping_dir/$name/treatment/"                
        done
else
    echo "Unknown sequencing type"
    exit 1
fi

echo "Step3 Peaks calling"
step3_peaks_calling_dir="$output_dir/step3_peaks_calling"
mkdir $step3_peaks_calling_dir

if [ "$sequencing_type" = "ATAC" ] || [ "$sequencing_type" = "ChIPPE" ] || [ "$sequencing_type" = "ChIPSE" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path   eap_image1 "$scripts_path/step3_peaks_calling.sh" $sequencing_type $genome_version $name "$output_dir/step2_mapping/" $step3_peaks_calling_dir
        done
else
    echo "Unknown sequencing type"
    exit 1
fi


echo "Step4 Peaks annotation"

step4_peaks_annotation_dir="$output_dir/step4_peaks_annotations"
mkdir $step4_peaks_annotation_dir
step4_link_peaks_to_genes_dir="$step4_peaks_annotation_dir/link_peaks_to_genes"
mkdir $step4_link_peaks_to_genes_dir
step4_peaks_annotation_dir1="$step4_peaks_annotation_dir/peaks_annotation"
mkdir $step4_peaks_annotation_dir1

docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap_image2 "$scripts_path/step4_peaks_annotations.sh" $genome_version "$output_dir/step3_peaks_calling/" $step4_peaks_annotation_dir

echo "Step5 Motif enrichment"

step5_motif_enrichment_dir="$output_dir/step5_motif_enrichment"
mkdir $step5_motif_enrichment_dir

if [ "$sequencing_type" = "ATAC" ] || [ "$sequencing_type" = "ChIPPE" ] || [ "$sequencing_type" = "ChIPSE" ];
    then
        for i in `tail -n +2 $metadata | cut -d ',' -f 1`;
            do arr=(`echo $i| tr ',' ' '`);
               name=${arr[0]};
               docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap_image2 "$scripts_path/step5_motif_enrichment.sh" $genome_version $name "$output_dir/step3_peaks_calling/" $step5_motif_enrichment_dir
        done
else
    echo "Unknown sequencing type"
    exit 1
fi

echo "Step6 Reads counting"

step6_reads_counting_dir="$output_dir/step6_reads_counting"
mkdir $step6_reads_counting_dir

docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap_image1 "$scripts_path/step6_reads_counting.sh" $sequencing_type $metadata $genome_version $typical_bin_size "$output_dir" $step6_reads_counting_dir


echo "Step7 Differential analysis"
step7_differential_analysis_dir="$output_dir/step7_differential_analysis"
mkdir $step7_differential_analysis_dir

docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap "$scripts_path/step7_differential_analysis.sh" $metadata $output_dir $step7_differential_analysis_dir $variable_of_interest

echo "Step8 Functional enrichment"

step8_functional_enrichment_dir="$output_dir/step8_functional_enrichment"
mkdir $step8_functional_enrichment_dir
docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap_image2 "$scripts_path/step8.1_functional_enrichment.sh" $genome_version "$output_dir/step7_differential_analysis" $step8_functional_enrichment_dir $scripts_path

docker run --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap "$scripts_path/step8.2_functional_enrichment.sh" $genome_version "$output_dir/step7_differential_analysis" $step8_functional_enrichment_dir $scripts_path


docker run  --rm -v $input_dir:$input_dir -v $output_dir:$output_dir -v $scripts_path:$scripts_path eap_image3 python "$scripts_path/create_report.py" --input_dir=$output_dir --name=$analysis_name --sequencing_type=$sequencing_type --outdir=$output_dir --cover_image="$scripts_path/logo_and_cover_image/cover.png" --logo_image="$scripts_path/logo_and_cover_image/logo.png"



