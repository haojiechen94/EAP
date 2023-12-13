#!/bin/bash
#/step2_reads_mapping_and_removing_duplicates.sh metainfo input_directory output_directory script_dir reference_directory ref_index ATAC|ChIPPE|ChIPSE cpu_number
#
#2022-05-05
#Haojie Chen/Zhijie Guo
meta_info=${1}
input_dir=${2}
output_dir=${3}
script_dir=${4}
ref_dir=${5}
ref_idx=${6}
seq_type=${7}
cpun=${8}

trim_dir="$input_dir/step1_fastqc_and_trim_galore/trim_galore"

bowtie_index="$ref_dir/bowtie_index/$ref_idx/$ref_idx"
annotation_file="$ref_dir/genomic_annotations/$ref_idx/$ref_idx.refGene.gtf"

if [ "$seq_type" = "ATAC" ]
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,3 |tr ',' ' '`)
        name=${file_info[0]}
        treatment_read1_name=(`echo ${file_info[1]}| tr '.' ' '`)
        treatment_read2_name=(`echo ${file_info[2]}| tr '.' ' '`)
        treatment_read1_fq_gz="$trim_dir/$name/${treatment_read1_name[0]}_val_1.fq"
        treatment_read2_fq_gz="$trim_dir/$name/${treatment_read2_name[0]}_val_2.fq"
        
        temp_dir="$output_dir/step2_mapping/$name"
        mkdir -m 775 -p $temp_dir
        mr="$temp_dir/$name.mr"
        mapstats="$temp_dir/$name.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -X 2000 -p $cpun --seed 99 --suppress 6,7 $bowtie_index -1 $treatment_read1_fq_gz -2 $treatment_read2_fq_gz 2> $mapstats|sed "s/\ .*\//\//g" > $mr
        
        dupstats="$temp_dir/$name.dupstats"
        drm_bed="$temp_dir/$name""_ss1_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --DNase --paired --final-ss=1 --stats=$dupstats $mr > $drm_bed
        python "$script_dir/bed_coverage_heatmap.py" --input=$drm_bed --ref=$annotation_file --shift_size=1 --outdir="$temp_dir/" --name=$name

elif [ "$seq_type" = "ChIPPE" ]
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,3,4,5 |tr ',' ' '`)
        name=${file_info[0]}
        treatment_read1_name=(`echo ${file_info[1]}| tr '.' ' '`)
        treatment_read2_name=(`echo ${file_info[2]}| tr '.' ' '`)
        control_read1_name=(`echo ${file_info[3]}| tr '.' ' '`)
        control_read2_name=(`echo ${file_info[4]}| tr '.' ' '`)
        treatment_read1_fq_gz="$trim_dir/$name/treatment/${treatment_read1_name[0]}_val_1.fq"
        treatment_read2_fq_gz="$trim_dir/$name/treatment/${treatment_read2_name[0]}_val_2.fq"
        control_read1_fq_gz="$trim_dir/$name/control/${control_read1_name[0]}_val_1.fq"
        control_read2_fq_gz="$trim_dir/$name/control/${control_read2_name[0]}_val_2.fq"
        
        temp_dir="$output_dir/step2_mapping/$name"
        mkdir -m 775 -p $temp_dir
        
        treatment_dir="$temp_dir/treatment"
        mkdir -m 775 -p $treatment_dir
        mr="$treatment_dir/$name""_treatment.mr"
        mapstats="$treatment_dir/$name""_treatment.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $cpun -X 500 --seed 99 --suppress 6,7 $bowtie_index -1 $treatment_read1_fq_gz -2 $treatment_read2_fq_gz 2> $mapstats|sed "s/\ .*\//\//g" > $mr
        
        dupstats="$treatment_dir/$name""_treatment.dupstats"
        drm_bed="$treatment_dir/$name""_treatment_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --paired --final-ss=100 --stats=$dupstats $mr > $drm_bed
        python "$script_dir/bed_coverage_heatmap.py" --input=$drm_bed --ref=$annotation_file --shift_size=100 --outdir="$treatment_dir/" --name=$name
        
        # treat control
        control_dir="$temp_dir/control"
        mkdir -m 775 -p $control_dir
        mr="$control_dir/$name""_control.mr"
        mapstats="$control_dir/$name""_control.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $cpun -X 500 --seed 99 --suppress 6,7 $bowtie_index -1 $control_read1_fq_gz -2 $control_read2_fq_gz 2> $mapstats|sed "s/\ .*\//\//g" > $mr
        
        dupstats="$control_dir/$name""_control.dupstats"
        drm_bed="$control_dir/$name""_control_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --paired --final-ss=100 --stats=$dupstats $mr > $drm_bed

elif [ "$seq_type" = "ChIPSE" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,4|tr ',' ' '`)
        name=${file_info[0]};
        treatment_name=(`echo ${file_info[1]}| tr '.' ' '`)
        control_name=(`echo ${file_info[2]}| tr '.' ' '`)
        treatment_read_fq_gz="$trim_dir/$name/treatment/${treatment_name[0]}_trimmed.fq"
        control_read_fq_gz="$trim_dir/$name/control/${control_name[0]}_trimmed.fq"
        
        
        temp_dir="$output_dir/step2_mapping/$name"
        mkdir -m 775 -p $temp_dir
        treatment_dir="$temp_dir/treatment"
        mkdir -m 775 -p $treatment_dir
        mr="$treatment_dir/$name""_treatment.mr"
        mapstats="$treatment_dir/$name""_treatment.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $cpun --seed 99 --suppress 6,7 $bowtie_index $treatment_read_fq_gz  2> $mapstats|sed "s/\ .*\//\//g" > $mr

        dupstats="$treatment_dir/$name""_treatment.dupstats"
        drm_bed="$treatment_dir/$name""_treatment_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --single --exp-ss=150 --final-ss=100 --stats=$dupstats $mr > $drm_bed
        python "$script_dir/bed_coverage_heatmap.py" --input=$drm_bed --ref=$annotation_file --shift_size=100 --outdir="$treatment_dir/" --name=$name
        
        # treat control
        control_dir="$temp_dir/control"
        mkdir -m 775 -p $control_dir
        mr="$control_dir/$name""_control.mr"
        mapstats="$control_dir/$name""_control.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $cpun --seed 99 --suppress 6,7 $bowtie_index $control_read_fq_gz > $mr 2> $mapstats

        dupstats="$control_dir/$name""_control.dupstats"
        drm_bed="$control_dir/$name""_control_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --single --exp-ss=150 --final-ss=100 --stats=$dupstats $mr > $drm_bed
else
    echo "Unknown sequencing type"
    exit 1
fi
