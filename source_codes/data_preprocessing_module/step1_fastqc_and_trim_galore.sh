#!/bin/bash
#./step1_fastqc_and_trim_galore.sh metainfo input_directory output_directory ATAC|ChIPPE|ChIPSE cpu_number
#
#2022-05-02
#Haojie Chen/Zhijie Guo
meta_info=${1}
input_dir=${2}
output_dir=${3}
seq_type=${4}
cpun=${5}

if [ "$seq_type" = "ATAC" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,3 |tr ',' ' '`)
        name=${file_info[0]};
        treatment_read1_fq_gz="$input_dir/${file_info[1]}";
        treatment_read2_fq_gz="$input_dir/${file_info[2]}";
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/";
        mkdir -m 775 -p $temp_dir
        treatment_read1_output_dir="$temp_dir/R1"
        mkdir -m 775 -p $treatment_read1_output_dir
        treatment_read2_output_dir="$temp_dir/R2"
        mkdir -m 775 -p $treatment_read2_output_dir
        
        
        fastqc $treatment_read1_fq_gz --outdir $treatment_read1_output_dir -t $cpun;
        fastqc $treatment_read2_fq_gz --outdir $treatment_read2_output_dir -t $cpun;

        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/"
        mkdir -m 775 -p $temp_dir
        trim_galore --dont_gzip --paired --trim1 -o $temp_dir --fastqc $treatment_read1_fq_gz $treatment_read2_fq_gz;
        

elif [ "$seq_type" = "ChIPPE" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,3,4,5 |tr ',' ' '`)
        name=${file_info[0]};
        treatment_read1_fq_gz="$input_dir/${file_info[1]}"
        treatment_read2_fq_gz="$input_dir/${file_info[2]}"
        control_read1_fq_gz="$input_dir/${file_info[3]}"
        control_read2_fq_gz="$input_dir/${file_info[4]}"
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/"
        mkdir -m 775 -p $temp_dir
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/treatment/"
        mkdir -m 775 -p $temp_dir
        treatment_read1_output_dir="$temp_dir/R1"
        mkdir -m 775 -p $treatment_read1_output_dir
        treatment_read2_output_dir="$temp_dir/R2"
        mkdir -m 775 -p $treatment_read2_output_dir        
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/control/"
        mkdir -m 775 -p $temp_dir                
        control_read1_output_dir="$temp_dir/R1"
        mkdir -m 775 -p $control_read1_output_dir
        control_read2_output_dir="$temp_dir/R2"
        mkdir -m 775 -p $control_read2_output_dir  
        
        fastqc $treatment_read1_fq_gz --outdir $treatment_read1_output_dir -t $cpun;
        fastqc $treatment_read2_fq_gz --outdir $treatment_read2_output_dir -t $cpun;
        fastqc $control_read1_fq_gz --outdir $control_read1_output_dir -t $cpun;
        fastqc $control_read2_fq_gz --outdir $control_read2_output_dir -t $cpun;

        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/"
        mkdir -m 775 -p $temp_dir
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/treatment/"
        mkdir -m 775 -p $temp_dir
        trim_galore --dont_gzip --paired --trim1 -o $temp_dir --fastqc $treatment_read1_fq_gz $treatment_read2_fq_gz;
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/control/"
        mkdir -m 775 -p $temp_dir
        trim_galore --dont_gzip --paired --trim1 -o $temp_dir --fastqc $control_read1_fq_gz $control_read2_fq_gz;
elif [ "$seq_type" = "ChIPSE" ];
    then
        file_info=(`echo $meta_info |cut -d ',' -f 1,2,4 |tr ',' ' '`)
        name=${file_info[0]};
        treatment_fq_gz="$input_dir/${file_info[1]}"
        control_fq_gz="$input_dir/${file_info[2]}"
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/"
        mkdir -m 775 -p $temp_dir
        treatment_output_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/treatment/"
        mkdir -m 775 -p $treatment_output_dir
        control_output_dir="$output_dir/step1_fastqc_and_trim_galore/fastqc/$name/control/"
        mkdir -m 775 -p $control_output_dir        

        fastqc $treatment_fq_gz --outdir $treatment_output_dir -t $cpun;
        fastqc $control_fq_gz --outdir $control_output_dir -t $cpun;

        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/"
        mkdir -m 775 -p $temp_dir
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/treatment/"
        mkdir -m 775 -p $temp_dir            
        trim_galore --dont_gzip -o $temp_dir --fastqc $treatment_fq_gz;
        temp_dir="$output_dir/step1_fastqc_and_trim_galore/trim_galore/$name/control/"
        mkdir -m 775 -p $temp_dir            
        trim_galore --dont_gzip -o $temp_dir --fastqc $control_fq_gz;
else
    echo "Unknown sequencing type"
    exit 1
fi
