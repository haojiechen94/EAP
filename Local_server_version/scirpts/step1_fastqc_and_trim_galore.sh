#!/bin/bash
#./step1_fastqc_and_trim_galore.sh ATAC name treatment_R1.fq.gz treatment_R2.fq.gz input_directory output_directory
#
#./step1_fastqc_and_trim_galore.sh ChIPPE name treatment_R1.fq.gz treatment_R2.fq.gz control_R1.fq.gz control_R2.fq.gz input_directory output_directory
#
#./step1_fastqc_and_trim_galore.sh ChIPSE name treatment.fq.gz control.fq.gz input_directory output_directory
#
#2022-04-22
#Haojie Chen

sequencing_type=$1
name=$2
echo "Sequencing type: $sequencing_type"
echo "Sample name: $name"

if [ "$sequencing_type" = "ATAC" ];
    then
        input_dir=$5
        treatment_read1_fq_gz="$input_dir/$3"
        treatment_read2_fq_gz="$input_dir/$4"

        output_dir=$6
        temp_dir="$output_dir/fastqc/$name/"
        mkdir $temp_dir
        treatment_read1_output_dir="$temp_dir/R1"
        mkdir $treatment_read1_output_dir
        treatment_read2_output_dir="$temp_dir/R2"
        mkdir $treatment_read2_output_dir

        fastqc $treatment_read1_fq_gz --outdir $treatment_read1_output_dir;
        fastqc $treatment_read2_fq_gz --outdir $treatment_read2_output_dir;

        temp_dir="$output_dir/trim_galore/$name/"
        mkdir $temp_dir
        trim_galore --dont_gzip --paired --trim1 -o $temp_dir --fastqc $treatment_read1_fq_gz $treatment_read2_fq_gz;

elif [ "$sequencing_type" = "ChIPPE" ];
    then
        input_dir=$7
        treatment_read1_fq_gz="$input_dir/$3"
        treatment_read2_fq_gz="$input_dir/$4"
        control_read1_fq_gz="$input_dir/$5"
        control_read2_fq_gz="$input_dir/$6"        

        output_dir=$8
        temp_dir="$output_dir/fastqc/$name/"
        mkdir $temp_dir
        temp_dir="$output_dir/fastqc/$name/treatment/"
        mkdir $temp_dir
        treatment_read1_output_dir="$temp_dir/R1"
        mkdir $treatment_read1_output_dir
        treatment_read2_output_dir="$temp_dir/R2"
        mkdir $treatment_read2_output_dir        
        temp_dir="$output_dir/fastqc/$name/control/"
        mkdir $temp_dir                
        control_read1_output_dir="$temp_dir/R1"
        mkdir $control_read1_output_dir
        control_read2_output_dir="$temp_dir/R2"
        mkdir $control_read2_output_dir  
        
        fastqc $treatment_read1_fq_gz --outdir $treatment_read1_output_dir;
        fastqc $treatment_read2_fq_gz --outdir $treatment_read2_output_dir;
        fastqc $control_read1_fq_gz --outdir $control_read1_output_dir;
        fastqc $control_read2_fq_gz --outdir $control_read2_output_dir;        

        temp_dir="$output_dir/trim_galore/$name/"
        mkdir $temp_dir
        temp_dir="$output_dir/trim_galore/$name/treatment/"
        mkdir $temp_dir            
        trim_galore --dont_gzip --paired --trim1 -o $temp_dir --fastqc $treatment_read1_fq_gz $treatment_read2_fq_gz;
        temp_dir="$output_dir/trim_galore/$name/control/"
        mkdir $temp_dir            
        trim_galore --dont_gzip --paired --trim1 -o $temp_dir --fastqc $control_read1_fq_gz $control_read2_fq_gz;        
elif [ "$sequencing_type" = "ChIPSE" ];
    then
        input_dir=$5
        treatment_fq_gz="$input_dir/$3"
        control_fq_gz="$input_dir/$4"

        output_dir=$6
        temp_dir="$output_dir/fastqc/$name/"
        mkdir $temp_dir
        treatment_output_dir="$output_dir/fastqc/$name/treatment/"
        mkdir $treatment_output_dir
        control_output_dir="$output_dir/fastqc/$name/control/"
        mkdir $control_output_dir        

        fastqc $treatment_fq_gz --outdir $treatment_output_dir;
        fastqc $control_fq_gz --outdir $control_output_dir;

        temp_dir="$output_dir/trim_galore/$name/"
        mkdir $temp_dir
        temp_dir="$output_dir/trim_galore/$name/treatment/"
        mkdir $temp_dir            
        trim_galore --dont_gzip -o $temp_dir --fastqc $treatment_fq_gz;
        temp_dir="$output_dir/trim_galore/$name/control/"
        mkdir $temp_dir            
        trim_galore --dont_gzip -o $temp_dir --fastqc $control_fq_gz;
else
    echo "Unknown sequencing type"
    exit 1
fi
