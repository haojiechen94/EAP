#!/bin/bash
#./step2.1_reads_mapping_and_removing_duplicates.sh ATAC genome_version name treatment_R1.fq.gz treatment_R2.fq.gz input_directory output_directory
#
#./step2.1_reads_mapping_and_removing_duplicates.sh ChIPPE genome_version name treatment_R1.fq.gz treatment_R2.fq.gz control_R1.fq.gz control_R2.fq.gz input_directory output_directory
#
#./step2.1_reads_mapping_and_removing_duplicates.sh ChIPSE genome_version name treatment.fq.gz control.fq.gz input_directory output_directory
#
#2022-04-22
#Haojie Chen

sequencing_type=$1
genone_version=$2
name=$3
threads_number=20
script_dir="/picb/rsgeno/chenhj/test_docker_images/scripts"

echo "Sequencing type: $sequencing_type"
echo "Sample name: $name"
echo "Genome version: $genone_version"

bowtie_index="/picb/rsgeno/chenhj/test_docker_images/bowtie_index/$genone_version/$genone_version"
annotation_file="/picb/rsgeno/chenhj/test_docker_images/genomic_annotations/$genone_version/$genone_version.refGene.gtf"

if [ "$sequencing_type" = "ATAC" ];
    then
        input_dir=$6
        arr=(`echo $4| tr '.' ' '`)
        treatment_read1_fq_gz="$input_dir/$name/${arr[0]}_val_1.fq"
        arr=(`echo $5| tr '.' ' '`)
        treatment_read2_fq_gz="$input_dir/$name/${arr[0]}_val_2.fq"

        output_dir=$7
        temp_dir="$output_dir/$name"
        mkdir $temp_dir
        mr="$temp_dir/$name.mr"
        mapstats="$temp_dir/$name.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -X 2000 -p $threads_number --seed 99 --suppress 6,7 $bowtie_index -1 $treatment_read1_fq_gz -2 $treatment_read2_fq_gz > $mr 2> $mapstats

        dupstats="$temp_dir/$name.dupstats"
        drm_bed="$temp_dir/$name""_ss1_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --DNase --paired --final-ss=1 --stats=$dupstats $mr > $drm_bed
elif [ "$sequencing_type" = "ChIPPE" ];
    then
        input_dir=$8
        arr=(`echo $4| tr '.' ' '`)
        treatment_read1_fq_gz="$input_dir/$name/treatment/${arr[0]}_val_1.fq"
        arr=(`echo $5| tr '.' ' '`)
        treatment_read2_fq_gz="$input_dir/$name/treatment/${arr[0]}_val_2.fq"
        arr=(`echo $6| tr '.' ' '`)
        control_read1_fq_gz="$input_dir/$name/control/${arr[0]}_val_1.fq"
        arr=(`echo $7| tr '.' ' '`)
        control_read2_fq_gz="$input_dir/$name/control/${arr[0]}_val_2.fq"

        output_dir=$9
        temp_dir="$output_dir/$name"
        mkdir $temp_dir

        treatment_dir="$temp_dir/treatment"
        mkdir $treatment_dir
        mr="$treatment_dir/$name""_treatment.mr"
        mapstats="$treatment_dir/$name""_treatment.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $threads_number --seed 99 --suppress 6,7 $bowtie_index -1 $treatment_read1_fq_gz -2 $treatment_read2_fq_gz > $mr 2> $mapstats

        dupstats="$treatment_dir/$name""_treatment.dupstats"
        drm_bed="$treatment_dir/$name""_treatment_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --paired --final-ss=100 --stats=$dupstats $mr > $drm_bed

        control_dir="$temp_dir/control"
        mkdir $control_dir
        mr="$control_dir/$name""_control.mr"
        mapstats="$control_dir/$name""_control.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $threads_number --seed 99 --suppress 6,7 $bowtie_index -1 $control_read1_fq_gz -2 $control_read2_fq_gz > $mr 2> $mapstats

        dupstats="$control_dir/$name""_control.dupstats"
        drm_bed="$control_dir/$name""_control_ss100_drm.bed"        
        python "$script_dir/bowtie_mr_check_dremove.py" --paired --final-ss=100 --stats=$dupstats $mr > $drm_bed
elif [ "$sequencing_type" = "ChIPSE" ];
    then
        input_dir=$6
        arr=(`echo $4| tr '.' ' '`)
        treatment_read_fq_gz="$input_dir/$name/treatment/${arr[0]}_trimmed.fq"
        arr=(`echo $5| tr '.' ' '`)
        control_read_fq_gz="$input_dir/$name/control/${arr[0]}_trimmed.fq"

        output_dir=$7
        temp_dir="$output_dir/$name"
        mkdir $temp_dir
        treatment_dir="$temp_dir/treatment"
        mkdir $treatment_dir
        mr="$treatment_dir/$name""_treatment.mr"
        mapstats="$treatment_dir/$name""_treatment.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $threads_number --seed 99 --suppress 6,7 $bowtie_index $treatment_read_fq_gz  > $mr 2> $mapstats

        dupstats="$treatment_dir/$name""_treatment.dupstats"
        drm_bed="$treatment_dir/$name""_treatment_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --single --exp-ss=150 --final-ss=100 --stats=$dupstats $mr > $drm_bed

        control_dir="$temp_dir/control"
        mkdir $control_dir
        mr="$control_dir/$name""_control.mr"
        mapstats="$control_dir/$name""_control.mapstats"
        bowtie -n 2 --nomaqround --maxbts 200 --chunkmbs 256 -k 1 -m 1 -p $threads_number --seed 99 --suppress 6,7 $bowtie_index $control_read_fq_gz > $mr 2> $mapstats

        dupstats="$control_dir/$name""_control.dupstats"
        drm_bed="$control_dir/$name""_control_ss100_drm.bed"
        python "$script_dir/bowtie_mr_check_dremove.py" --single --exp-ss=150 --final-ss=100 --stats=$dupstats $mr > $drm_bed                
else
    echo "Unknown sequencing type"
    exit 1
fi