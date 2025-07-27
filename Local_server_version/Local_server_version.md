## Local deployment of EAP

We provide the scripts and pre-configured Docker images for the local deployment of EAP, enabling users to run the analysis pipelines without relying on cloud infrastructure.

Step1 Download docker images from Zenodo:.

Step2 Download reference data from: .

Step3 Download scripts from Github[https://github.com/haojiechen94/EAP/tree/main/Local_server_version/scirpts].

Step4 Running analysis on your data using eap_pipeline.sh.

Here is the excution command:
"./eap_pipeline.sh ATAC|ChIPPE|ChIPSE genome_version typical_bin_size metadata.csv input_directory output_diectory variable_of_interest analysis_name scripts_path"

ATAC|ChIPPE|ChIPSE: Choosing a sequencing library type for your data.

genome_version:   Now support human (hg19 and hg38) and mouse (mm9 and mm10).

typical_bin_size: For narrow peaks, suggesting 1000 for ATAC-seq and 2000 for H3K27ac-seq, larger values for broad peaks.

metadata.csv: metadata file descripts the details of study design and sample information.

input_directory/output_diectory: Directory for input data and output data.

variable_of_interest: One of the columns of metadata, performing differential analyssi on this variable.

scripts_path: The directory of scirpts.
