# EAP: a cloud-based platform for scalable and comprehensive ChIP/ATAC-seq data analysis
# Url: https://www.biosino.org/epigenetics
## Schematic overview of EAP architecture and its analysis function
![workflow](https://github.com/haojiechen94/EAP/blob/main/images/a.png)
Inputs to EAP include raw sequencing data (in FASTQ format) and metadata (describes study design and sample phenotypes, in CSV format). EAP consists of Data preprocessing module and Downstream analysis module. Data preprocessing module performs quality control, read mapping, peaks calling, creates a summary report for filtering poor quality samples and generates analysis-ready count tables, which are then used as inputs for Downstream analysis module. The Downstream analysis module contains various analytical tools for ChIP/ATAC-seq data analyses, which produce publication-ready results (figures and tables). 


Url: https://www.biosino.org/epigenetics
