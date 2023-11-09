# EAP: a cloud-based platform for scalable and comprehensive ChIP/ATAC-seq data analysis
# Url: https://www.biosino.org/epigenetics ; [Run](https://www.biosino.org/epigenetics) analysis on your own data!
## Schematic overview of EAP architecture and its analysis function
![workflow](https://github.com/haojiechen94/EAP/blob/main/images/a.png)

Inputs to EAP include raw sequencing data (in FASTQ format) and metadata (describes study design and sample phenotypes, in CSV format). EAP consists of Data preprocessing module and Downstream analysis module. Data preprocessing module performs quality control, read mapping, peaks calling, creates a summary report for filtering poor quality samples and generates analysis-ready count tables, which are then used as inputs for Downstream analysis module. The Downstream analysis module contains various analytical tools for ChIP/ATAC-seq data analyses, which produce publication-ready results (figures and tables). 

## Data preprocessing module
![Data preprocessing module](https://github.com/haojiechen94/EAP/blob/main/images/b.png)

Workflow of Data preprocessing module, requires the upload of raw sequencing data and metadata using the Cloud Gene-Client tool. Upon completion, analysis-ready count tables and a summary quality control report are generated for downstream analysis.

## Downstream analysis module
![Downstream analysis module](https://github.com/haojiechen94/EAP/blob/main/images/cd.png)

Two common research scenarios are depicted: (c) In scenarios where ChIP/ATAC-seq datasets have clearly defined sample labels, differential analysis can be used to detect differential signals between samples with different labels. This can be followed by differential TF motif enrichment analysis and differential TF activity analysis to explore the TFs associated with differential binding or open chromatin sites. (d) For datasets without pre-defined sample labels or with highly sophisticated sample labels (covering multiple cellular states or disease types), hypervariable analysis can be applied to identify hypervariable ChIP/ATAC-seq signals across the samples. These signals can then be used for clustering analysis to uncover the underlying heterogeneity structure among the samples. Samples can be grouped into different clusters, and signature genes scoring analysis can annotate these clusters based on established gene sets. Supervised analysis tools can also be used to detect binding/open chromatin sites or transcriptional regulators specific to each cluster.

## Data Set Browser
![Data Set Browser](https://github.com/haojiechen94/EAP/blob/main/images/e.png)

EAP has been successfully applied to analyze ChIP/ATAC-seq data from various cancer epigenomic studies. The processed datasets are available on the Data Set Browser in EAP, and the platform offers an interactive interface for easy visualization of TF activity scores in each dataset. This module allows users to investigate the role of interest of transcriptional regulators in oncogenesis by choosing an appropriate data set.

More details see [Document](https://github.com/haojiechen94/EAP/blob/main/doc/Help%20document.pdf)

Url: https://www.biosino.org/epigenetics
