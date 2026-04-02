# bioinformatics-projects
A collection of bioinformatics projects focusing on genomic data analysis, including gene expression analysis, data processing, and visualization using Python/R.

## Project 1: Gene Expression Analysis of Breast Tissue

### Overview
This project analyzes gene expression differences between normal breast tissue and breast cancer samples using a public dataset from GEO (GSE20437).

### Objective
Identify genes that are differentially expressed between control and cancer samples.

### Methods
- Loaded and processed gene expression data
- Defined control vs cancer sample groups
- Calculated mean expression for each group
- Computed expression differences
- Identified top upregulated and downregulated genes
- Visualized distribution of expression differences

### Results
- Identified genes with large expression differences between conditions
- Observed overall distribution of gene expression changes

![Distribution Plot](figures/diff_distribution.png)

### Tools Used
- Python
- pandas
- matplotlib

### Dataset
- GEO accession: GSE20437