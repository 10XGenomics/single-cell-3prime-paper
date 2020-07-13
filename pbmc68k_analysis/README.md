# Single Cell RNA-seq Secondary Analysis of 68k PBMCs

The R scripts here were used to generate the analysis of 68k PBMC single cell data, which was described in the manuscript, Zheng et al., "Massively parallel digital transcriptional profiling of single cells". 

For any questions about the data and the analysis, please contact support@10xgenomics.com.

## Data

- scRNA-seq data of 68k PBMCs from Donor A
- scRNA-seq data of bead-enriched sub-populations of PBMCs from Donor A 
 
All of the data were processed by the Cell Ranger 1.1 pipeline. The processed data used for the single cell RNA-seq secondary analysis are stored in three R data files:

* [pbmc68k_data.rds](https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/pbmc68k_data.rds) (77MB): consists of the gene expression profiles of ~68k PBMCs
* [all_pure_pbmc_data.rds](https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/all_pure_pbmc_data.rds) (1.5GB): consists of the gene expression profiles of 10 bead-enriched PBMC samples
* [all_pure_select_11types.rds](https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/all_pure_select_11types.rds) (687KB): consists of the gene expression and meta-information of the 11 sub-populations of PBMC identified from the 10 samples

## Scripts

The following R scripts are used to cluster and classify the 68k PBMC scRNA-seq data with that of purified PBMCs. The analysis involves dimension reduction via PCA, visualization via t-SNE, and cluster analysis via k-means clustering. 

* __main_process_pure_pbmc.R__: main script to analyze the 11 bead-enriched sub-populations of PBMC samples and generate reference expression profiles
* __main_process_68k_pbmc.R__:  main script to cluster and classify the 68k PBMC scRNA-seq data with that of purified PBMCs

Functions for basic analysis and figure generation in the scripts above are also provided in:

* __util.R__: general utility functions used by the main scripts for the gene expression analysis 
* __select_pure_pbmc.R__: specific functions used to identify pure cell types in the scRNA-seq data of bead-enriched sub-populations of PBMCs

## Reference

Massively parallel digital transcriptional profiling of single cells, in review, 2016.
