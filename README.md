# Haradhvala_et_al_2022

This repository contains code used to produce the results in Haradhvala, Leick, Maurer, Gohil et al. Nature Medicine 2022. 

## Data download

Gene expression matrices for this project are available at GEO (accession GSE197268). Raw data is available on dbGaP (accession phs002922). Clinical and other metadata are available in Supplementary tables 1-2 of our paper. Cell-level metadata and clustered anndata objects are available at the link below.

https://drive.google.com/drive/folders/1vw7J8HqUX22ICZmJ0UjAYEBpVjRJ9U9-?usp=sharing

## Analysis

### Pre-processing and clustering.

The preprocessing/ subfolder contains the code run to produce clustered anndata objects. Due to several steps involving randomization the results if re-run in another environment may be very subtly different.

- IP_clustering.py produces an infusion-product-specific clustering of cells
- pbmc_clustering.py clusters PBMC populations from baseline and day 7 post-treatment (both CAR+ and CAR-)
- global_cluster.py clusters all cells in our dataset together
- retreatment_cluster.py clusters just the six samples from a patient treated with two infusions of tisa-cel

### Pseudobulk differential expression

The pseudobulk folder contains code to perform the limma-voom pseudobulk differential expression tests in our paper. 
- prepare_pseudobulk_matrices.py add counts from cells of the sample subtype/patient to create pseudobulk expression matrices
- run_limm_CAR_vs_nonCAR tests differences in CAR-expressing and non-expressing cells in the infusion products and at D7
- run_limma_D7_vs_ip.R test changes in expression occuring between infusion products and day 7
- run_limma_R_vs_NR.R tests for genes differentially expressed by responders and non-responders in our dataset

### T-cell annotation

We used the following notebooks in the Tcell_analyses/ subfolder to annotate various properties of the T-cells in our dataset

- T_cell_cycle.py scores T-cells for cell cycle
- T_CD4_CD8_classification.ipynb categorizes T-cells into CD4 and CD8 subsets
- T_sorting.ipynb performed "pseudo-flow" sorting of our T-cells into classical differentiation subsets

### Figures

Jupyter notebooks figures/Figure[1-6].ipynb can be run to produce figure panels related to scRNA analysis. 

