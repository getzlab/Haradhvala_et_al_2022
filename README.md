# Haradhvala_et_al_2022

This repository contains code used to produce the results in Haradhvala, Leick, Maurer, Gohil et al. Nature Medicine 2022. 

## Data download

Gene expression matrices for this project are available at GEO (accession GSE197268). Raw data is available on dbGaP (accession pending). Clinical and other metadata are available in Supplementary tables 1-2 of our paper.

## Analysis

### Pre-processing and clustering.

The preprocessing/ subfolder contains the code run to produce clustered anndata objects. Due to several steps involving randomization the results if re-run in another environment may be very subtly different, the exact annotated versions used in our analyses are available on GEO.

- IP_clustering.py produces an infusion-product-specific clustering of cells
- pbmc_clustering.py clusters PBMC populations from baseline and day 7 post-treatment (both CAR+ and CAR-)
- global_cluster.py clusters all cells in our dataset together
- retreatment_cluster.py clusters just the six samples from a patient treated with two infusions of tisa-cel

### T-cell annotation

We used the following notebooks in the Tcell_analyses/ subfolder to annotate various properties of the T-cells in our dataset

- T_cell_cycle.ipynb scores T-cells for cell cycle
- T_CD_classification.ipynb categorizes T-cells into CD4 and CD8 subsets
- T_sorting.ipynb performed "pseudo-flow" sorting of our T-cells into classical differentiation subsets

### Figures

Jupyter notebooks figures/Figure[1-6].ipynb can be run to produce figure panels related to scRNA analysis. 

