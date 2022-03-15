# Haradhvala_et_al_2022

Scripts to generate figures presented in paper are stored in figures/
Analyses

Data aggregation
Global pre-processing
CD4/CD8 classification

## Pseudo-bulk differential expression
### Data preparation
### Run limma

## Cell cycle analysis
## CD45RA/RO classification
## Cell subset annotation

## CAR-T sub-clustering

## Re-treatment patient clustering

## IP clustering

Figure 1
Figure 2
Figure 3
Figure 4
Figure 5
Figure 6

```mermaid
flowchart LR
    agg[Data Agg] --> global[Global clustering]
    global ---> fig1[Figure 1]

    global ---> classification[T cell classification]
    classification ---> preppseudo[Pseudobulk calculation]
    
    subgraph ide1 [Pseudobulk]
    preppseudo ---> limma[Limma]
    limma ---> fig2[Figure2]
    end

    classification ---> cellcycle[Cell Cycle Analysis]
    cellcycle ---> fig3[Figure 3]

    cellcycle ---> subclust[CAR-T subclustering]
    subclust ---> fig4[Figure 4]
    
    agg ---> ip[Infusion product clustering]
    ip ---> fig5[Figure 5]
    deng[Deng et al. 2020 re-analysis] ---> fig5
    
    retreat[Re-treatment analysis] ---> fig6[Figure 6]

```

