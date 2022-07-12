import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, ttest_ind
import numpy as np
import pandas as pd
sc.set_figure_params(dpi=120)
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

raw_h5ad_file = "../data/CART_may2021_raw.v3.h5ad"
cell_cycle_genes_file = "../data/regev_lab_cell_cycle_genes.txt"
subtypes_file = "../ALLT_CD_classifications.txt"

adata = read_h5ad_gs(raw_h5ad_file)

adata.obs['subtype'] = pd.read_csv(subtypes_file,sep='\t',index_col=0)
idx = adata.obs['subtype'].isin(['CD4 T','CD8 T'])
adata = adata[idx]

# Subset to genes with sufficient counts for scoring
counts = np.array(adata.X.sum(axis=0)).reshape(-1)
gidx = counts>=10000
adata = adata[:,gidx]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)

cell_cycle_genes = pd.read_csv(cell_cycle_genes_file,header=None)[0]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.obs[['S_score','G2M_score','phase']].to_csv('../data/T_cellcycle_annotations.csv')










