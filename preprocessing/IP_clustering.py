import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

input_h5ad = "../data/CART_raw_data.h5ad"
cell_cycle_genes = '../data/regev_lab_cell_cycle_genes.txt'

adata = sc.read(input_h5ad)

adata = adata[adata.obs['timepoint']=='Infusion']

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 15, :]

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

## Regress out cell cycle
adata_cc = adata.copy()
sc.pp.scale(adata_cc)
cell_cycle_genes = pd.read_csv(cell_cycle_genes,header=None)[0]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata_cc, s_genes=s_genes, g2m_genes=g2m_genes)

adata.obs=adata_cc.obs

sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=3, min_disp=0.2)

adata.var.loc[['XIST','RPS4Y1','RPS4Y2','Yescarta','Kymriah'],'highly_variable'] = False
adata.var.loc[adata.var.index.str.match('TR.V|IG.V'),'highly_variable'] = False
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt','S_score','G2M_score'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.external.pp.harmony_integrate(adata,'channel')
sc.pp.neighbors(adata,use_rep='X_pca_harmony')
sc.tl.umap(adata)

sc.tl.leiden(adata,resolution=1.0,key_added = 'leiden_1.0')
sc.tl.leiden(adata,resolution=2.0,key_added='leiden_2.0')

adata.obs['cell_type'] = 'T-conv'
adata.obs.loc[adata.obs['leiden_1.0'] == '12','cell_type'] = 'T-reg'
adata.obs.loc[adata.obs['leiden_1.0'] == '14','cell_type'] = 'Myeloid'
adata.obs.loc[adata.obs['leiden_2.0'] == '24','cell_type'] = 'NK'

adata.write("../data/IP_specific_clustering.h5ad")
adata.obs.to_csv("../data/IP_specific_clustering_obs.csv")






