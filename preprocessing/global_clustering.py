import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scrublet as scr
sc.set_figure_params(dpi=120)
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

input_h5ad = "gs://ibm-cart-0/analysis/may2021/CART_may2021_raw.v3.h5ad"

# Load data
adata = read_h5ad_gs(input_h5ad)

## Doublet filtering
# Scrublet not reliable for IP since they're almost all T-cells, so we'll just use the results on PBMC channels
obs_pbmc = pd.read_csv('PBMC_specific_clustering_obs.csv',index_col=0)
idx = (adata.obs['timepoint']=='Infusion')|(adata.obs.index.isin(obs_pbmc.index))
adata = adata[idx]

## Other filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 15, :]

## Normalization / HVG
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=3, min_disp=0.2)
adata.var.loc[['XIST','RPS4Y1','RPS4Y2','Yescarta','Kymriah'],'highly_variable'] = False
adata.var.loc[adata.var.index.str.match('TR.V|IG.V'),'highly_variable'] = False

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

# PCA / clustering
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

sc.external.pp.harmony_integrate(adata,'batch')
adata.obsm['X_pca_harmony_sub']=adata.obsm['X_pca_harmony'][:,0:18]
sc.pp.neighbors(adata,use_rep='X_pca_harmony_sub')
sc.tl.umap(adata)

sc.tl.leiden(adata,resolution=1,key_added='leiden_1.0')


cell_types = ['CD4 T','Infusion T','Infusion T','Monocyte','CD4 T','CD8 T','CD8 T',
            'CD8 T','Monocyte','NK','Infusion T','Monocyte','Monocyte','Monocyte','CD8 T','mDC','Infusion T',
            'Platelet','B','Monocyte']

adata.obs['cell_type'] = [cell_types[int(i)] for i in adata.obs['leiden_1.0']]

# Rescue pDCs not isolated in 1.0 clustering
adata.obs.loc[adata.obs['leiden_2.0']=='37','cell_type'] = 'pDC'

adata.write("global_clustering.h5ad")
adata.obs.to_csv("global_clustering_obs.csv")








