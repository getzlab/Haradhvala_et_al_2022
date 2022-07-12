import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scrublet as scr
import pickle
import os

# Location of retreatment samples
data_dir = "../data/GSE197268/"
dirs = os.listdir(data_dir)
#retreat_dirs = [x for x in dirs if x.endswith("retreatment")]
retreat_dirs = [x for x in dirs if x.startswith("Patient29")]

# Other samples
full_h5ad = "../data/CART_raw_data.h5ad"


adatas=list()
scrubs = dict()

for sample in retreat_dirs:
    print(sample)
    
    adata = sc.read_10x_mtx(data_dir + sample)
    adata.var_names_make_unique()
    
    adata.obs['channel']=sample
    adata.obs['sample']=sample
    adata.obs['patient']='Patient29'
    
    # doublet annotation
    scrub = scr.Scrublet(adata.X)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
    
    scrubs[sample] = scrub
    
    adata.obs['treatment'] = 'second' if sample.endswith("retreatment") else 'first'

    adatas.append(adata)

adata = adatas[0].concatenate(adatas[1:])

# Filtering
adata = adata[adata.obs['predicted_doublets']==False,:]

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 15, :]

adata.write("../data/retreatment_patient_raw.h5ad")
 
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata_cc = adata.copy()
sc.pp.scale(adata_cc)
cell_cycle_genes = pd.read_csv('../data/regev_lab_cell_cycle_genes.txt',header=None)[0]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata_cc, s_genes=s_genes, g2m_genes=g2m_genes)

adata.obs=adata_cc.obs


# Variable gene selection
sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=3, min_disp=0.2)
adata.raw = adata
adata.var.loc[['Kymriah'],'highly_variable'] = False
adata.var.loc[adata.var.index.str.match('TR.V|IG.V'),'highly_variable'] = False
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt','S_score','G2M_score'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata,n_pcs=17)
sc.tl.umap(adata)

# CD4/8 classification
from sklearn.neighbors import NearestNeighbors

nn=20

nbrs = NearestNeighbors(n_neighbors=nn, algorithm='ball_tree').fit(adata.obsm['X_pca'])
distances, indices = nbrs.kneighbors(adata.obsm['X_pca'])
markers = ['CD4','CD8A','CD3D']

for gene in markers:
    adata.obs['knn'+gene] = np.array(adata.raw[:,gene].X[indices.reshape(-1)].reshape(-1,nn).mean(axis=1)).reshape(-1)

cd4_thresh = .2
cd8_thresh=1
cd3_thresh=1

adata.obs['subtype'] = 'Unknown'
idx = (adata.obs['knnCD4']>cd4_thresh)&(adata.obs['knnCD8A']<cd8_thresh)&(adata.obs['knnCD3D']>cd3_thresh)
adata.obs.loc[idx,'subtype']='CD4 T'

idx = (adata.obs['knnCD4']<cd4_thresh)&(adata.obs['knnCD8A']>cd8_thresh)&(adata.obs['knnCD3D']>cd3_thresh)
adata.obs.loc[idx,'subtype']='CD8 T'

adata.obs['CAR'] = np.array((adata.raw[:,'Kymriah'].X>0).todense()).reshape(-1)

# Clustering
sc.tl.leiden(adata,resolution=0.3,key_added='leiden_0.3')
sc.tl.leiden(adata,resolution=1.2,key_added='leiden_1.2')
sc.tl.leiden(adata,resolution=3.0,key_added='leiden_3.0')

adata.write("../data/retreatment_patient_unnamed.h5ad")

names = ['CD4 T IP','CD8 T','CD8 T','CD4 T','Myeloid','$\gamma\delta$ T','Ery','Platelet']
adata.obs['cell_type'] = [names[int(i)] for i in adata.obs['leiden_0.3']]

adata.obs.loc[adata.obs['leiden_3.0']=='34','cell_type'] = 'NK'
adata.obs.loc[adata.obs['leiden_3.0']=='40','cell_type'] = 'T-M doublets'
adata.obs.loc[adata.obs['leiden_1.2']=='17','cell_type'] = 'CD8 T IP'

#adata = adata[~adata.obs['cell_type'].str.match('.*doublet')]

adata.write("../data/retreatment_patient_clustering.h5ad")
adata.obs.to_csv("../data/retreatment_patient_clustering_obs.csv")


