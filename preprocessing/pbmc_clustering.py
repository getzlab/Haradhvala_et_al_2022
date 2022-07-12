import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scrublet as scr

input_h5ad = "../data/CART_raw_data.h5ad"

adata = sc.read(input_h5ad)
adata = adata[adata.obs['timepoint']!='Infusion']

scrubs = dict()
for channel,g in adata.obs.groupby('channel'):
    print(channel)
    
    scrub = scr.Scrublet(adata[g.index,:].X)
    adata.obs.loc[g.index,'doublet_scores'], adata.obs.loc[g.index,'predicted_doublets'] = scrub.scrub_doublets()
    
    scrub.plot_histogram()
    
    scrubs[channel] = scrub
    
#adata.obs['predicted_doublets'] = adata.obs['predicted_doublets']=="True"

## Overide cases where scrublet picked a bad cutoff
channel = 'A4_2'
idx = adata.obs['channel']==channel
adata.obs.loc[idx,'predicted_doublets'] = scrubs[channel].call_doublets(threshold=.45)

channel = 'F7'
idx = adata.obs['channel']==channel
adata.obs.loc[idx,'predicted_doublets'] = scrubs[channel].call_doublets(threshold=.99)

channel = 'H1'
idx = adata.obs['channel']==channel
adata.obs.loc[idx,'predicted_doublets'] = scrubs[channel].call_doublets(threshold=.6)

channel = 'Pilot_A1'
idx = adata.obs['channel']==channel
adata.obs.loc[idx,'predicted_doublets'] = scrubs[channel].call_doublets(threshold=.3)

channel = 'Pilot_A2'
idx = adata.obs['channel']==channel
adata.obs.loc[idx,'predicted_doublets'] = scrubs[channel].call_doublets(threshold=.35)

# Since we only have 4, D14 samples not included in primary analysis
idx = adata.obs['timepoint'] == 'D14'
adata = adata[~idx]

## Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 15, :]

adata = adata[adata.obs['predicted_doublets']==False]

adata.obs = adata.obs.drop(columns=['predicted_doublets'])

adata.write("../data/PBMC_filtered.h5ad")

## Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


## Highly variable selection
sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=3, min_disp=0.2)

# Blacklist genes with sex bias, CAR transcripts, and TCR/BCR variable regions from clustering analyses
adata.var.loc[['XIST','RPS4Y1','RPS4Y2','Yescarta','Kymriah'],'highly_variable'] = False
adata.var.loc[adata.var.index.str.match('TR.V|IG.V'),'highly_variable'] = False

adata.raw = adata

adata = adata[:, adata.var.highly_variable]


sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

sc.external.pp.harmony_integrate(adata,'channel')
sc.pp.neighbors(adata,use_rep='X_pca_harmony')
sc.tl.umap(adata)

sc.tl.leiden(adata,resolution=2.5,key_added='leiden_2.5')

# Filter a few small doublet clusters with markers from multiple cell types
doublet_clusters = ['30','29','35','20']
adata = adata[~adata.obs['leiden_2.5'].isin(doublet_clusters)]

idx = np.array(adata.raw[:,'HBB'].X.todense()>5).reshape(-1)
adata = adata[~idx]
adata.raw[:,'HBB'].X[:] = 0

## Label cell types
markers = ['CD3D','NKG7','LYZ','MS4A1','PF4','LILRA4','CD1C']
X = pd.DataFrame(adata.raw[:,markers].X.todense(),index=adata.obs.index,columns=markers)

X['leiden'] = adata.obs['leiden_2.5']
Xm=X.groupby('leiden').mean()
Xm['cell_type'] = ""
Xm.loc[Xm['LILRA4']>1,'cell_type'] = 'pDC'
Xm.loc[Xm['PF4']>2,'cell_type'] = 'Platelet'

Xm.loc[Xm['NKG7']>2,'cell_type'] = 'NK'
Xm.loc[Xm['CD3D']>.75,'cell_type'] = 'T'
Xm.loc[Xm['LYZ']>2,'cell_type'] = 'Myeloid'
Xm.loc[Xm['CD1C']>1,'cell_type'] = 'mDC'
Xm.loc[Xm['MS4A1']>1,'cell_type'] = 'B'

adata.obs['cell_type'] = adata.obs['leiden_2.5'].map(Xm['cell_type'])

adata.write("../data/PBMC_specific_clustering.h5ad")
adata.obs.to_csv("../data/PBMC_specific_clustering_obs.csv")



