import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Calculates pseudobulk count matrix
# min_cells: Minimum number of cells for a sample to be included
# cpm_thresh: Minimum mean counts-per-million for a gene to be included
# min_pat: Require at least this many patients ot express the gene with at least min_pat_cpm expression to include
def get_pseudobulk(adata,id_vars,clin_vars,min_cells=25,cpm_thresh=0,min_pat=5,min_pat_cpm=10):
    
    adata.var['cpm'] = np.maximum(.1,1e6*np.array(adata.X.sum(axis=0)).reshape(-1)/adata.X.sum())

    gene_use = adata.var.index[adata.var['cpm']>=cpm_thresh]
    
    N = adata.obs.value_counts(id_vars)
    
    Xsum = pd.DataFrame(index=N.index,columns=gene_use)
    for ind,g in adata.obs.groupby(id_vars):
        Xsum.loc[ind,:] = np.array(adata[g.index,gene_use].X.sum(axis=0)).reshape(-1)
    
    clin = adata.obs[clin_vars + id_vars].drop_duplicates().set_index(id_vars)

    Xsum = clin.join(Xsum)
    
    # Apply min cells filter
    use_idx = N>=min_cells
    Xsum=Xsum[use_idx].reset_index()

    # Apply min patient filter
    genes = Xsum.columns[(len(id_vars)+len(clin_vars)):]
    cpm = 1e6*(Xsum.loc[:,genes].T / Xsum.loc[:,genes].T.sum()).T
    nexp = (cpm>min_pat_cpm).sum()
    Xsum = Xsum[list(id_vars)+list(clin_vars)+list(genes[nexp>=min_pat])]

    return(Xsum)


adata = sc.read("../data/CART_raw_data.h5ad")


# Identify gender (note we use transcriptomics rather than clinical data since one patient received an allogenic stem cell transplant)
adata.obs['eX'] = np.array(adata[:,'XIST'].X.todense()).reshape(-1)
adata.obs['eY'] = np.array(adata[:,'RPS4Y1'].X.todense()).reshape(-1)
S = adata.obs.groupby('barcode')[['eX','eY']].sum()
S['gender'] = 'female'
S.loc[(S['eY']+1)/(S['eX']+1) > .5,'gender'] = 'male'
adata.obs=adata.obs.join(S[['gender']],on='barcode',how='left')


adata.obs['subtype'] = pd.read_csv('../data/ALLT_subtypes_annotated_obs.txt')
adata.obs['type']=adata.obs['timepoint'].astype(str)
adata.obs.loc[adata.obs['CAR']&adata.obs['timepoint'].str.match('D7'),'type'] = 'D7-CAR-T'
adata.obs.loc[~adata.obs['CAR']&(adata.obs['timepoint']=="D7-CART"),'type'] = 'Unknown'

# Filter to T-cells
keep_idx = (adata.obs['timepoint']=='Infusion')|adata.obs['subtype'].str.match('.*T')
adata = adata[keep_idx]

id_vars = ['barcode','type','subtype']
clin_vars = ['response','product','gender','batch']

# Just T-cells, with known CAR status
idx = ~adata.obs['subtype'].isna()&~(adata.obs['type']=='Unknown')
Xsum = get_pseudobulk(adata[idx],id_vars,clin_vars)

Xsum.to_csv('../data/ALL_T_cells_pseudobulk.v2.csv',index=False)

# Make version also that differentiaties CAR vs nonCAR in the IPs
id_vars = ['barcode','type']
clin_vars = ['response','product','gender','batch']
idx = ~(adata.obs['type']=='Unknown')
Xsum = get_pseudobulk(adata[idx],id_vars,clin_vars,cpm_thresh=0)

Xsum.to_csv('../data/ALL_T_cells_pseudobulk_nosubtypesplit.v3.csv',index=False)

