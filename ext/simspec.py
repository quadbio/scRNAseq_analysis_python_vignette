import numpy as np
import pandas as pd
import anndata
import scanpy as sc

from scipy import sparse
from scipy.stats import rankdata
from sklearn.preprocessing import scale as scaledata
from sklearn.decomposition import PCA
from tqdm import tqdm

# Calculate correlation between columns in matrix X and those in matrix Y, with X and Y being dense or sparse
def corSparse(X, Y = None):
  if Y is None: Y = X
  
  n = X.shape[0]
  muX = np.ravel(X.mean(0))
  muY = np.ravel(Y.mean(0))
  covmat = ( X.T.dot(Y) - (n * muX[:,np.newaxis].dot(muY[:,np.newaxis].T)) ) / (n-1)
  sdvecX = np.ravel(np.sqrt(((X.power(2)).sum(0) - n*(muX**2)) / (n-1)) if sparse.issparse(X) else np.sqrt(((X**2).sum(0) - n*(muX**2)) / (n-1)))
  sdvecY = np.ravel(np.sqrt(((Y.power(2)).sum(0) - n*(muY**2)) / (n-1)) if sparse.issparse(Y) else np.sqrt(((Y**2).sum(0) - n*(muY**2)) / (n-1)))
  cormat = covmat / sdvecX[:,np.newaxis].dot(sdvecY[:,np.newaxis].T)
  
  return np.array(cormat)

# For each column of the input matrix X, convert the values into the ranks (for calculation of Spearman correlation, for example)
def rankMatrix(X):
  if sparse.issparse(X):
    idx_row, idx_col, dat = sparse.find(X)
    df = pd.DataFrame({'i' : idx_row, 'j' : idx_col, 'x' : dat}).sort_values('j')
    df['r'] = np.concatenate([ rankdata(x) + (df.shape[0]-len(x)) - (1+(df.shape[0]-len(x)))/2 for x in np.split(df.x, np.unique(df.j, return_index=True)[1][1:]) ])
    ranked = sparse.csr_matrix((df['r'], (df['i'], df['j'])), shape = X.shape)
  else:
    ranked = rankdata(X, method = "average", axis=0)
  
  return ranked

# Similar to rankMatrix but to use "dense" method for ranking non-zero values
def rankMatrix_dense(X):
  if sparse.issparse(X):
    idx_row, idx_col, dat = sparse.find(X)
    df = pd.DataFrame({'i' : idx_row, 'j' : idx_col, 'x' : dat}).sort_values('j')
    df['r'] = np.concatenate([ rankdata(x, method = "dense") for x in np.split(df['x'], np.unique(df['j'], return_index=True)[1][1:]) ])
    ranked = sparse.csr_matrix((df['r'], (df['i'], df['j'])), shape = X.shape)
  else:
    ranked = rankdata(X, method = "average", axis=0)
  
  return ranked

# Similar to rankMatrix, but only ranking non-zero values with zero values remaining zero
def rankMatrix_nonzero(X):
  if sparse.issparse(X):
    idx_row, idx_col, dat = sparse.find(X)
    df = pd.DataFrame({'i' : idx_row, 'j' : idx_col, 'x' : dat}).sort_values('j')
    value_split = np.split(df.x, np.unique(df.j, return_index=True)[1][1:])
    df['r'] = np.concatenate([ rankdata(x) for x in value_split ])
    ranked = sparse.csr_matrix((df['r'], (df['i'], df['j'])), shape = X.shape)
  else:
    ranked = rankdata(X, method = "average", axis=0)
  
  return ranked

## summarize single-cell data to pseudocell-level
def group_vec_to_ident_mat(group, norm = True):
  i = np.where(group.notnull())[0]
  j = group[i].astype(int)
  mat_ident = sparse.csr_matrix((np.repeat(1, len(i)), (i, j)), shape=(len(group), len(np.unique(j))))
  if norm:
      num_cells_per_group = np.ravel(mat_ident.sum(axis=0))
      mat_ident = mat_ident @ sparse.diags(1 / num_cells_per_group)
  
  return mat_ident

def summarize_numeric_matrix(mat, group, use_mean = True): # mat.shape = cell*feature; len(group) = cell
  mat_ident = group_vec_to_ident_mat(group, norm = use_mean)
  mat_summ = mat_ident.transpose() @ mat
  return mat_summ

# Reference Similarity Spectrum (RSS)
def ref_sim_spectrum(adata, 
                     ref, # rows are samples and cols are features
                     vars_ref = None,
                     layer=None,
                     method='pearson',
                     scale=True,
                     return_sim_only=False,
                     copy=False
                    ):
    if isinstance(ref, pd.DataFrame):
        ref = anndata.AnnData(X = np.array(ref),
                              var = pd.DataFrame(index = ref.columns if vars_ref is None else vars_ref),
                              obs = pd.DataFrame(index = ref.index))
    elif not isinstance(ref, anndata.AnnData):
        ref = anndata.AnnData(X = ref,
                              var = pd.DataFrame(index = vars_ref),
                              obs = np.arange(ref.shape[0]))
    
    shared_genes = np.intersect1d(adata.var_names, ref.var_names)
    X = adata[:,shared_genes].X.T if layer is None or layer not in adata.layers.keys() else adata[:,shared_genes].layers[layer].T
    refX = ref[:,shared_genes].X.T
    
    if method == 'spearman':
        X = rankMatrix_nonzero(X)
        refX = rankMatrix_nonzero(refX)
    corr = corSparse(X, refX)
    if scale:
        corr = scaledata(corr, axis = 1)
    corr[np.isnan(corr)] = 0
    
    if return_sim_only:
        return corr
    
    rep_name = ('X' if layer is None or layer not in adata.layers.keys() else layer)+'_rss'
    if copy:
        adata = adata.copy()
        adata.obsm[rep_name] = corr
        return adata
    else:
        adata.obsm[rep_name] = corr
        return

# Cluster Similarity Spectrum (CSS)
def cluster_sim_spectrum(adata,
                         batch = 'batch',
                         use_rep='X_pca',
                         layer=None,
                         n_neighbors=15,
                         n_pcs=None,
                         method_clustering='louvain',
                         resolution_clustering=1,
                         highly_variable=True,
                         method_corr='pearson',
                         verbose=True,
                         return_corr_only=False,
                         copy=False
                        ):
    if copy:
        adata_raw = adata.copy()
    else:
        adata_raw = adata
    if highly_variable and 'highly_variable' in adata.var.columns:
        adata = adata[:,adata.var['highly_variable']]
    
    if verbose:
        print('Splitting samples...')
    adata_batch = [ adata[adata.obs[batch] == x,:].copy() for x in adata.obs[batch].unique() ]
    
    if verbose:
        print('Clustering each sample...')
    for ad in tqdm(adata_batch):
        sc.pp.neighbors(ad, use_rep=use_rep, n_neighbors=n_neighbors, n_pcs=n_pcs)
        if method_clustering=='louvain':
            sc.tl.louvain(ad, resolution=resolution_clustering, key_added='cluster')
        if method_clustering=='leiden':
            sc.tl.leiden(ad, resolution=resolution_clustering, key_added='cluster')

    if verbose:
        print('Estimate average transcriptomic profiles for clusters...')
    avg_expr_cl = list()
    for ad in tqdm(adata_batch):
        avg_expr = summarize_numeric_matrix(ad.X if layer is None or layer not in ad.layers.keys() else ad.layers[layer], ad.obs['cluster'])
        avg_expr_cl.append(avg_expr)

    if verbose:
        print('Calculate similarities...')
    sims = list()
    for avg_expr in tqdm(avg_expr_cl):
        sim = ref_sim_spectrum(adata, ref=avg_expr, vars_ref=adata.var_names, layer=layer, method=method_corr, scale=True, return_sim_only=True)
        sims.append(sim)

    if verbose:
        print('Generate final embedding...')
    css = np.concatenate(sims, axis=1)

    if return_corr_only:
        return css
    rep_name = ('X' if layer is None or layer not in adata.layers.keys() else layer)+'_css'
    adata_raw.obsm[rep_name] = css.copy()
    if copy:
        return adata_raw
    else:
        return

# PCA on the embeddings
def run_PCA(adata,
            use_rep,
            n_pcs = 20,
            copy=False
           ):
    embed = adata.obsm[use_rep].copy()
    embed_pca = PCA(n_components=n_pcs).fit_transform(embed)
    rep_name = use_rep + '_pca'
    if copy:
        adata = adata.copy()
        adata.obsm[rep_name] = embed_pca
        return adata
    else:
        adata.obsm[rep_name] = embed_pca
        return
