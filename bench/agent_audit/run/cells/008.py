import scanpy.neighbors._common as nc, scanpy.neighbors as snb
_orig=nc._get_indices_distances_from_sparse_matrix
def _patched(*args,**kw):
    i,d=_orig(*args,**kw); return np.ascontiguousarray(i),np.ascontiguousarray(d)
nc._get_indices_distances_from_sparse_matrix=_patched
snb._get_indices_distances_from_sparse_matrix=_patched
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))