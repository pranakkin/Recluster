def _getnnz(self,axis=None):
    if axis is None: return self.nnz
    if axis==1: return np.diff(self.tocsr().indptr)
    return np.bincount(self.tocsr().indices,minlength=self.shape[1])
sp.csr_array.getnnz=_getnnz
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))