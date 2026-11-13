import scipy.sparse as sp
if not hasattr(sp.csr_array,'getnnz'):
    sp.csr_array.getnnz=lambda self,axis=None: self.count_nonzero(axis=axis) if axis is not None else self.nnz
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))