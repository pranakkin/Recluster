import scipy, sklearn;print(scipy.__version__, sc.__version__, sklearn.__version__)
a.obsm['X_pca_harmony']=np.ascontiguousarray(Z,dtype=np.float32)
try:
    sc.pp.neighbors(a,use_rep='X_pca_harmony')
except Exception as e: print(e)