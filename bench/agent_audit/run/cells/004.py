import harmonypy
ho=harmonypy.run_harmony(a.obsm['X_pca'],a.obs,'label',max_iter_harmony=20)
Z=np.asarray(ho.Z_corr); Z=Z if Z.shape[0]==a.n_obs else Z.T
a.obsm['X_pca_harmony']=Z
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))