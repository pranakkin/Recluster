mk=['CD3E','CD4','IL7R','CCR7','SELL','CD8A','GZMK','GNLY','NKG7','GZMB','MS4A1','CD79A','CD14','LYZ','S100A8','FCGR3A','MS4A7','FCER1A','CLEC10A','LILRA4','IL3RA','PPBP','PF4','HBB','CD34','FOXP3','MKI67','ISG15']
mk=[g for g in mk if g in a.var_names];print(set(['CD3E','CD4','IL7R','CCR7','SELL','CD8A','GZMK','GNLY','NKG7','MS4A1','CD14','LYZ','FCGR3A','FCER1A','LILRA4','PPBP','HBB','FOXP3','MKI67'])-set(mk))
df=sc.get.obs_df(a,mk+['leiden']).groupby('leiden').mean().round(2);print(df.to_string())
sc.tl.rank_genes_groups(a,'leiden',method='wilcoxon')
print(pd.DataFrame(a.uns['rank_genes_groups']['names']).head(8).to_string())