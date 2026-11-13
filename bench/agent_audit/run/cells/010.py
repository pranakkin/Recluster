g=['CD3D','CD8B','CD4','IL7R','CCR7','HSPA1A','DNAJB1','CACYBP','FOS','JUN','S100A4','CREM','GZMA','CCL5','TRDC','KLRB1','MALAT1']
g=[x for x in g if x in a.var_names]
print(sc.get.obs_df(a,g+['leiden']).groupby('leiden',observed=True).mean().round(2).to_string())
print(a.obs.groupby('leiden',observed=True)[['n_genes_by_counts','total_counts']].median())
print(pd.crosstab(a.obs.leiden,a.obs.replicate))