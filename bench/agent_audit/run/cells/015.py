import matplotlib;matplotlib.use('Agg')
sc.pl.umap(a,color=['cell_type','label','ISG15'],show=False,wspace=0.5);import matplotlib.pyplot as plt;plt.savefig('umap.png',dpi=120,bbox_inches='tight')
pd.concat({c:r for c,r in res.items()},names=['cell_type','gene']).to_csv('pseudobulk_DE_stim_vs_ctrl.csv')
ab.to_csv('abundance_stim_vs_ctrl.csv',index=False)
del a.uns['rank_genes_groups']
a.write_h5ad('annotated.h5ad');print(a.obs.cell_type.value_counts())