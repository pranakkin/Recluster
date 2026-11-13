st=a[a.obs.label=='stim']
X=sc.get.obs_df(st,['CXCL10','CCL8','ISG15','IFIT1','CD14','LYZ']+['cell_type'])
print(X.groupby('cell_type',observed=True).agg(lambda v:(v>0).mean()).round(2).to_string())
up={c:set(r.index[(r.padj<0.05)&(r.log2FoldChange>1)]) for c,r in res.items()}
shared=set.intersection(*up.values());print(len(shared),sorted(shared)[:40])
mono_only=(up['CD14 Mono']|up['CD16 Mono']|up['cDC'])-set.union(*[up[k] for k in up if k not in['CD14 Mono','CD16 Mono','cDC']]);print(len(mono_only))