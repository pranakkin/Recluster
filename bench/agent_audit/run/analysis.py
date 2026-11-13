# --- cell 001 ---
import scanpy as sc, pandas as pd, numpy as np
a=sc.read_h5ad('pbmc_ifnb.h5ad');print(a);print(a.obs.head());print(pd.crosstab(a.obs.replicate,a.obs.label));print(a.X[:5,:20].toarray() if hasattr(a.X,'toarray') else a.X[:5,:20])

# --- cell 002 ---
a.var_names_make_unique()
a.var['mt']=a.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],inplace=True,percent_top=None)
print(a.obs[['n_genes_by_counts','total_counts','pct_counts_mt']].describe())
print(a.var.mt.sum())
import importlib;print([m for m in ['harmonypy','scvi','bbknn'] if importlib.util.find_spec(m)])

# --- cell 003 ---
print((a.obs.n_genes_by_counts<200).sum(), (a.obs.n_genes_by_counts>2000).sum())
a=a[a.obs.n_genes_by_counts>=200].copy()
sc.pp.filter_genes(a,min_cells=10)
a.layers['counts']=a.X.copy()
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000,batch_key='label')
sc.pp.pca(a,n_comps=30,mask_var='highly_variable')
sc.external.pp.harmony_integrate(a,'label')
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,key_added='leiden',flavor='igraph',n_iterations=2)
print(a.shape, pd.crosstab(a.obs.leiden,a.obs.label))

# --- cell 004 ---
import harmonypy
ho=harmonypy.run_harmony(a.obsm['X_pca'],a.obs,'label',max_iter_harmony=20)
Z=np.asarray(ho.Z_corr); Z=Z if Z.shape[0]==a.n_obs else Z.T
a.obsm['X_pca_harmony']=Z
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))

# --- cell 005 ---
import scipy, sklearn;print(scipy.__version__, sc.__version__, sklearn.__version__)
a.obsm['X_pca_harmony']=np.ascontiguousarray(Z,dtype=np.float32)
try:
    sc.pp.neighbors(a,use_rep='X_pca_harmony')
except Exception as e: print(e)

# --- cell 006 ---
import scipy.sparse as sp
if not hasattr(sp.csr_array,'getnnz'):
    sp.csr_array.getnnz=lambda self,axis=None: self.count_nonzero(axis=axis) if axis is not None else self.nnz
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))

# --- cell 007 ---
def _getnnz(self,axis=None):
    if axis is None: return self.nnz
    if axis==1: return np.diff(self.tocsr().indptr)
    return np.bincount(self.tocsr().indices,minlength=self.shape[1])
sp.csr_array.getnnz=_getnnz
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))

# --- cell 008 ---
import scanpy.neighbors._common as nc, scanpy.neighbors as snb
_orig=nc._get_indices_distances_from_sparse_matrix
def _patched(*args,**kw):
    i,d=_orig(*args,**kw); return np.ascontiguousarray(i),np.ascontiguousarray(d)
nc._get_indices_distances_from_sparse_matrix=_patched
snb._get_indices_distances_from_sparse_matrix=_patched
sc.pp.neighbors(a,use_rep='X_pca_harmony');sc.tl.umap(a);sc.tl.leiden(a,resolution=0.8,flavor='igraph',n_iterations=2)
print(pd.crosstab(a.obs.leiden,a.obs.label))

# --- cell 009 ---
mk=['CD3E','CD4','IL7R','CCR7','SELL','CD8A','GZMK','GNLY','NKG7','GZMB','MS4A1','CD79A','CD14','LYZ','S100A8','FCGR3A','MS4A7','FCER1A','CLEC10A','LILRA4','IL3RA','PPBP','PF4','HBB','CD34','FOXP3','MKI67','ISG15']
mk=[g for g in mk if g in a.var_names];print(set(['CD3E','CD4','IL7R','CCR7','SELL','CD8A','GZMK','GNLY','NKG7','MS4A1','CD14','LYZ','FCGR3A','FCER1A','LILRA4','PPBP','HBB','FOXP3','MKI67'])-set(mk))
df=sc.get.obs_df(a,mk+['leiden']).groupby('leiden').mean().round(2);print(df.to_string())
sc.tl.rank_genes_groups(a,'leiden',method='wilcoxon')
print(pd.DataFrame(a.uns['rank_genes_groups']['names']).head(8).to_string())

# --- cell 010 ---
g=['CD3D','CD8B','CD4','IL7R','CCR7','HSPA1A','DNAJB1','CACYBP','FOS','JUN','S100A4','CREM','GZMA','CCL5','TRDC','KLRB1','MALAT1']
g=[x for x in g if x in a.var_names]
print(sc.get.obs_df(a,g+['leiden']).groupby('leiden',observed=True).mean().round(2).to_string())
print(a.obs.groupby('leiden',observed=True)[['n_genes_by_counts','total_counts']].median())
print(pd.crosstab(a.obs.leiden,a.obs.replicate))

# --- cell 011 ---
ann={'0':'CD14 Mono','1':'CD16 Mono','2':'T (stress/HSP-high)','3':'CD4 T memory','4':'CD4 T naive','6':'CD4 T naive','5':'pDC','7':'NK','8':'CD8 T','9':'B','10':'cDC','11':'Megakaryocyte/platelet'}
a.obs['cell_type']=a.obs.leiden.map(ann).astype('category')
print(importlib.util.find_spec('pydeseq2'))
ct=pd.crosstab([a.obs.replicate,a.obs.label],a.obs.cell_type)
prop=ct.div(ct.sum(1),axis=0)
from scipy.stats import wilcoxon
rows=[]
for c in ct.columns:
    p=prop[c].unstack('label')
    lfc=np.log2((p.stim+1e-3)/(p.ctrl+1e-3))
    rows.append([c,ct[c].xs('ctrl',level=1).sum(),ct[c].xs('stim',level=1).sum(),p.ctrl.mean()*100,p.stim.mean()*100,lfc.median(),(lfc>0).sum(),wilcoxon(p.stim,p.ctrl).pvalue])
ab=pd.DataFrame(rows,columns=['cell_type','n_ctrl','n_stim','mean%ctrl','mean%stim','median_log2FC','donors_up','wilcoxon_p'])
from statsmodels.stats.multitest import multipletests
ab['padj']=multipletests(ab.wilcoxon_p,method='fdr_bh')[1];print(ab.round(3).to_string())

# --- cell 012 ---
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import warnings,logging;warnings.filterwarnings('ignore')
inf=DefaultInference(n_cpus=4)
res={}
for c in a.obs.cell_type.cat.categories:
    s=a[a.obs.cell_type==c]
    keys=s.obs.replicate.astype(str)+'|'+s.obs.label.astype(str)
    vc=keys.value_counts(); ok=vc[vc>=10].index
    donors=[d for d in s.obs.replicate.unique() if f'{d}|ctrl' in ok and f'{d}|stim' in ok]
    if len(donors)<3: print(c,'skipped, donors',len(donors));continue
    samples=[f'{d}|{l}' for d in donors for l in ['ctrl','stim']]
    M=np.vstack([np.asarray(s.layers['counts'][(keys==k).values].sum(0)).ravel() for k in samples])
    cnt=pd.DataFrame(M,index=samples,columns=s.var_names).astype(int)
    cnt=cnt.loc[:,(cnt>=10).sum(0)>=len(donors)//2+1]
    meta=pd.DataFrame({'donor':[x.split('|')[0] for x in samples],'condition':[x.split('|')[1] for x in samples]},index=samples)
    dds=DeseqDataSet(counts=cnt,metadata=meta,design='~donor+condition',inference=inf,quiet=True);dds.deseq2()
    st=DeseqStats(dds,contrast=['condition','stim','ctrl'],inference=inf,quiet=True);st.summary()
    r=st.results_df.dropna(subset=['padj']);r['n_donors']=len(donors);res[c]=r
    print(c,len(donors),'tested',len(r),'up',((r.padj<0.05)&(r.log2FoldChange>1)).sum(),'down',((r.padj<0.05)&(r.log2FoldChange<-1)).sum())

# --- cell 013 ---
for c,r in res.items():
    sig=r[r.padj<0.05]
    print(f'\n{c}: UP',', '.join(f'{g}({v:.1f})' for g,v in sig.log2FoldChange.nlargest(10).items()))
    print('   DOWN',', '.join(f'{g}({v:.1f})' for g,v in sig.log2FoldChange.nsmallest(8).items()))
core=['ISG15','IFI6','IFIT1','IFIT3','MX1','LY6E','CXCL10','CCL8','CCL2','IFITM3','APOBEC3A','OAS1','IFI27']
tab=pd.DataFrame({c:r.log2FoldChange.reindex(core) for c,r in res.items()}).round(1);print(tab.to_string())

# --- cell 014 ---
st=a[a.obs.label=='stim']
X=sc.get.obs_df(st,['CXCL10','CCL8','ISG15','IFIT1','CD14','LYZ']+['cell_type'])
print(X.groupby('cell_type',observed=True).agg(lambda v:(v>0).mean()).round(2).to_string())
up={c:set(r.index[(r.padj<0.05)&(r.log2FoldChange>1)]) for c,r in res.items()}
shared=set.intersection(*up.values());print(len(shared),sorted(shared)[:40])
mono_only=(up['CD14 Mono']|up['CD16 Mono']|up['cDC'])-set.union(*[up[k] for k in up if k not in['CD14 Mono','CD16 Mono','cDC']]);print(len(mono_only))

# --- cell 015 ---
import matplotlib;matplotlib.use('Agg')
sc.pl.umap(a,color=['cell_type','label','ISG15'],show=False,wspace=0.5);import matplotlib.pyplot as plt;plt.savefig('umap.png',dpi=120,bbox_inches='tight')
pd.concat({c:r for c,r in res.items()},names=['cell_type','gene']).to_csv('pseudobulk_DE_stim_vs_ctrl.csv')
ab.to_csv('abundance_stim_vs_ctrl.csv',index=False)
del a.uns['rank_genes_groups']
a.write_h5ad('annotated.h5ad');print(a.obs.cell_type.value_counts())
