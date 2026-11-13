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