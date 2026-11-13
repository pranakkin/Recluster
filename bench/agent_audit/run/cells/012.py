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