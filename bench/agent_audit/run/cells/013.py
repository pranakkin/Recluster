for c,r in res.items():
    sig=r[r.padj<0.05]
    print(f'\n{c}: UP',', '.join(f'{g}({v:.1f})' for g,v in sig.log2FoldChange.nlargest(10).items()))
    print('   DOWN',', '.join(f'{g}({v:.1f})' for g,v in sig.log2FoldChange.nsmallest(8).items()))
core=['ISG15','IFI6','IFIT1','IFIT3','MX1','LY6E','CXCL10','CCL8','CCL2','IFITM3','APOBEC3A','OAS1','IFI27']
tab=pd.DataFrame({c:r.log2FoldChange.reindex(core) for c,r in res.items()}).round(1);print(tab.to_string())