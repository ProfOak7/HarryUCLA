from pathlib import Path
import csv, re, math
import numpy as np
from scipy.stats import ttest_ind
ROOT=Path(__file__).resolve().parent

def read(name,delimiter=','):
    with (ROOT/name).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter=delimiter))
def number(x): return float(x) if x not in ('',None) else np.nan
def label(x):
    s=str(x or '').strip()
    return '' if not s or re.fullmatch(r'\d+(\.0+)?|nan|none',s,re.I) else s

def load(tissue):
    rows=read('BALF_Normalized_named_reconstructed.csv' if tissue=='BALF' else 'Annotation.csv')
    names=[k for k in rows[0] if k.startswith('Grp' if tissue=='BALF' else 'TMT_Grp')]
    ids=[r['Master Protein'] for r in rows]
    groups=np.array([int(re.search(r'Grp([1-4])',k).group(1)) for k in names])
    data=np.array([[number(r[k]) for k in names] for r in rows])
    assert len(set(ids))==len(ids) and len(names)==30
    return ids,groups,data,rows

COMPS=[('A',(2,),(1,),'Air + Ex vs Air + Rest','grp21pi'),('B',(3,),(1,),'CS + Rest vs Air + Rest','grp31pi'),('C',(4,),(3,),'CS + Ex vs CS + Rest','grp43pi'),('D',(4,),(2,),'CS + Ex vs Air + Ex','grp42pi'),('CS',(3,4),(1,2),'CS vs Air (pooled)','grp12vs34pi'),('Ex',(2,4),(1,3),'Exercise vs Rest (pooled)','grp13vs24pi')]
HEADERS=['Master protein accession','Archived protein label','Annotation status','n numerator','n denominator','Mean numerator','Mean denominator','Fold change','Log2 fold change','Raw p-value','Signed pi-score','BH q-value','Bonferroni p-value','Volcano cutoff','|pi| > 1.0','|pi| > 1.3','BH q < 0.05']

def bh(p):
    idx=np.flatnonzero(np.isfinite(p)); order=idx[np.argsort(p[idx])]; q=np.full_like(p,np.nan); q[order]=np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1].clip(0,1); return q

def build():
    tables={}; summaries=[]; checks=[]
    for tissue in ['BALF','Serum']:
        ids,g,m,rows=load(tissue)
        ann={}; legacy={}
        if tissue=='Serum': ann={r['Master Protein']:label(r['Proteins']) for r in rows}
        for code,a,b,title,stem in COMPS:
            if tissue=='BALF' and code in ['CS','Ex']:continue
            old=read(stem+(' (2)' if tissue=='BALF' else '')+'.csv'); legacy[code]=old
            for r in old:
                k=r['Master Protein']; v=label(r['Proteins'])
                if v:
                    assert not ann.get(k) or ann[k]==v,(k,ann.get(k),v)
                    ann[k]=v
        print(tissue, 'annotations',sum(bool(ann.get(k)) for k in ids),'of',len(ids))
        for code,a,b,title,stem in COMPS:
            x=m[:,np.isin(g,a)]; y=m[:,np.isin(g,b)]
            nx=np.isfinite(x).sum(1); ny=np.isfinite(y).sum(1); mx=np.nanmean(x,axis=1); my=np.nanmean(y,axis=1); fc=mx/my; lfc=np.log2(fc)
            p=ttest_ind(x,y,axis=1,equal_var=False,nan_policy='omit').pvalue; pi=lfc*-np.log10(p); q=bh(p); bon=np.minimum(1,p*np.isfinite(p).sum())
            sel=np.abs(pi)>1; strict=np.abs(pi)>1.3; vol=(p<.05)&(np.abs(lfc)>=.5); fdr=q<.05
            assert np.isfinite(p).all()
            direction=np.where(vol,np.where(lfc>0,'Higher','Lower'),'Below cutoff')
            name=tissue+'_'+code
            if code in legacy:
                old=legacy[code]; index={k:i for i,k in enumerate(ids)}; ix=[index[r['Master Protein']] for r in old]
                assert {k for k,s in zip(ids,sel) if s}=={r['Master Protein'] for r in old}
                pi_col=next(k for k in old[0] if k.endswith('pi score'))
                for arr,col in [(p,'ttest_pvalue'),(lfc,'Log2(Fold Change)'),(pi,pi_col)]:np.testing.assert_allclose(arr[ix],[float(r[col]) for r in old],rtol=1e-6,atol=1e-8)
                valid='Exact archived hit set; statistics matched'
            else:valid='Recomputed pooled contrast; no archived list available'
            # Verify against the exact revised Figure 3 full tables, not just selected rows.
            if tissue=='BALF' and code in 'ABCD':
                path=ROOT/'BALF_Figure3_revision'/'results'/f'{code}_Grp{a[0]}_vs_Grp{b[0]}_all_proteins.csv'
                if path.exists():
                    ref=read(str(path.relative_to(ROOT))); assert [r['Master Protein'] for r in ref]==ids
                    for arr,col in [(p,'raw p'),(lfc,'log2 fold change'),(pi,'pi score')]:np.testing.assert_allclose(arr,[float(r[col]) for r in ref],rtol=1e-9,atol=1e-12)
            order=sorted(range(len(ids)),key=lambda i:(not sel[i],-abs(pi[i]),ids[i]))
            def yn(x):return 'Yes' if x else 'No'
            out=[]
            for i in order:
                k=ids[i]; v=ann.get(k,''); out.append([k,v or None,'Archived label' if v else 'Accession only',int(nx[i]),int(ny[i]),float(mx[i]),float(my[i]),float(fc[i]),float(lfc[i]),float(p[i]),float(pi[i]),float(q[i]),float(bon[i]),str(direction[i]),yn(sel[i]),yn(strict[i]),yn(fdr[i])])
            tables[name]=out
            summaries.append([name,tissue,title,'+'.join(map(str,a)),'+'.join(map(str,b)),int(x.shape[1]),int(y.shape[1]),len(ids),int(vol.sum()),int(sel.sum()),int(strict.sum()),int(fdr.sum()),valid])
            print(summaries[-1])
    return tables,summaries
if __name__=='__main__':
    tables,summaries=build()
    out=ROOT/'supplement_csv';out.mkdir(exist_ok=True)
    for name,rows in tables.items():
        with (out/(name+'.csv')).open('w',newline='',encoding='utf-8') as f:
            w=csv.writer(f);w.writerow(HEADERS);w.writerows(rows)
    with (out/'Comparison_summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['Sheet','Compartment','Contrast (numerator vs denominator)','Numerator groups','Denominator groups','Nominal n numerator','Nominal n denominator','Proteins tested','Volcano hits','|pi| > 1.0','|pi| > 1.3','BH q < 0.05','Validation']);w.writerows(summaries)
