#!/usr/bin/env python3
"""
Generate the state-matched, filtered cross-species comparison + figures from a
master metrics CSV produced by hinge_pipeline.py.

Usage:
    python3 hinge_report.py --master master_metrics_rerun.csv --outdir <dir> \
            [--exp experimental_metrics.csv]
"""
import warnings; warnings.filterwarnings('ignore')
import argparse, os
import numpy as np, pandas as pd
from scipy.stats import kruskal
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
try:
    from sklearn.mixture import GaussianMixture; import diptest; HAVE=True
except Exception: HAVE=False

ORDER=['human','mouse','cerevisiae','pombe']
COL={'human':'#2c7fb8','mouse':'#41ab5d','cerevisiae':'#f16913','pombe':'#cb181d'}
SPREAD_THR=13.0; GLOBAL_IPSAE=0.66

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--master',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--exp',default=None)
    a=ap.parse_args(); os.makedirs(a.outdir,exist_ok=True)
    m=pd.read_csv(a.master)
    m=m[m.species.isin(ORDER)].copy()
    m['state']=np.where(m.iface_spread>=SPREAD_THR,'closed','open')
    m['edens']=m.prodigy_dG/m.bsa*1000.0
    # per-species reliable (upper ipSAE mode; unimodal -> keep all)
    rel=np.zeros(len(m),bool)
    for sp in ORDER:
        idx=m[m.species==sp].index; x=m.loc[idx,'ipsae'].values
        keep=np.ones(len(x),bool)
        if HAVE and len(x)>10:
            try:
                dp,p=diptest.diptest(x)
                if p<0.05:
                    g=GaussianMixture(2,n_init=6,random_state=0).fit(x.reshape(-1,1))
                    hi=np.argmax(g.means_.ravel()); keep=g.predict(x.reshape(-1,1))==hi
            except Exception: pass
        rel[[i for i,k in zip(idx,keep) if k]]=True
    m['reliable_ps']=rel; m['reliable_glob']=m.ipsae>=GLOBAL_IPSAE
    m.to_csv(os.path.join(a.outdir,'master_metrics.csv'),index=False)

    mets=['bsa','n_contacts','n_saltbridge','prodigy_dG','edens']
    sub=m[(m.state=='closed')&(m.reliable_ps)]
    tbl=sub.groupby('species')[mets].median().reindex(ORDER)
    tbl['n']=sub.groupby('species').size().reindex(ORDER)
    tbl.to_csv(os.path.join(a.outdir,'closed_reliable_medians.csv'))
    lines=['# Rerun analysis summary','',f'Total models: {len(m)}',
           'Per-species counts: '+', '.join(f'{s}={int((m.species==s).sum())}' for s in ORDER),'',
           '## Closed-state, per-species-reliable medians','',tbl.round(2).to_string(),'','## Kruskal-Wallis (closed, reliable)']
    for c in mets:
        groups=[sub[sub.species==s][c].values for s in ORDER if len(sub[sub.species==s])>2]
        if len(groups)>=2:
            H,p=kruskal(*groups); lines.append(f'  {c}: H={H:.1f} p={p:.2e}')
    open(os.path.join(a.outdir,'SUMMARY.md'),'w').write('\n'.join(lines)+'\n')

    # Fig1 landscape
    fig,ax=plt.subplots(figsize=(9,5))
    for i,sp in enumerate(ORDER):
        y=m[m.species==sp].iface_spread.values; ax.scatter(np.random.normal(i,0.08,len(y)),y,s=12,color=COL[sp],alpha=0.5,edgecolors='none')
        if len(y): ax.hlines(np.median(y),i-0.25,i+0.25,color='black',lw=2)
    ax.axhline(SPREAD_THR,ls='--',color='grey'); ax.set_xticks(range(4)); ax.set_xticklabels(ORDER)
    ax.set_ylabel('Interface spatial spread (Å)'); ax.set_title('Conformational landscape (rerun)',fontweight='bold')
    plt.tight_layout(); plt.savefig(os.path.join(a.outdir,'fig1_state_landscape.png'),dpi=140); plt.close()

    # Fig2 closed strength
    panels=[('bsa','Buried surface area (Å²)'),('n_saltbridge','Salt bridges'),('prodigy_dG','PRODIGY ΔG'),('edens','Energy density')]
    fig,axes=plt.subplots(1,4,figsize=(18,4.6))
    for ax,(c,lab) in zip(axes,panels):
        data=[sub[sub.species==sp][c].values for sp in ORDER]
        if all(len(d)>0 for d in data):
            bp=ax.boxplot(data,patch_artist=True,widths=0.6,showfliers=False)
            for pt,sp in zip(bp['boxes'],ORDER): pt.set_facecolor(COL[sp]); pt.set_alpha(0.75)
            for md in bp['medians']: md.set_color('black')
        ax.set_xticks(range(1,5)); ax.set_xticklabels([s[:4] for s in ORDER],fontsize=9); ax.set_title(lab,fontweight='bold',fontsize=10.5); ax.grid(axis='y',alpha=0.25)
    fig.suptitle('Closed-state interface strength, per-species reliable (rerun)',fontweight='bold',fontsize=13)
    plt.tight_layout(rect=[0,0,1,0.94]); plt.savefig(os.path.join(a.outdir,'fig2_closed_strength.png'),dpi=140); plt.close()

    # Fig3 filter comparison
    regimes=[('ALL',m),('per-species reliable',m[m.reliable_ps]),('global ipSAE',m[m.reliable_glob]),('CLOSED+reliable',m[(m.state=='closed')&m.reliable_ps])]
    fig,axes=plt.subplots(1,2,figsize=(15,5)); w=0.2
    for ax,(metric,lab) in zip(axes,[('bsa','Median BSA (Å²)'),('prodigy_dG','Median ΔG')]):
        for j,(rn,dd) in enumerate(regimes):
            vals=[dd[dd.species==s][metric].median() if len(dd[dd.species==s]) else np.nan for s in ORDER]
            ax.bar(np.arange(4)+j*w-1.5*w,vals,w,label=rn,alpha=0.9)
        ax.set_xticks(range(4)); ax.set_xticklabels(ORDER); ax.set_title(lab,fontweight='bold'); ax.grid(axis='y',alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle('Filtering & state-matching effect (rerun)',fontweight='bold'); plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(os.path.join(a.outdir,'fig3_filter_comparison.png'),dpi=140); plt.close()
    print('report written to',a.outdir)

if __name__=='__main__': main()
