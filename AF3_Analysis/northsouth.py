#!/usr/bin/env python3
"""For OPEN models, localise the retained sub-interface: mean Smc1 & Smc3 interface
residue index of the (dominant) contact patch. Bimodality across open models =
two open sub-states (north-open vs south-open). Resumable.
Usage: python3 northsouth.py --master <1000csv> --root <rerun_root> --out <csv> --budget 38
"""
import warnings; warnings.filterwarnings('ignore')
import os,sys,csv,time,argparse
import numpy as np, pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch
from sklearn.cluster import KMeans
cif=MMCIFParser(QUIET=True)

def retained_patch(model):
    aA=[a for a in model['A'].get_atoms()]; aB=[a for a in model['B'].get_atoms()]
    ns=NeighborSearch(aA+aB); mids=[]; info=[]
    for a,b in ns.search_all(5.0):
        ra=a.get_parent(); rb=b.get_parent()
        if ra.get_parent().id==rb.get_parent().id: continue
        mids.append((a.coord+b.coord)/2)
        # smc1 (A) residue index, smc3 (B) residue index
        ai = ra.id[1] if a.get_parent().get_parent().id=='A' else rb.id[1]
        bi = rb.id[1] if b.get_parent().get_parent().id=='B' else ra.id[1]
        info.append((ai,bi))
    if len(mids)<8: return None
    mids=np.array(mids); info=np.array(info)
    km=KMeans(2,n_init=4,random_state=0).fit(mids); lab=km.labels_
    # dominant patch = larger cluster
    dom=0 if (lab==0).sum()>=(lab==1).sum() else 1
    sel=info[lab==dom]
    return float(sel[:,0].mean()), float(sel[:,1].mean()), int((lab==dom).sum())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--master',required=True); ap.add_argument('--root',required=True)
    ap.add_argument('--out',required=True); ap.add_argument('--budget',type=float,default=38.0)
    a=ap.parse_args()
    m=pd.read_csv(a.master); m=m[m.iface_spread<13.0]  # open models
    done=set()
    if os.path.exists(a.out):
        for r in pd.read_csv(a.out).itertuples(): done.add((r.species,int(r.seed),int(r.sample)))
    new=not os.path.exists(a.out); fh=open(a.out,'a',newline=''); w=csv.writer(fh)
    if new: w.writerow(['species','seed','sample','smc1_idx','smc3_idx','patch_n'])
    t0=time.time(); n=0
    for r in m.itertuples():
        if (r.species,r.seed,r.sample) in done: continue
        if time.time()-t0>a.budget: break
        complex=f'{r.species}_cohesin_hinge_harmonized'
        tag=f'{complex}_seed-{r.seed}_sample-{r.sample}'
        p=os.path.join(a.root,r.species,'results',f'seed_{r.seed}',complex,f'seed-{r.seed}_sample-{r.sample}',f'{tag}_model.cif')
        try:
            res=retained_patch(cif.get_structure('m',p)[0])
            if res: w.writerow([r.species,r.seed,r.sample,round(res[0],1),round(res[1],1),res[2]]); fh.flush(); n+=1
        except Exception as e: sys.stderr.write(f'ERR {r.species}{r.seed}.{r.sample}: {e}\n')
    fh.close(); print(f'processed={n} total={len(pd.read_csv(a.out))}')

if __name__=='__main__': main()
