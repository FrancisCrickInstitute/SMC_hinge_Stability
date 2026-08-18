#!/usr/bin/env python3
"""Canonical Cβ–Cβ contact map (Cα for glycine): inter-chain residue pairs with
Cβ–Cβ distance <= 8 Å (and <= 5 Å). Parallel, resumable.
Usage: python3 cb_contacts.py --root <root> --out <csv> --budget 38 --workers 16
"""
import warnings; warnings.filterwarnings('ignore')
import os, sys, csv, time, argparse, multiprocessing as mp
import numpy as np
from Bio.PDB import MMCIFParser
from scipy.spatial.distance import cdist
cif=MMCIFParser(QUIET=True)
SPECIES=['human','mouse','cerevisiae','pombe']

def cbeta(chain):
    coords=[]
    for r in chain:
        if r.id[0]!=' ': continue
        atom='CB' if ('CB' in r) else ('CA' if 'CA' in r else None)
        if atom: coords.append(r[atom].coord)
    return np.array(coords)

def metrics(model):
    A=cbeta(model['A']); B=cbeta(model['B'])
    D=cdist(A,B)
    return dict(cb_contacts_8A=int((D<=8.0).sum()),
                cb_contacts_5A=int((D<=5.0).sum()),
                cb_iface_res_8A=int((D<=8.0).any(1).sum()+(D<=8.0).any(0).sum()))

def build_jobs(root):
    jobs=[]
    for sp in SPECIES:
        cx=f'{sp}_cohesin_hinge_harmonized'; base=os.path.join(root,sp,'results')
        for N in range(1,51):
            for k in range(5):
                tag=f'{cx}_seed-{N}_sample-{k}'
                jobs.append((sp,N,k,os.path.join(base,f'seed_{N}',cx,f'seed-{N}_sample-{k}',f'{tag}_model.cif')))
    return jobs

def work(job):
    sp,N,k,model=job
    try:
        if not os.path.exists(model): return {'__skip__':1}
        row=dict(species=sp,seed=N,sample=k); row.update(metrics(cif.get_structure('m',model)[0])); return row
    except Exception as e: return {'__err__':f'{sp} {N}.{k}: {e}'}

FIELDS=['species','seed','sample','cb_contacts_5A','cb_contacts_8A','cb_iface_res_8A']
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True); ap.add_argument('--out',required=True)
    ap.add_argument('--budget',type=float,default=38.0); ap.add_argument('--workers',type=int,default=16)
    a=ap.parse_args()
    done=set()
    if os.path.exists(a.out):
        import pandas as pd
        for r in pd.read_csv(a.out).itertuples(): done.add((r.species,int(r.seed),int(r.sample)))
    jobs=[j for j in build_jobs(a.root) if (j[0],j[1],j[2]) not in done]
    new=not os.path.exists(a.out); fh=open(a.out,'a',newline=''); w=csv.DictWriter(fh,fieldnames=FIELDS)
    if new: w.writeheader()
    t0=time.time(); n=0; err=0
    with mp.Pool(a.workers) as pool:
        for row in pool.imap_unordered(work, jobs, chunksize=1):
            if '__skip__' in row: pass
            elif '__err__' in row: err+=1; sys.stderr.write('ERR '+row['__err__']+'\n')
            else: w.writerow(row); fh.flush(); n+=1
            if time.time()-t0>a.budget: pool.terminate(); break
    fh.close(); import pandas as pd
    print(f'processed={n} err={err} total={len(pd.read_csv(a.out))} remaining~={len(jobs)-n}')
if __name__=='__main__': main()
