#!/usr/bin/env python3
"""Parallel deterministic-path driver (all diffusion samples). ~4x faster.
Usage: python3 fast_rerun_par.py --root <root> --out <csv> --budget 38 [--workers 4]
"""
import warnings; warnings.filterwarnings('ignore')
import os, sys, csv, time, argparse, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hinge_pipeline import metrics_for, FIELDS

SPECIES=['human','mouse','cerevisiae','pombe']
def build_jobs(root, maxseed=50, nsample=5):
    """Deterministic paths, NO stat checks (avoids slow network stats)."""
    jobs=[]
    for sp in SPECIES:
        complex=f'{sp}_cohesin_hinge_harmonized'; base=os.path.join(root,sp,'results')
        for N in range(1,maxseed+1):
            sd=os.path.join(base,f'seed_{N}',complex)
            for k in range(nsample):
                tag=f'{complex}_seed-{N}_sample-{k}'; d=os.path.join(sd,f'seed-{N}_sample-{k}')
                jobs.append((sp,N,k,os.path.join(d,f'{tag}_model.cif'),
                             os.path.join(d,f'{tag}_confidences.json'),
                             os.path.join(d,f'{tag}_summary_confidences.json')))
    return jobs

def work(job):
    sp,N,k,model,conf,summ=job
    try:
        if not os.path.exists(summ): return {'__skip__':1}
        row=dict(species=sp,seed=N,sample=k); row.update(metrics_for(model,conf,summ)); return row
    except Exception as e:
        return {'__err__':f'{sp} {N}.{k}: {e}'}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True); ap.add_argument('--out',required=True)
    ap.add_argument('--budget',type=float,default=38.0); ap.add_argument('--workers',type=int,default=4)
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
            if time.time()-t0>a.budget:
                pool.terminate(); break
    fh.close()
    import pandas as pd
    tot=len(pd.read_csv(a.out)); left=len(jobs)-n
    print(f'processed={n} err={err} total={tot} remaining~={max(left,0)}')

if __name__=='__main__': main()
