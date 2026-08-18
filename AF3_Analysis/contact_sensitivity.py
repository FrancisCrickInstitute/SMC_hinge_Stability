#!/usr/bin/env python3
"""Contact-cutoff sensitivity: compute contacts/interface residues at 5 Å (atom-pair,
as in the main pipeline) and at 8 Å (residue-level, min heavy-atom distance) for all
1000 harmonised models. Parallel, resumable.
Usage: python3 contact_sensitivity.py --root <root> --out <csv> --budget 38 --workers 16
"""
import warnings; warnings.filterwarnings('ignore')
import os, sys, csv, time, argparse, multiprocessing as mp
from Bio.PDB import MMCIFParser, NeighborSearch
cif=MMCIFParser(QUIET=True)
SPECIES=['human','mouse','cerevisiae','pombe']

def metrics(model):
    aA=[a for a in model['A'].get_atoms()]; aB=[a for a in model['B'].get_atoms()]
    ns=NeighborSearch(aA+aB)
    c5_atom=0; rp={}; irA5=set(); irB5=set(); irA8=set(); irB8=set()
    for a,b in ns.search_all(8.0):
        ca=a.get_parent().get_parent().id; cb=b.get_parent().get_parent().id
        if ca==cb: continue
        d=a-b
        ra=a.get_parent(); rb=b.get_parent()
        # order as (A-residue, B-residue)
        if ca=='A': keyA=ra.id[1]; keyB=rb.id[1]
        else: keyA=rb.id[1]; keyB=ra.id[1]
        k=(keyA,keyB); rp[k]=min(rp.get(k,99), d)
        if d<=8.0: irA8.add(keyA); irB8.add(keyB)
        if d<=5.0:
            c5_atom+=1; irA5.add(keyA); irB5.add(keyB)
    rp5=sum(1 for v in rp.values() if v<=5.0)
    rp8=sum(1 for v in rp.values() if v<=8.0)
    return dict(c5_atom=c5_atom, respairs_5A=rp5, respairs_8A=rp8,
                iface_res_5A=len(irA5)+len(irB5), iface_res_8A=len(irA8)+len(irB8))

def build_jobs(root):
    jobs=[]
    for sp in SPECIES:
        cx=f'{sp}_cohesin_hinge_harmonized'; base=os.path.join(root,sp,'results')
        for N in range(1,51):
            sd=os.path.join(base,f'seed_{N}',cx)
            for k in range(5):
                tag=f'{cx}_seed-{N}_sample-{k}'
                jobs.append((sp,N,k,os.path.join(sd,f'seed-{N}_sample-{k}',f'{tag}_model.cif')))
    return jobs

def work(job):
    sp,N,k,model=job
    try:
        if not os.path.exists(model): return {'__skip__':1}
        row=dict(species=sp,seed=N,sample=k); row.update(metrics(cif.get_structure('m',model)[0])); return row
    except Exception as e: return {'__err__':f'{sp} {N}.{k}: {e}'}

FIELDS=['species','seed','sample','c5_atom','respairs_5A','respairs_8A','iface_res_5A','iface_res_8A']
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
