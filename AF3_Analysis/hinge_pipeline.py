#!/usr/bin/env python3
"""
Resumable, self-contained per-model interface pipeline for the cohesin hinge study.
Auto-discovers AF3 models under a root dir (works for both alphastream and AF3.0.2
split-per-seed layouts). Computes: interface ipTM, interface min-PAE, ipTM, pTM,
ranking_score, ipSAE, BSA, interface residues, contacts, H-bond-like, salt bridges,
PRODIGY dG/Kd, and the conformational-state descriptors (spread, elongation).

Usage (call repeatedly until it prints 100%):
    python3 hinge_pipeline.py --root <results_root> --out master_metrics_rerun.csv --budget 40

Species is inferred from the file path (human/mouse/cerevisiae|cereviae/pombe).
"""
import warnings; warnings.filterwarnings('ignore')
import os, re, sys, csv, json, time, argparse, glob
import numpy as np
import freesasa; freesasa.setVerbosity(freesasa.silent)
from Bio.PDB import MMCIFParser, PDBIO, Select, NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from prodigy_prot.modules.prodigy import Prodigy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipsae_lib import compute_ipsae

NEG={('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2')}
POS={('LYS','NZ'),('ARG','NH1'),('ARG','NH2'),('ARG','NE'),('HIS','ND1'),('HIS','NE2')}
cif=MMCIFParser(QUIET=True); pdbp=PDBParser(QUIET=True); io=PDBIO()

def species_of(path):
    p=path.lower()
    if 'human' in p: return 'human'
    if 'mouse' in p: return 'mouse'
    if 'cerevisiae' in p or 'cereviae' in p: return 'cerevisiae'
    if 'pombe' in p: return 'pombe'
    return 'unknown'

def discover(root):
    """Yield (species, seed, sample, model_cif, conf_json) tuples."""
    summ=[]
    for dp,_,fns in os.walk(root):
        for fn in fns:
            if fn.endswith('.json') and 'summary_confidences' in fn:
                summ.append(os.path.join(dp,fn))
    for s in summ:
        model=s.replace('summary_confidences','model').replace('.json','.cif')
        conf=s.replace('summary_confidences','confidences')
        if not os.path.exists(model): continue
        m=re.search(r'seed[-_]?(\d+)', s); seed=int(m.group(1)) if m else -1
        m2=re.search(r'sample[-_]?(\d+)', s); sample=int(m2.group(1)) if m2 else 0
        yield species_of(s), seed, sample, model, conf, s

class Sel(Select):
    def __init__(s,c): s.c=c
    def accept_chain(s,ch): return ch.id==s.c
def sasa(f): return freesasa.calc(freesasa.Structure(f)).totalArea()

def metrics_for(model, conf, summ):
    d=json.load(open(summ))
    cpi=d.get('chain_pair_iptm'); cpp=d.get('chain_pair_pae_min')
    out=dict(iptm=d.get('iptm'), ptm=d.get('ptm'), ranking_score=d.get('ranking_score'),
             interface_iptm=(cpi[0][1] if cpi else None),
             interface_pae_min=(min(cpp[0][1],cpp[1][0]) if cpp else None))
    T=f'/tmp/_{os.getpid()}_'
    s=cif.get_structure('m',model); model0=s[0]
    chains=[c.id for c in model0][:2]; cA,cB=chains[0],chains[1]
    io.set_structure(s)
    io.save(T+'cx.pdb'); io.save(T+'A.pdb',Sel(cA)); io.save(T+'B.pdb',Sel(cB))
    out['bsa']=round(sasa(T+'A.pdb')+sasa(T+'B.pdb')-sasa(T+'cx.pdb'),1)
    aA=list(model0[cA].get_atoms()); aB=list(model0[cB].get_atoms())
    ns=NeighborSearch(aA+aB); resA=set(); resB=set(); ncont=nhb=nsb=0; mids=[]
    for a,b in ns.search_all(5.0):
        ca=a.get_parent().get_parent().id; cb=b.get_parent().get_parent().id
        if ca==cb: continue
        ncont+=1; mids.append((a.coord+b.coord)/2); dd=a-b
        ra=a.get_parent(); rb=b.get_parent()
        (resA if ca==cA else resB).add(ra.id[1]); (resA if cb==cA else resB).add(rb.id[1])
        if dd<3.5 and a.name[0] in 'NO' and b.name[0] in 'NO': nhb+=1
        ka=(ra.resname,a.name); kb=(rb.resname,b.name)
        if dd<4.0 and ((ka in NEG and kb in POS) or (ka in POS and kb in NEG)): nsb+=1
    out.update(iface_res_A=len(resA), iface_res_B=len(resB), n_contacts=ncont,
               n_hbond_like=nhb, n_saltbridge=nsb)
    mids=np.array(mids)
    if len(mids)>=10:
        c=mids-mids.mean(0); ev=np.sort(np.linalg.eigvalsh(np.cov(c.T)))[::-1]
        out['iface_spread']=round(float(np.sqrt(ev.sum())),2)
        out['iface_elong']=round(float(np.sqrt(ev[0]/max(ev[1],1e-6))),3)
    else:
        out['iface_spread']=0.0; out['iface_elong']=0.0
    st=pdbp.get_structure('c',T+'cx.pdb')
    prod=Prodigy(st[0],selection=[cA,cB],temp=25.0); prod.predict()
    out['prodigy_dG']=round(prod.ba_val,2); out['prodigy_Kd']=prod.kd_val
    out['ipsae']=round(compute_ipsae(conf),4) if os.path.exists(conf) else None
    return out

FIELDS=['species','seed','sample','iptm','ptm','ranking_score','interface_iptm','interface_pae_min',
        'ipsae','bsa','iface_res_A','iface_res_B','n_contacts','n_hbond_like','n_saltbridge',
        'iface_spread','iface_elong','prodigy_dG','prodigy_Kd']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True); ap.add_argument('--out',required=True)
    ap.add_argument('--budget',type=float,default=40.0)
    a=ap.parse_args()
    items=list(discover(a.root))
    done=set()
    if os.path.exists(a.out):
        import pandas as pd
        for r in pd.read_csv(a.out).itertuples(): done.add((r.species,r.seed,r.sample))
    new=not os.path.exists(a.out); fh=open(a.out,'a',newline=''); w=csv.DictWriter(fh,fieldnames=FIELDS)
    if new: w.writeheader()
    t0=time.time(); n=0
    for sp,seed,sample,model,conf,summ in items:
        if (sp,seed,sample) in done: continue
        if time.time()-t0>a.budget: break
        try:
            row=dict(species=sp,seed=seed,sample=sample); row.update(metrics_for(model,conf,summ))
            w.writerow(row); fh.flush(); n+=1
        except Exception as e:
            sys.stderr.write(f'ERR {sp} {seed}.{sample}: {e}\n')
    fh.close()
    import pandas as pd
    tot=len(pd.read_csv(a.out)) if os.path.exists(a.out) else 0
    print(f'discovered={len(items)} processed_this_run={n} total_done={tot}')

if __name__=='__main__': main()
