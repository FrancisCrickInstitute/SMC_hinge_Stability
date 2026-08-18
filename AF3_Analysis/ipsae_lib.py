import json, numpy as np
def compute_ipsae(conf_path, pae_cutoff=15.0):
    """Interface Score of Aligned Errors (Dunbrack), from AF3 confidences.json PAE matrix."""
    d=json.load(open(conf_path))
    pae=np.array(d['pae'],dtype=np.float32)
    ch=np.array(d['token_chain_ids'])
    chains=list(dict.fromkeys(ch.tolist()))
    best=0.0
    for c1 in chains:
        for c2 in chains:
            if c1==c2: continue
            i_idx=np.where(ch==c1)[0]; j_idx=np.where(ch==c2)[0]
            sub=pae[np.ix_(i_idx,j_idx)]
            best_i=0.0
            for r in range(sub.shape[0]):
                valid=sub[r]<pae_cutoff
                n0=int(valid.sum())
                if n0<1: continue
                d0=max(1.24*np.cbrt(max(n0,19)-15)-1.8,1.0)
                ptm=np.mean(1.0/(1.0+(sub[r][valid]/d0)**2))
                if ptm>best_i: best_i=ptm
            if best_i>best: best=best_i
    return float(best)
