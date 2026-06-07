import random, math, re
from collections import defaultdict
from Bio.Data import CodonTable
try: from seqfold import fold
except: fold = None
try:
    from numba import cuda, jit
    import numpy as np
    CUDA_AVAILABLE = cuda.is_available()
except: CUDA_AVAILABLE = False

def get_codon_table(table_id=1):
    try:
        t = CodonTable.unambiguous_dna_by_id[table_id].forward_table
        stops = CodonTable.unambiguous_dna_by_id[table_id].stop_codons
    except:
        t = CodonTable.unambiguous_dna_by_id[1].forward_table
        stops = CodonTable.unambiguous_dna_by_id[1].stop_codons
    table = dict(t)
    for stop in stops: table[stop] = '*'
    return table

def reverse_codon_table(table_id=1):
    rev = defaultdict(list)
    for c, aa in get_codon_table(table_id).items(): rev[aa].append(c)
    return rev

def calculate_cai(dna, freq, table_id=1):
    if not dna: return 0.0
    dfreq = {k.replace('U','T'): v for k,v in freq.items()}
    w, rev = {}, reverse_codon_table(table_id)
    for aa, codons in rev.items():
        m = max([dfreq.get(c, 0) for c in codons]) if codons else 0
        for c in codons: w[c] = dfreq.get(c, 0) / m if m > 0 else 0.01
    s, n = 0, 0
    for i in range(0, len(dna), 3):
        c = dna[i:i+3]
        if len(c) == 3: s += math.log(max(w.get(c, 0.01), 0.01)); n += 1
    return math.exp(s / n) if n > 0 else 0.0

def calculate_mfe(dna):
    if not fold or not dna: return 0.0
    if len(dna) > 500:
        if len(dna) > 5000 and CUDA_AVAILABLE:
            try: return calculate_mfe_gpu(dna)
            except: pass
        s, n = 0, 0
        for i in range(0, len(dna) - 50, 50):
            s += sum(st.e for st in fold(dna[i:i+100].replace('T','U'))); n += 1
        return round(s, 2) if n > 0 else 0.0
    return round(sum(st.e for st in fold(dna.replace('T','U'))), 2)

@cuda.jit
def _nussinov_kernel(seqs, res, W):
    idx = cuda.grid(1)
    if idx < seqs.shape[0]:
        s = seqs[idx]
        for l in range(1, W):
            for i in range(W - l):
                j = i + l
                r = max(res[idx, (i+1)*W + j], res[idx, i*W + (j-1)])
                p = 1 if ((s[i]==1 and s[j]==4) or (s[i]==4 and s[j]==1) or (s[i]==2 and s[j]==3) or (s[i]==3 and s[j]==2)) else 0
                r = max(r, res[idx, (i+1)*W + (j-1)] + p)
                for k in range(i + 1, j): r = max(r, res[idx, i*W + k] + res[idx, (k+1)*W + j])
                res[idx, i*W + j] = r

def calculate_mfe_gpu(dna, W=100, S=50):
    m = {'A':1, 'C':2, 'G':3, 'T':4, 'U':4}
    wins = [ [m.get(b,0) for b in dna[i:i+W].upper()] for i in range(0, len(dna)-W+1, S) ]
    if not wins: return 0.0
    nw = len(wins)
    s_gpu = cuda.to_device(np.array(wins, dtype=np.int32))
    r_gpu = cuda.to_device(np.zeros((nw, W*W), dtype=np.int32))
    tp = 32
    _nussinov_kernel[(nw+tp-1)//tp, tp](s_gpu, r_gpu, W)
    res = r_gpu.copy_to_host()
    return round(sum(res[i, W-1] for i in range(nw)) * -2.0, 2)

def scan_cryptic_signals(dna):
    sigs = {"PolyA (AATAAA)": "AATAAA", "PolyA (ATTAAA)": "ATTAAA", "Cryptic Splice Donor": "CAGGTAAGT", "Cryptic Splice Acceptor": "TTCAG[GA]"}
    return {n: [m.start() for m in re.finditer(p, dna)] for n,p in sigs.items()}

def scan_restriction_sites(dna, slist):
    db = {'EcoRI':'GAATTC', 'BamHI':'GGATCC', 'HindIII':'AAGCTT', 'NcoI':'CCATGG', 'XhoI':'CTCGAG', 'NdeI':'CATATG', 'SalI':'GTCGAC', 'XbaI':'TCTAGA', 'SpeI':'ACTAGT', 'PstI':'CTGCAG'}
    res = {s: [] for s in slist}
    for s in slist:
        p, i = db.get(s, s), 0
        while True:
            i = dna.find(p, i)
            if i == -1: break
            res[s].append(i); i += 1
    return res

def optimize_frequency(aa, freq, gc_min, gc_max, sites, table_id=1):
    dfreq, rev = {k.replace('U','T'):v for k,v in freq.items()}, reverse_codon_table(table_id)
    return ''.join([sorted(rev.get(a,['']), key=lambda c: dfreq.get(c,0), reverse=True)[0] for a in aa])

def optimize_harmonization(aa, freq, acc, table_id=1):
    dfreq, rev = {k.replace('U','T'):v for k,v in freq.items()}, reverse_codon_table(table_id)
    res = ''
    for a in aa:
        cands = rev.get(a, [])
        if not cands: continue
        w = [dfreq.get(c,0) for c in cands]
        res += random.choices(cands, weights=w, k=1)[0] if sum(w)>0 else cands[0]
    return res

def optimize_probabilistic(aa, freq, n, gc_min, gc_max, sites=[], table_id=1, avoid_signals=False):
    dfreq, rev = {k.replace('U','T'):v for k,v in freq.items()}, reverse_codon_table(table_id)
    res = []
    for _ in range(n * 2):
        s = ''
        for a in aa:
            cands = rev.get(a, [])
            if not cands: continue
            w = [dfreq.get(c,0) for c in cands]
            s += random.choices(cands, weights=w, k=1)[0] if sum(w)>0 else random.choice(cands)
        gc = (s.count('G')+s.count('C'))/len(s)*100 if s else 0
        if gc < gc_min or gc > gc_max: continue
        if sites and any(len(p)>0 for p in scan_restriction_sites(s, sites).values()): continue
        if avoid_signals and any(len(p)>0 for p in scan_cryptic_signals(s).values()): continue
        res.append((s, calculate_cai(s, freq, table_id)))
        if len(res) >= n: break
    if not res:
        for _ in range(n):
            s = ''.join([random.choices(rev[a], weights=[dfreq.get(c,0) for c in rev[a]])[0] for a in aa if a in rev])
            res.append((s, calculate_cai(s, freq, table_id)))
    return sorted(res, key=lambda x: x[1], reverse=True)

def sliding_window_cai(dna, freq, w=10, table_id=1):
    res = [calculate_cai(dna[i:i+w*3], freq, table_id) for i in range(0, len(dna)-w*3+1, 3)]
    while len(res) < len(dna)//3: res.append(res[-1] if res else 0)
    return res

def sliding_window_gc(dna, w=30):
    res = [(dna[i:i+w].count('G')+dna[i:i+w].count('C'))/w*100 for i in range(len(dna)-w+1)]
    while len(res) < len(dna): res.append(res[-1] if res else 0)
    return res
