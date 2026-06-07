from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import os, glob
from snapgene_reader import snapgene_file_to_dict

VECTORS = {}

def load_all_vectors(directory="data/vectors"):
    global VECTORS
    VECTORS = {"Common Vectors": {"pUC19 (E. coli)": {"sequence": "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGA", "insertion_site": 396, "is_circular": True, "features": [("AmpR", 150, 250, +1, "CDS"), ("ori", 600, 700, +1, "rep_origin")]}}}
    if not os.path.exists(directory): return
    for fpath in glob.glob(os.path.join(directory, "**/*.dna"), recursive=True):
        try:
            name, cat = os.path.basename(fpath).replace(".dna", ""), os.path.basename(os.path.dirname(fpath))
            if cat not in VECTORS: VECTORS[cat] = {}
            d = snapgene_file_to_dict(fpath)
            s = d['seq']
            site = len(s)//2
            for p in ["GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "TCTAGA"]:
                idx = s.find(p)
                if idx != -1: site = idx; break
            feats = [(f.get('name', f.get('type', 'feat')), f['start'], f['end'], 1 if f.get('strand')=='+' else (-1 if f.get('strand')=='-' else 0), f.get('type', 'misc')) for f in d['features']]
            VECTORS[cat][name] = {"sequence": s, "insertion_site": site, "is_circular": True, "features": feats}
        except: pass

load_all_vectors()

def assemble_vector(cat, name, dna, id_name="Optimized Insert"):
    if cat not in VECTORS or name not in VECTORS[cat]:
        v = next((v for c in VECTORS.values() for n, v in c.items() if name in n), None)
        if not v: return None, None
    else: v = VECTORS[cat][name]
    site = v["insertion_site"]
    fseq = v["sequence"][:site] + dna + v["sequence"][site:]
    feats = []
    for fn, fs, fe, fst, ft in v["features"]:
        if fs >= site: fs += len(dna); fe += len(dna)
        elif fs < site and fe > site: fe += len(dna)
        c = "#ccebc5"
        if "promoter" in ft.lower() or "ITR" in fn: c = "#fbb4ae"
        elif "polyA" in ft.lower() or "terminator" in ft.lower(): c = "#decbe4"
        elif "resistance" in fn.lower() or "CDS" in ft: c = "#b3cde3"
        elif "primer" in ft.lower(): c = "#fed9a6"
        elif "restriction" in ft.lower() or "enzyme" in fn.lower(): c = "#ffffcc"
        elif "ORF" in ft or "gene" in ft.lower(): c = "#b3e2cd"
        feats.append(GraphicFeature(start=fs, end=fe, strand=fst, color=c, label=fn))
    feats.append(GraphicFeature(start=site, end=site+len(dna), strand=+1, color="#00f5ff", label=id_name))
    rec = CircularGraphicRecord(sequence_length=len(fseq), features=feats) if v.get("is_circular") else GraphicRecord(sequence_length=len(fseq), sequence=fseq, features=feats)
    rec.sequence = fseq
    return fseq, rec

def plot_labels_enhanced(rec, ax, lvl):
    lbf = [f for f in rec.features if f.label]
    if not lbf: return
    max_l = max(lvl.values()) if lvl else 0
    br, lbls = rec.radius + (max_l + 2.5) * rec.feature_level_height, []
    for f in lbf:
        a = rec.position_to_angle(f.x_center)
        ar = np.deg2rad(a)
        fr = rec.radius + (lvl.get(f, 0) + 0.5) * rec.feature_level_height
        lbls.append({'f': f, 'fx': fr*np.cos(ar), 'fy': fr*np.sin(ar)-rec.radius, 'lx': br*np.cos(ar), 'ly': br*np.sin(ar)-rec.radius, 'a': a, 't': f.label})
    L, R = sorted([l for l in lbls if l['lx'] < 0], key=lambda x: x['ly']), sorted([l for l in lbls if l['lx'] >= 0], key=lambda x: x['ly'])
    def spread(g):
        for _ in range(50):
            c = False
            for i in range(len(g)-1):
                if (g[i+1]['ly'] - g[i]['ly']) < 0.25:
                    m = (g[i]['ly'] + g[i+1]['ly'])/2
                    g[i]['ly'], g[i+1]['ly'], c = m-0.125, m+0.125, True
            if not c: break
    spread(L); spread(R)
    for l in L+R:
        ha = 'right' if l['lx'] < 0 else 'left'
        na = l['a'] % 360
        if (80 < na < 100) or (260 < na < 280): ha = 'center'
        ax.plot([l['fx'], l['lx']], [l['fy'], l['ly']], color='#ffffff', lw=0.6, alpha=0.5, zorder=1)
        ax.text(l['lx'], l['ly'], l['t'], ha=ha, va='center', fontsize=9, fontweight='bold', color='#ffffff', path_effects=[path_effects.withStroke(linewidth=2, foreground="black")], zorder=3)

def render_vector_map(rec, out="vector_map.png"):
    circ = isinstance(rec, CircularGraphicRecord)
    plt.rcParams.update({'text.color': '#ffffff', 'axes.labelcolor': '#ffffff'})
    fig, ax = plt.subplots(1, figsize=(14, 14 if circ else 4))
    fig.patch.set_facecolor('#0a0d14'); ax.set_facecolor('#0a0d14')
    if circ:
        orig = [f.label for f in rec.features]
        for f in rec.features: f.label = None
        ax, (lvl, _) = rec.plot(ax=ax, with_ruler=False)
        for f, l in zip(rec.features, orig): f.label = l
        plot_labels_enhanced(rec, ax, lvl)
    else:
        rec.plot(ax=ax, with_ruler=True, ruler_color='#ffffff')
        for t in ax.texts:
            t.set_color('#ffffff'); t.set_weight('bold'); t.set_fontsize(10)
            t.set_path_effects([path_effects.withStroke(linewidth=2, foreground="black")])
    if not circ:
        ax.spines['bottom'].set_color('#e2e8f0'); ax.tick_params(axis='x', colors='#e2e8f0')
    else: ax.axis('off')
    plt.tight_layout()
    fig.savefig(out, dpi=300, facecolor='#0a0d14', transparent=False); plt.close(fig)
    return out

def render_sequence_zoom(rec, site, win=20, out="sequence_zoom.png"):
    plt.rcParams.update({'text.color': '#ffffff', 'axes.labelcolor': '#ffffff'})
    s, e = max(0, site-win), min(rec.sequence_length, site+win+20)
    cf = [GraphicFeature(start=max(0, f.start-s), end=min(e-s, f.end-s), strand=f.strand, color=f.color, label=f.label) for f in rec.features if f.start < e and f.end > s]
    crec = GraphicRecord(sequence_length=e-s, sequence=rec.sequence[s:e], features=cf)
    fig, ax = plt.subplots(1, figsize=(14, 4))
    fig.patch.set_facecolor('#0a0d14'); ax.set_facecolor('#0a0d14')
    crec.plot(ax=ax, with_ruler=True, ruler_color='#ffffff'); crec.plot_sequence(ax=ax)
    for t in ax.texts:
        t.set_color('#ffffff'); t.set_weight('bold')
        t.set_path_effects([path_effects.withStroke(linewidth=2, foreground="black")])
    ax.spines['bottom'].set_color('#e2e8f0'); ax.tick_params(axis='x', colors='#e2e8f0')
    plt.tight_layout(); fig.savefig(out, dpi=300, facecolor='#0a0d14', transparent=False); plt.close(fig)
    return out
