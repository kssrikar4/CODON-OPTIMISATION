import streamlit as st
import streamlit.components.v1 as components
import time, os, base64, json, glob
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import codon_loader, optimizer, report_generator, vector_integration

st.set_page_config(layout='wide', page_title='Codon Optimisation')

@st.cache_data
def load_db(path):
    raw = codon_loader.load_cocoputs(path)
    freq = codon_loader.normalize_to_frequency(raw)
    return freq, codon_loader.get_organism_list(freq)

st.markdown("""<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600&family=Rajdhani:wght@500;700&family=Share+Tech+Mono&display=swap');
    .stApp, .main, [data-testid="stAppViewContainer"] { background-color: #0a0d14 !important; color: #e2e8f0; font-family: 'Inter', sans-serif; }
    [data-testid="stSidebar"] { min-width: 350px; max-width: 350px; }
    h1, h2, h3, h4, h5, h6 { font-family: 'Rajdhani', sans-serif !important; color: #00f5ff !important; text-transform: uppercase; }
    .section-header { font-family: 'Rajdhani', sans-serif; color: #39ff14; font-weight: 700; border-bottom: 1px solid rgba(57,255,20,0.2); padding-bottom: 5px; margin: 20px 0 10px 0; }
    .terminal-box { background: #05070a; border: 1px solid #1e293b; padding: 15px; font-family: 'Share Tech Mono', monospace; color: #39ff14; overflow-y: auto; max-height: 300px; }
</style>""", unsafe_allow_html=True)

for k in ['optimization_complete', 'optimized_sequence', 'optimization_metrics', 'optimization_logs', 'found_sites', 'cryptic_signals']:
    if k not in st.session_state: st.session_state[k] = False if 'complete' in k else ('' if 'sequence' in k else {})

with st.sidebar:
    st.markdown('<div style="text-align: center; padding-bottom: 15px;"><h1 style="margin:0; font-size:1.8em; line-height:1.2;"><span style="color:#00f5ff;">CODON</span><span style="color:#39ff14;">OPTIMISATION</span></h1><div style="color: #64748b; font-size: 0.8em; font-family:\'Share Tech Mono\', monospace;">Sequence Engineering Platform</div><hr style="border: 0; height: 1px; background-image: linear-gradient(to right, rgba(0,0,0,0), rgba(0,245,255,0.75), rgba(0,0,0,0));"></div>', unsafe_allow_html=True)
    st.markdown('<div class="section-header">DATABASE CONFIGURATION</div>', unsafe_allow_html=True)
    st.markdown('<div style="background:rgba(0,245,255,0.05); padding:10px; border-left:3px solid #39ff14; font-size:0.85em; margin-bottom:15px;">Powered by CoCoPUTs / HIVE-CUTs. Please download data from FDA DNA HIVE and place it in the data/cocoputs/ directory.</div>', unsafe_allow_html=True)
    db_files = glob.glob('data/cocoputs/*.tsv')
    if db_files:
        sel_db = st.selectbox('Select Database File:', [os.path.basename(f) for f in db_files])
        freq_db, org_list = load_db(os.path.join('data/cocoputs/', sel_db))
    else: freq_db, org_list = {}, []
    if not org_list: st.error("Add TSV to data/cocoputs/"); st.stop()
    
    header = "TARGET ORGANISM" if "Refseq" in sel_db else "TARGET TISSUE"
    st.markdown(f'<div class="section-header">{header}</div>', unsafe_allow_html=True)
    
    def_idx = 0
    try:
        matches = [i for i, x in enumerate(org_list) if "Homo sapiens" in x]
        if matches: def_idx = matches[0]
    except: pass
    
    sel_org = st.selectbox('Search Taxonomy or Tissue:', org_list, index=def_idx, label_visibility='collapsed')
    run_btn = st.button('RUN OPTIMIZATION', type='primary', use_container_width=True)

tab1, tab2, tab3, tab4, tab5 = st.tabs(['SEQUENCE INPUT', 'OPTIMIZATION LOG', 'ANALYSIS', 'VECTOR INTEGRATION', 'EXPORT ASSETS'])

with tab1:
    col1, col2 = st.columns([6, 4])
    with col1:
        st.markdown('<div class="section-header">PROTEIN SEQUENCE INPUT</div>', unsafe_allow_html=True)
        raw_seq = st.text_area('Input Amino Acid Sequence:', value='MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN', height=200, label_visibility='collapsed')
        clean_seq = ''.join([c.upper() for c in raw_seq if c.isalpha()])
        if raw_seq.startswith('>'): clean_seq = ''.join([c.upper() for c in ''.join(raw_seq.split('\n')[1:]) if c.isalpha()])
        if all(c in 'ATGCUN' for c in clean_seq) and clean_seq: st.warning("High DNA character density detected. Please enter an Amino Acid sequence.")
        valid = set('ACDEFGHIKLMNPQRSTVWY')
        inv = set(clean_seq) - valid
        if inv: st.error(f"Invalid characters: {inv}")
        else: 
            st.markdown(f"<div style='color:#39ff14; font-family:\"Share Tech Mono\";'>Confirmed: {len(clean_seq)} AA residues loaded.</div>", unsafe_allow_html=True)
            comp = {'Hydrophobic': sum((clean_seq.count(a) for a in 'VILMFWC')), 'Polar': sum((clean_seq.count(a) for a in 'STYNQ')), 'Charged+': sum((clean_seq.count(a) for a in 'RHK')), 'Charged-': sum((clean_seq.count(a) for a in 'DE')), 'Special': sum((clean_seq.count(a) for a in 'GP'))}
            colors, svg, total = ['#00f5ff', '#39ff14', '#ff3b5c', '#ffb300', '#9d4edd'], '', len(clean_seq)
            if total > 0:
                x = 0
                for i, (k, v) in enumerate(comp.items()):
                    w = (v / total) * 100
                    if w > 0:
                        svg += f'<rect x="{x}%" y="0" width="{w}%" height="10" fill="{colors[i]}"><title>{k}: {v}</title></rect>'
                        x += w
                st.markdown(f'<svg width="100%" height="10" style="border-radius:5px; margin-top:10px;">{svg}</svg>', unsafe_allow_html=True)
    with col2:
        st.markdown('<div class="section-header">OPTIMIZATION PARAMETERS</div>', unsafe_allow_html=True)
        strat = st.selectbox('Strategy:', ['Codon Frequency Replacement', 'Codon Harmonization', 'Probabilistic Sampling'])
        st.caption(f"*Strategy: {('Replaces every AA with the most abundant codon.' if 'Frequency' in strat else 'Matches usage profile of a reference gene.' if 'Harmonization' in strat else 'Samples codons based on frequencies while avoiding constraints.')}*")
        gc_min, gc_max = st.slider('Target GC% Range:', 20, 80, (40, 60))
        sites = st.multiselect('Avoid Restriction Sites:', ['EcoRI', 'BamHI', 'HindIII', 'NcoI', 'XhoI', 'NdeI', 'SalI', 'XbaI', 'SpeI', 'PstI'], default=[])
        avoid_sigs = st.checkbox("Avoid Cryptic Signals", value=False)
        ref, n_cand = '', 5
        if 'Harmonization' in strat: ref = st.text_input('Reference GenBank Accession ID:')
        elif 'Probabilistic' in strat: n_cand = st.number_input('Number of Candidates:', 1, 20, 5)

with tab2:
    st.markdown('<div class="section-header">EXECUTION LOG</div>', unsafe_allow_html=True)
    if st.session_state.optimization_complete:
        st.markdown(f'<div class="terminal-box">{"<br>".join(st.session_state.optimization_logs)}</div>', unsafe_allow_html=True)
    elif not run_btn: st.info("Run optimization to generate logs.")

if run_btn:
    st.session_state.optimization_complete, st.session_state.optimization_logs = False, []
    with tab2:
        log = st.empty()
        def add_log(m):
            t = time.strftime('%M:%S'); st.session_state.optimization_logs.append(f'[{t}] {m}')
            log.markdown(f'<div class="terminal-box">{"<br>".join(st.session_state.optimization_logs)}</div>', unsafe_allow_html=True)
        
        add_log("Initializing engine...")
        org_data = freq_db.get(sel_org, list(freq_db.values())[0])
        ofreq, tid = org_data['frequencies'], org_data['table_id']
        
        if 'Frequency' in strat: 
            add_log("Executing Codon Frequency Replacement..."); dna = optimizer.optimize_frequency(clean_seq, ofreq, gc_min, gc_max, sites, tid)
        elif 'Harmonization' in strat: 
            add_log(f"Executing Codon Harmonization against {ref}..."); dna = optimizer.optimize_harmonization(clean_seq, ofreq, ref, tid)
        else:
            add_log(f"Executing Probabilistic Sampling ({n_cand} candidates)...")
            cands = optimizer.optimize_probabilistic(clean_seq, ofreq, n_cand, gc_min, gc_max, sites, tid, avoid_sigs)
            dna = cands[0][0]
        
        add_log("Scanning optimized sequence...")
        st.session_state.found_sites = optimizer.scan_restriction_sites(dna, sites)
        st.session_state.cryptic_signals = optimizer.scan_cryptic_signals(dna)
        add_log("Calculating metrics...")
        mfe = optimizer.calculate_mfe(dna)
        cai = optimizer.calculate_cai(dna, ofreq, tid)
        gc_act = (dna.count('G')+dna.count('C'))/len(dna)*100
        st.session_state.optimized_sequence = dna
        st.session_state.optimization_metrics = {'CAI':round(cai,3), 'GC%':round(gc_act,1), 'LENGTH':len(dna), 'MFE':mfe}
        st.session_state.optimization_complete = True
        st.rerun()

if st.session_state.optimization_complete:
    dna, m = st.session_state.optimized_sequence, st.session_state.optimization_metrics
    org_data = freq_db.get(sel_org, list(freq_db.values())[0])
    ofreq, tid = org_data['frequencies'], org_data['table_id']
    
    with tab3:
        cols = st.columns(4)
        for i, (k, v) in enumerate(m.items()): cols[i].metric(k, v)
        
        codons = [dna[i:i+3] for i in range(0, len(dna), 3) if len(dna[i:i+3]) == 3]
        c_counts = {c: codons.count(c) for c in set(codons)}
        df_bar = pd.DataFrame({'Codon': list(c_counts.keys()), 'Frequency': list(c_counts.values())})
        fig1 = px.bar(df_bar, x='Codon', y='Frequency', title='Optimized Codon Usage').update_layout(template='plotly_dark', plot_bgcolor='#0a0d14', paper_bgcolor='#0a0d14', font_color='#00f5ff')
        st.plotly_chart(fig1, use_container_width=True)
        
        c1, c2 = st.columns(2)
        with c1:
            st.markdown('<div class="section-header">SLIDING WINDOW CAI</div>', unsafe_allow_html=True)
            cais = optimizer.sliding_window_cai(dna, ofreq, 10, tid)
            fig2 = px.line(y=cais).update_layout(template='plotly_dark', plot_bgcolor='#0a0d14', paper_bgcolor='#0a0d14', font_color='#39ff14')
            st.plotly_chart(fig2, use_container_width=True)
        with c2:
            st.markdown('<div class="section-header">SLIDING WINDOW GC%</div>', unsafe_allow_html=True)
            gcs = optimizer.sliding_window_gc(dna)
            fig3 = px.line(y=gcs).update_layout(template='plotly_dark', plot_bgcolor='#0a0d14', paper_bgcolor='#0a0d14', font_color='#00f5ff')
            st.plotly_chart(fig3, use_container_width=True)
            
        st.markdown('<div class="section-header">AMINO ACID COMPOSITION</div>', unsafe_allow_html=True)
        aa_comp = {aa: clean_seq.count(aa) for aa in set(clean_seq)}
        df_aa = pd.DataFrame({'AA': list(aa_comp.keys()), 'Count': list(aa_comp.values())}).sort_values('Count')
        fig4 = px.bar(df_aa, x='Count', y='AA', orientation='h').update_layout(template='plotly_dark', plot_bgcolor='#0a0d14', paper_bgcolor='#0a0d14', font_color='#39ff14')
        st.plotly_chart(fig4, use_container_width=True)

        st.markdown('<div class="section-header">RARE CODON MAP</div>', unsafe_allow_html=True)
        cells = ''
        for i in range(0, len(dna), 3):
            c = dna[i:i + 3]
            f = ofreq.get(c, 0)
            color = '#39ff14' if f > 0.2 else '#ffb300' if f > 0.05 else '#ff3b5c'
            cells += f'<div style="width:10px; height:20px; background-color:{color}; display:inline-block; margin:1px;" title="{c}: {f:.3f}"></div>'
        st.markdown(f'<div style="background:rgba(0,0,0,0.5); padding:10px; border-radius:5px;">{cells}</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-header">SEQUENCE INTEGRITY & ISSUES</div>', unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        with c1:
            st.subheader("Restriction Sites")
            for s, p in st.session_state.found_sites.items():
                if p: st.warning(f"{s} found at {p}")
                else: st.success(f"{s} successfully avoided")
        with c2:
            st.subheader("Cryptic Signals")
            for s, p in st.session_state.cryptic_signals.items():
                if p: st.warning(f"{s} found at {p}")
                else: st.success(f"{s} successfully avoided")

        st.markdown('<div class="section-header">INTEGRATED SITE MAP</div>', unsafe_allow_html=True)
        fig6 = go.Figure()
        fig6.add_trace(go.Scatter(x=[0, len(dna)], y=[0, 0], mode='lines', line=dict(color='#64748b', width=4), name='Backbone'))
        for s, pos in st.session_state.found_sites.items():
            if pos: fig6.add_trace(go.Scatter(x=pos, y=[2]*len(pos), mode='markers+text', marker=dict(symbol='triangle-down', size=15, color='#ff3b5c'), text=[s]*len(pos), textposition='top center', name=s))
        for s, pos in st.session_state.cryptic_signals.items():
            if pos: fig6.add_trace(go.Scatter(x=pos, y=[1.5]*len(pos), mode='markers+text', marker=dict(symbol='circle', size=12, color='#ffb300'), text=[s]*len(pos), textposition='top center', name=s))
        fig6.update_layout(template='plotly_dark', plot_bgcolor='#0a0d14', paper_bgcolor='#0a0d14', yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, 3]), height=400)
        st.plotly_chart(fig6, use_container_width=True)

    with tab4:
        st.markdown('<div class="section-header">IN-SILICO VECTOR INTEGRATION</div>', unsafe_allow_html=True)
        c1, c2 = st.columns([1, 2])
        with c1:
            st.info("Gibson Assembly / Restriction Cloning Simulation.")
            cat = st.selectbox("Category", sorted(vector_integration.VECTORS.keys()))
            vec = st.selectbox("Plasmid", sorted(vector_integration.VECTORS[cat].keys()))
            if st.button("CONSTRUCT PLASMID", type='primary', use_container_width=True):
                with st.spinner("Processing..."):
                    fseq, rec = vector_integration.assemble_vector(cat, vec, dna)
                    st.session_state.map = vector_integration.render_vector_map(rec)
                    st.session_state.zoom = vector_integration.render_sequence_zoom(rec, vector_integration.VECTORS[cat][vec]["insertion_site"])
                    st.session_state.vector_final_seq = fseq
        with c2:
            if 'map' in st.session_state:
                st.markdown("##### Global Map")
                st.image(st.session_state.map, use_column_width=True)
                st.markdown("##### Zoomed View")
                st.image(st.session_state.zoom, use_column_width=True)
                st.download_button("DOWNLOAD FULL VECTOR FASTA", st.session_state.vector_final_seq, f"{vec}_opt.fasta", use_container_width=True)

    with tab5:
        st.markdown('<div class="section-header">EXPORT ASSETS</div>', unsafe_allow_html=True)
        e1, e2 = st.columns(2)
        with e1:
            st.download_button('DOWNLOAD OPTIMIZED FASTA', dna, 'optimized.fasta', use_container_width=True)
            b64 = base64.b64encode(dna.encode()).decode()
            copy_html = f"<button onclick=\"navigator.clipboard.writeText(atob('{b64}'))\" style=\"width:100%; padding:8px; background:transparent; border:1px solid #00f5ff; color:#00f5ff; font-family:'Rajdhani', sans-serif; cursor:pointer;\">COPY TO CLIPBOARD</button>"
            components.html(copy_html, height=50)
        with e2:
            if st.button('PREPARE PDF REPORT', use_container_width=True):
                with st.spinner("Generating..."):
                    figs = [fig1, fig2, fig3, fig4, fig6]
                    path = report_generator.generate_pdf_report(clean_seq, dna, m, figs, sel_org, None)
                    with open(path, 'rb') as f: st.download_button('DOWNLOAD PDF NOW', f, 'report.pdf', use_container_width=True, type='primary')
