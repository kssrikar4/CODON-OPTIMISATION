import os, tempfile, textwrap
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor

def generate_pdf_report(aa, dna, metrics, figs, org, organ, out='report.pdf'):
    c = canvas.Canvas(out, pagesize=letter)
    w, h = letter
    def page(): c.setFillColor(HexColor('#0a0d14')); c.rect(0, 0, w, h, fill=1)
    page()
    c.setFillColor(HexColor('#00f5ff')); c.setFont('Courier-Bold', 24); c.drawString(50, h-50, 'CODON OPTIMISATION: REPORT')
    c.setFont('Courier', 12); c.setFillColor(HexColor('#e2e8f0'))
    y = h-90
    for txt in [f'Target: {org}', f'Length: {len(aa)} AA', f"CAI: {metrics.get('cai','N/A')}", f"GC: {metrics.get('gc','N/A')}%"]:
        c.drawString(50, y, txt); y -= 20
    if 'mfe' in metrics: c.drawString(50, y, f"MFE: {metrics.get('mfe')} kcal/mol"); y -= 20
    y -= 20; c.setFillColor(HexColor('#39ff14')); c.drawString(50, y, 'Optimized DNA:'); y -= 20
    c.setFillColor(HexColor('#e2e8f0')); c.setFont('Courier', 10)
    for line in textwrap.wrap(dna, width=80):
        if y < 50: c.showPage(); page(); c.setFillColor(HexColor('#e2e8f0')); c.setFont('Courier', 10); y = h-50
        c.drawString(50, y, line); y -= 12
    y -= 30
    with tempfile.TemporaryDirectory() as tmp:
        for i, fig in enumerate(figs):
            img = os.path.join(tmp, f"{i}.png")
            try:
                fig.write_image(img, width=800, height=400, scale=2)
                if y < 250: c.showPage(); page(); y = h-50
                c.drawImage(img, 50, y-220, width=500, height=220, preserveAspectRatio=True); y -= 240
            except: pass
    c.save(); return out
