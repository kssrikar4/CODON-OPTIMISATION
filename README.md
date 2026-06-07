# CODON OPTIMISATION

![Python](https://img.shields.io/badge/python-3.10%2B-blue?style=flat-square&logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/Streamlit-1.33.0-FF4B4B?style=flat-square&logo=streamlit&logoColor=white)
![Polars](https://img.shields.io/badge/Polars-Data%20Engine-CD7F32?style=flat-square&logo=polars&logoColor=white)
![Plotly](https://img.shields.io/badge/Plotly-5.20.0-3F4F75?style=flat-square&logo=plotly&logoColor=white)
![Biopython](https://img.shields.io/badge/Biopython-1.83-blueviolet?style=flat-square)
![SnapGene Reader](https://img.shields.io/badge/SnapGene%20Reader-0.1.23-00A896?style=flat-square)
![seqfold](https://img.shields.io/badge/seqfold-0.10.0-orange?style=flat-square)
![CUDA](https://img.shields.io/badge/CUDA-Accelerated-76B900?style=flat-square&logo=nvidia&logoColor=white)

A high-performance sequence engineering platform for synthetic biology. Codon Optimisation provides a streamlined pipeline for advanced codon optimization, plasmid mapping, and sequence integrity analysis.

## 1. Architecture

- **Optimization Engine (`optimizer.py`)**: Implements deterministic and stochastic optimization algorithms with heuristic constraints (GC content, restriction sites, cryptic signals).
- **Vector Integration Engine (`vector_integration.py`)**: A geometric rendering system for plasmid maps featuring iterative label collision avoidance and dynamic MCS detection.
- **Data Pipeline (`codon_loader.py`)**: A Polars-based ingestion system for CoCoPUTs / HIVE-CUTs databases, supporting thousands of taxonomic and tissue-specific usage profiles.
- **Reporting System (`report_generator.py`)**: A minimalist PDF generator.

## 2. Optimization Strategies

### Codon Frequency Replacement (CFR)
Replaces every amino acid with the codon that has the highest relative frequency in the target organism's usage table. This maximizes translational speed but increases the risk of restriction site motifs and mRNA secondary structure stability issues.

### Codon Harmonization
Matches the codon usage profile of the target organism to that of a reference gene (via GenBank ID). This strategy aims to preserve the translational "rhythm" and co-translational folding of the protein by using rare codons where the source organism uses them.

### Probabilistic Sampling (Recommended)
Uses a stochastic sampling algorithm to select codons proportional to their frequency. It iteratively validates candidates against GC constraints, restriction site filters, and cryptic signal scanners, selecting the candidate with the highest Codon Adaptation Index (CAI).

## 3. Sequence Integrity Analysis

### Restriction Site Scanner
Automated detection of standard cloning motifs (EcoRI, BamHI, HindIII, etc.). The system provides explicit "Success" or "Warning" feedback for every selected site, ensuring no unintended motifs remain in the optimized DNA.

### Cryptic Signal Filter
Scans for eukaryotic and prokaryotic regulatory signals that can interfere with expression:
- **PolyA Signals**: AATAAA, ATTAAA
- **Splice Donors**: CAGGTAAGT
- **Splice Acceptors**: TTCAG[GA]

### mRNA Secondary Structure (MFE)
Calculates the Minimum Free Energy (MFE) of the mRNA transcript using a Nussinov-based folding algorithm. For long sequences (>5000bp), the engine supports CUDA-accelerated computation if a compatible GPU is detected.

## 4. Vector Integration & Mapping

The engine dynamically loads SnapGene (.dna) files from `data/vectors/` and performs in-silico cloning simulation.

### Geometric Layout Engine
To handle dense feature sets (ORFs, primers, promoters), the mapping engine implements a **Vertical Spreading Algorithm**:
1. **Hemisphere Segregation**: Labels are split into left/right hemispheres to prevent cross-path leader lines.
2. **Radial Ordering**: Features are sorted by genome position; leader lines are routed straight to prevent intersections.
3. **Collision Avoidance**: Overlapping label bounding boxes are iteratively pushed apart along the Y-axis until a minimum vertical padding is reached.
4. **Visibility Logic**: Labels are rendered with a 2px black stroke/outline to ensure legibility over high-contrast features.

## 5. Setup & Usage

### Prerequisites
- Python 3.10+

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/kssrikar4/CODON-OPTIMISATION.git
    cd CODON-OPTIMISATION
    ```
    
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Download usage data from [FDA DNA HIVE CoCoPUTs](https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs).
4. Place `.tsv` files in `data/cocoputs/`.
5. Add SnapGene `.dna` files to categorized subfolders in `data/vectors/`.

### Execution
Run the local server:
```bash
streamlit run app.py
```

## 6. Data Credits
- **Usage Profiles**: CoCoPUTs / HIVE-CUTs (FDA / NCBI).
- **Vector Sequences**: Sourced from [SnapGene Plasmid Files](https://www.snapgene.com/plasmids).

### Acknowledgements

I express my sincere gratitude to the developers and communities behind the core open-source libraries and platforms that made this project possible:

* **Streamlit:** For providing the framework to build this interactive web application entirely in Python, streamlining UI development.
  * [Documentation](https://docs.streamlit.io/) | [GitHub](https://github.com/streamlit/streamlit)

* **Plotly & Matplotlib:** For their graphing and visualization libraries that power the high-quality, interactive rendering engine of this application.
  * [Plotly Python](https://plotly.com/python/) | [Matplotlib](https://matplotlib.org/)

* **Biopython & Domain-Specific Tooling:** For the computational biology utilities, sequence manipulation, secondary structure prediction (`seqfold`), and visual mapping tools (`dna_features_viewer`, `snapgene-reader`).
  * [Biopython Resource](https://biopython.org/) | [DNA Features Viewer GitHub](https://github.com/EdJoJob/DnaFeaturesViewer)

* **Pandas & Polars:** For high-performance data manipulation, structured dataframe analysis, and data handling backend capabilities.
  * [Pandas](https://pandas.pydata.org) | [Polars Documentation](https://pola.rs)

* **ReportLab & FPDF2:** For enabling programmatic PDF generation and reporting functionality.

* **Google Gemini:** For invaluable assistance in code analysis, optimization, and refining the architectural development process of this project.
  * [Gemini](https://geminicli.com)

Their contributions to the open-source ecosystem have made the development of `Codon Optimisation` feasible.

[Mozilla Public License Version 2.0](LICENSE) - Feel free to use and modify
