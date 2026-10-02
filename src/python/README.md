# Python scripts — prime-gap autocorrelation

## Requirements

    pip install -r requirements.txt

## Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `fit_c.py` | Weighted fit of $c$ in $|\rho_1| \sim (\log p)^{-c}$ (§7.2) | `correlogram_v9_1.csv` | `fit_c_summary.txt` |
| `compute_psi.py` | Direct `theta(P)` by prime enumeration, `psi(10^9)-10^9` (Appendix A, Route 2) | — | console report |
| `delta_mean_verify.py` | Telescoping identity, $\psi(10^9)-10^9$ (Appendix A) | `theta_powers_terms.csv` | `delta_mean_verification.csv` |
| `make_figures.py` | Figures 1–5 | 4 CSV files | `fig1.pdf` … `fig5.pdf` |
| `make_appendix_b.py` | Synthetic non-stationarity control (Appendix B) | — | `appendix_b_table.tex`, `appendix_b_acf.csv` |
| `hl_pairs_v9_1.py` | HL + IE pair-ratio prototype (§8) | — | `hl_pairs_v9_1_python.csv` |

## Reproduction order

    cd src/python
    python fit_c.py
    python delta_mean_verify.py
python compute_psi.py
    python make_figures.py
    python make_appendix_b.py
    python hl_pairs_v9_1.py

## Note on `hl_pairs_v9_1.py`

The published Table 9 numbers come from the C++ pipeline
(`task9_ie_model` in `wheel_analysis_v9_1_complete.cpp`). This Python
script is a pedagogical, self-contained prototype. Small numerical
differences (<0.1%) may occur.