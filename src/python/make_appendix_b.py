# make_appendix_b.py
# Appendix B: Synthetic demonstration of the non-stationarity artifact.
# The noise variance is tuned so that max_raw / threshold ≈ 7,
# matching the real data in Table 4c.

import numpy as np
import pandas as pd
from math import sqrt

try:
    from scipy.stats import norm
    z = norm.ppf(1 - 0.05 / (2 * 1024))   # ≈ 4.0612
except ImportError:
    z = 4.0612


def acf_global(x, K):
    """Global sample autocorrelation, lags 1..K, global mean/variance."""
    x = np.asarray(x, dtype=float)
    m = len(x)
    xc = x - x.mean()
    denom = np.dot(xc, xc) / m
    out = np.empty(K)
    for k in range(1, K + 1):
        num = np.dot(xc[:-k], xc[k:]) / (m - k)
        out[k - 1] = num / denom
    return out


# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
n         = 1_000_000      # synthetic sequence length
K         = 1024           # number of lags
seed      = 42
sigma_gap = 6.0            # tuned so that max_raw / threshold ≈ 7 (cf. Table 4c)

rng = np.random.default_rng(seed)

# Synthetic prime positions ~ n log n, so log p_n ~ log n + log log n
idx  = np.arange(2, n + 2)
logp = np.log(idx * np.log(idx))

# Non-stationary sequence: deterministic trend + i.i.d. noise
eta   = rng.normal(0.0, sigma_gap, size=n)
h_raw = logp + eta         # non-stationary (trend + noise)
h_det = h_raw - logp       # exactly detrended => pure i.i.d. noise

raw_acf = acf_global(h_raw, K)
det_acf = acf_global(h_det, K)

se  = 1.0 / sqrt(n)
thr = z * se


def stats(acf):
    pos  = int(np.sum(acf > 0))
    fp   = int(np.argmax(acf > 0) + 1) if pos > 0 else -1
    mx   = float(acf.max())
    mk   = int(np.argmax(acf) + 1)
    tail = float(acf[159:].mean())   # mean over k in [160, K]
    return dict(pos=pos, fp=fp, mx=mx, mk=mk, ratio=mx / thr, tail=tail)


rs = stats(raw_acf)
ds = stats(det_acf)

# ----------------------------------------------------------------------
# Save CSV of full ACFs
# ----------------------------------------------------------------------
pd.DataFrame({
    'k': np.arange(1, K + 1),
    'raw': raw_acf,
    'detrended': det_acf,
}).to_csv('appendix_b_acf.csv', index=False)

# ----------------------------------------------------------------------
# Write LaTeX table
# ----------------------------------------------------------------------
with open('appendix_b_table.tex', 'w', encoding='utf-8') as f:
    f.write(r'\begin{table}[h]\centering' + '\n')
    f.write(r'\begin{tabular}{lcc}\toprule' + '\n')
    f.write(r'Quantity & Raw $h_n$ & Detrended $\tilde h_n$ \\ \midrule' + '\n')
    f.write(f'Total lags & {K} & {K} \\\\\n')
    f.write(f'Positive lags & {rs["pos"]} & {ds["pos"]} \\\\\n')
    f.write(f'First positive lag & {rs["fp"]} & {ds["fp"]} \\\\\n')
    f.write(f'Maximum positive & {rs["mx"]:.3e} & {ds["mx"]:.3e} \\\\\n')
    f.write(f'Bonferroni threshold & {thr:.3e} & {thr:.3e} \\\\\n')
    f.write(f'Max / threshold & {rs["ratio"]:.3f} & {ds["ratio"]:.3f} \\\\\n')
    f.write(f'Tail mean ($k\\in[160,{K}]$) & {rs["tail"]:.3e} & {ds["tail"]:.3e} \\\\\n')
    f.write(r'\bottomrule\end{tabular}' + '\n')
    f.write(
        r'\caption{Synthetic non-stationarity control. The raw sequence '
        r'$h_n = \log p_n + \eta_n$ with $\eta_n \sim \mathcal N(0, 6^2)$ '
        r'i.i.d.\ exhibits a positive plateau whose maximum exceeds the '
        r'Bonferroni threshold by a factor comparable to the raw prime-gap '
        r'correlogram (Table~4c, ratio $= 7.45$). After exact detrending '
        r'$\tilde h_n = h_n - \log p_n$, the plateau disappears entirely. '
        r'This confirms that the raw positive plateau is a generic '
        r'consequence of non-stationarity, not a specific property of the '
        r'prime gaps.}' + '\n'
    )
    f.write(r'\label{tab:B1}\end{table}' + '\n')

# ----------------------------------------------------------------------
# Console summary
# ----------------------------------------------------------------------
print(f'raw : pos={rs["pos"]:4d}/{K}  '
      f'first_pos={rs["fp"]:4d}  '
      f'max={rs["mx"]:+.3e}  '
      f'max/thr={rs["ratio"]:.3f}  '
      f'tail={rs["tail"]:+.3e}')
print(f'det : pos={ds["pos"]:4d}/{K}  '
      f'first_pos={ds["fp"]:4d}  '
      f'max={ds["mx"]:+.3e}  '
      f'max/thr={ds["ratio"]:.3f}  '
      f'tail={ds["tail"]:+.3e}')