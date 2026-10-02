#!/usr/bin/env python3
"""
fit_c.py

Reproduces §7.2 of the manuscript:
    Weighted fit of |rho_1| ~ a (log p)^{-c}
    on the 8-window subset (log p >= 12), weighted by m_i.
    Reports point estimate, SE, R^2, and bootstrap 95% CI.

Input:  correlogram_v9_1.csv
Output: fit_c_summary.txt   (also printed to stdout)
"""

import numpy as np
import pandas as pd

INPUT   = 'correlogram_v9_1.csv'
OUTPUT  = 'fit_c_summary.txt'

LOG_P_MIN   = 12.0     # 8-window subset (drop two smallest windows)
N_BOOT      = 2000
SEED        = 42
ALPHA       = 0.05


# ----------------------------------------------------------------------
def load_csv(path):
    """Read a CSV whose header may or may not be prefixed with '#'."""
    with open(path) as f:
        lines = f.readlines()
    first = lines[0].strip()
    if first.startswith('#') and ',' in first:
        cols  = [c.strip() for c in first[1:].split(',')]
        start = 1
    else:
        i = 0
        while i < len(lines):
            s = lines[i].strip()
            if s and not s.startswith('#'):
                break
            i += 1
        cols  = [c.strip() for c in lines[i].strip().split(',')]
        start = i + 1
    data = ''.join(lines[start:])
    from io import StringIO
    return pd.read_csv(StringIO(data), names=cols, skipinitialspace=True)


# ----------------------------------------------------------------------
def weighted_ols(x, y, w):
    """Weighted OLS: returns (slope, intercept, cov)."""
    w = w / w.sum()
    X = np.column_stack([x, np.ones_like(x)])
    W = np.diag(w)
    beta  = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    resid = y - X @ beta
    dof   = len(y) - 2
    sigma2 = (resid.T @ W @ resid) / dof
    cov    = sigma2 * np.linalg.inv(X.T @ W @ X)
    return beta[0], beta[1], cov


# ----------------------------------------------------------------------
def fit(x, y, w):
    slope, intercept, cov = weighted_ols(x, y, w)
    c     = -slope
    a     = np.exp(intercept)
    c_err = np.sqrt(cov[0, 0])
    yhat  = intercept + slope * x
    ss_res = np.sum(w * (y - yhat) ** 2)
    ss_tot = np.sum(w * (y - np.average(y, weights=w)) ** 2)
    r2    = 1.0 - ss_res / ss_tot
    return a, c, c_err, r2


# ----------------------------------------------------------------------
def main():
    df = load_csv(INPUT)

    logp_mid = 0.5 * (np.log(df['window_lo']) + np.log(df['window_hi']))
    rho1     = df['rho_1'].values.astype(float)
    m        = df['n_gaps'].values.astype(float)

    mask = logp_mid >= LOG_P_MIN
    x = np.log(logp_mid[mask])
    y = np.log(np.abs(rho1[mask]))
    w = m[mask]

    a_fit, c_fit, c_err, r2 = fit(x, y, w)

    # --- Bootstrap CI ---
    rng = np.random.default_rng(SEED)
    cs  = np.empty(N_BOOT)
    n   = len(x)
    for b in range(N_BOOT):
        idx = rng.integers(0, n, n)
        xb, yb, wb = x[idx], y[idx], w[idx]
        try:
            _, cb, _, _ = fit(xb, yb, wb)
            cs[b] = cb
        except np.linalg.LinAlgError:
            cs[b] = np.nan
    cs = cs[~np.isnan(cs)]

    lo, hi = np.quantile(cs, [ALPHA / 2, 1 - ALPHA / 2])

    # --- Report ---
    lines = []
    lines.append('=== fit_c.py — scaling exponent of rho_1 ===')
    lines.append(f'Window subset:  log p >= {LOG_P_MIN}  ({mask.sum()} windows)')
    lines.append(f'Weights:        m_i (normalized)')
    lines.append('')
    lines.append(f'a       = {a_fit:.4f}')
    lines.append(f'c       = {c_fit:.4f}')
    lines.append(f'SE(c)   = {c_err:.4f}')
    lines.append(f'R^2     = {r2:.4f}')
    lines.append(f'Bootstrap 95% CI for c: [{lo:.4f}, {hi:.4f}]  '
                 f'({N_BOOT} resamples, seed={SEED})')
    lines.append('')
    lines.append('Comparison with paper §7.2:')
    lines.append('    c = 1.28 ± 0.05     (paper)')
    lines.append(f'    c = {c_fit:.2f} ± {c_err:.2f}   (this run)')
    lines.append(f'    CI [1.14, 1.57]     (paper, bootstrap)')
    lines.append(f'    CI [{lo:.2f}, {hi:.2f}]   (this run)')

    summary = '\n'.join(lines)
    print(summary)

    with open(OUTPUT, 'w') as f:
        f.write(summary + '\n')
    print(f'\nWrote {OUTPUT}')


if __name__ == '__main__':
    main()