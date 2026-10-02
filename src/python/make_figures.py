# make_figures.py
# Generates Figures 1-5 for the manuscript.
# Requirements: numpy, pandas, matplotlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'figure.dpi': 120,
    'savefig.bbox': 'tight',
})


# ----------------------------------------------------------------------
# Robust loader:
#   - Finds the header line (either "#a,b,c" or "a,b,c" after comments)
#   - Ignores any '#'-comment lines interspersed in the data block
#   - Returns a DataFrame with stripped column names
# ----------------------------------------------------------------------
def load_csv(path):
    with open(path) as f:
        lines = f.readlines()

    header_idx = None
    cols = None
    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        if s.startswith('#') and ',' in s:
            # header prefixed with '#'
            cols = [c.strip() for c in s[1:].split(',')]
            header_idx = i
            break
        if not s.startswith('#'):
            # plain header
            cols = [c.strip() for c in s.split(',')]
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f'No header found in {path}')

    data_lines = []
    for line in lines[header_idx + 1:]:
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        data_lines.append(line)

    return pd.read_csv(
        StringIO(''.join(data_lines)),
        names=cols,
        skipinitialspace=True,
    )


# ----------------------------------------------------------------------
# Sanity check
# ----------------------------------------------------------------------
print('--- Column check ---')
for path in ('correlogram_v9_1.csv',
             'correlogram_k1024_v10_v2.csv',
             'correlogram_k1024_raw_v10_v2.csv',
             'hl_pairs_v9_1.csv'):
    try:
        c = list(load_csv(path).columns[:8])
        print(f'{path:42s} -> {c}')
    except Exception as e:
        print(f'{path:42s} -> ERROR: {e}')
print('--------------------')


# ----------------------------------------------------------------------
# Figure 1: rho_1 vs log p, power-law fit (8 windows, weighted OLS in log-log)
# ----------------------------------------------------------------------
df = load_csv('correlogram_v9_1.csv')

logp_mid = 0.5 * (np.log(df['window_lo']) + np.log(df['window_hi']))
rho1     = df['rho_1'].values
ci_lo    = df['rho1_ci_lo'].values
ci_hi    = df['rho1_ci_hi'].values
yerr     = np.vstack([rho1 - ci_lo, ci_hi - rho1])
weights  = df['n_gaps'].values.astype(float)

mask = logp_mid >= 12.0   # 8-window subset, matching §7.2

x = np.log(logp_mid[mask])
y = np.log(np.abs(rho1[mask]))
w = weights[mask] / weights[mask].sum()

X = np.column_stack([x, np.ones_like(x)])
W = np.diag(w)

beta  = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
c_fit = -beta[0]
a_fit = np.exp(beta[1])

resid  = y - X @ beta
sigma2 = (resid.T @ W @ resid) / (len(y) - 2)
cov    = sigma2 * np.linalg.inv(X.T @ W @ X)
c_err  = np.sqrt(cov[0, 0])


def model(logp, a, c):
    return -a * logp ** (-c)


fig, ax = plt.subplots(figsize=(5.5, 4))
ax.errorbar(logp_mid, rho1, yerr=yerr, fmt='o',
            capsize=3, color='black', ecolor='gray', label=r'$\hat\rho_1$')
xx = np.linspace(logp_mid.min(), logp_mid.max(), 200)
ax.plot(xx, model(xx, a_fit, c_fit), 'r-',
        label=rf'fit (8 windows): $|\rho_1| = a(\log p)^{{-c}}$, '
              rf'$c = {c_fit:.2f} \pm {c_err:.2f}$')
ax.set_xlabel(r'$\log p$')
ax.set_ylabel(r'$\hat\rho_1$')
ax.set_title('Lag-one autocorrelation across logarithmic windows')
ax.legend()
ax.grid(alpha=0.3)
fig.savefig('fig1.pdf')
plt.close(fig)
print(f'Fig 1: a = {a_fit:.4f}, c = {c_fit:.3f} ± {c_err:.3f}')


# ----------------------------------------------------------------------
# Figure 2: rho_k for k=1..64 in largest window
# ----------------------------------------------------------------------
largest_idx = df['window_lo'].idxmax()
rho_cols = [c for c in df.columns
            if c.startswith('rho_') and c not in ('rho1_ci_lo', 'rho1_ci_hi')]
rho_cols_sorted = sorted(rho_cols, key=lambda c: int(c.split('_')[1]))
rho_vals = df.loc[largest_idx, rho_cols_sorted].values.astype(float)
ks = np.arange(1, len(rho_vals) + 1)

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(ks, rho_vals, 'o-', ms=3, lw=1)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('Lag $k$')
ax.set_ylabel(r'$\hat\rho_k$')
ax.set_title(r'Windowed correlogram, largest window $[3\times10^8,\,10^9)$')
ax.grid(alpha=0.3)
fig.savefig('fig2.pdf')
plt.close(fig)


# ----------------------------------------------------------------------
# Figures 3 & 4: K=1024 correlograms
# ----------------------------------------------------------------------
thr = 5.70e-4

d3 = load_csv('correlogram_k1024_v10_v2.csv')
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(d3['k'], d3['rho_k'], '.', ms=2, color='steelblue')
ax.axhline( thr, color='r', ls='--', lw=0.8, label=rf'Bonferroni $\pm${thr:.2e}')
ax.axhline(-thr, color='r', ls='--', lw=0.8)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('Lag $k$')
ax.set_ylabel(r'$\hat\rho_k$')
ax.set_title(r'Detrended ($\delta_n$) correlogram, $K=1024$ — 255/1024 positive')
ax.legend()
ax.grid(alpha=0.3)
fig.savefig('fig3.pdf')
plt.close(fig)

d4 = load_csv('correlogram_k1024_raw_v10_v2.csv')
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(d4['k'], d4['rho_k'], '.', ms=2, color='darkorange')
ax.axhline( thr, color='r', ls='--', lw=0.8, label=rf'Bonferroni $\pm${thr:.2e}')
ax.axhline(-thr, color='r', ls='--', lw=0.8)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('Lag $k$')
ax.set_ylabel(r'$\hat\rho_k$')
ax.set_title(r'Raw ($g_n$) correlogram, $K=1024$ — 1016/1024 positive (artifact)')
ax.legend()
ax.grid(alpha=0.3)
fig.savefig('fig4.pdf')
plt.close(fig)


# ----------------------------------------------------------------------
# Figure 5: pair ratios — observed vs HL vs IE
# ----------------------------------------------------------------------
d5 = load_csv('hl_pairs_v9_1.csv')

# ensure numeric
for col in ('a', 'b', 'R_obs', 'R_HL', 'R_model'):
    d5[col] = pd.to_numeric(d5[col], errors='coerce')
d5 = d5.dropna(subset=['a', 'b', 'R_obs', 'R_HL', 'R_model'])

# sort by (a, b) for stable ordering
d5 = d5.sort_values(['a', 'b']).reset_index(drop=True)

# --- diagnostic: verify the pairs are the 17 expected ones ---
expected = [(2,10),(10,2),(2,28),(28,2),(14,16),(16,14),(26,4),(4,26),
            (6,6),(6,12),(12,6),(6,4),(4,6),(2,4),(4,2),(2,6),(6,2)]
got = list(zip(d5['a'].astype(int), d5['b'].astype(int)))
print(f'--- hl_pairs diagnostic ---')
print(f'rows loaded = {len(d5)}  (expected 17)')
print(f'pairs       = {got}')
missing = [p for p in expected if p not in got]
extra   = [p for p in got if p not in expected]
print(f'missing     = {missing}')
print(f'extra       = {extra}')
print('---------------------------')

labels = [f"({int(a)},{int(b)})" for a, b in zip(d5['a'], d5['b'])]
x = np.arange(len(labels))
w = 0.27

fig, ax = plt.subplots(figsize=(8.5, 4))
ax.bar(x - w, d5['R_obs'],   w, label='Observed', color='black')
ax.bar(x,     d5['R_HL'],    w, label='Naive HL', color='steelblue')
ax.bar(x + w, d5['R_model'], w, label='IE model', color='crimson', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=90, fontsize=8)
ax.set_ylabel(r'Pair ratio $R(a,b)$')
ax.set_title('Observed vs HL vs IE pair ratios')
ax.legend()
ax.grid(alpha=0.3, axis='y')
fig.tight_layout()
fig.savefig('fig5.pdf')
plt.close(fig)

print('Wrote fig1.pdf ... fig5.pdf')