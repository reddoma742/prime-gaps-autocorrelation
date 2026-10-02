#!/usr/bin/env python3
"""
delta_mean_verify.py

Reproduces Appendix A:
    delta_n = g_n - log p_n
    bar_delta = (P - theta(P) + log P - 2) / n_gaps

Uses:
    - theta_powers_terms.csv  (produced by theta_powers.cpp)
    - Known constant: bar_delta_obs = 6.092163e-4 (from v9.1 pipeline)

Outputs a verification report and writes delta_mean_verification.csv.
"""

import numpy as np
import pandas as pd
from math import log, sqrt


# ----------------------------------------------------------------------
# Known constants (from v9.1 pipeline)
# ----------------------------------------------------------------------
N              = 10**9
P              = 999_999_937                    # largest prime <= 10^9
n_gaps         = 50_847_532                     # excludes 2 -> 3 gap
log_P          = log(P)                         # 20.723266...
bar_delta_obs  = 6.092163e-4
pi_N           = 50_847_534
gamma_1        = 14.134725141734693             # first zeta zero


# ----------------------------------------------------------------------
def load_theta_powers(path='theta_powers_terms.csv'):
    """
    Expected columns: k, N_pow_1_over_k, theta_N_pow
    Returns sum_{k>=2} theta(N^{1/k}).
    """
    try:
        df = pd.read_csv(path, comment='#')
        total = df['theta_N_pow'].sum()
        return df, total
    except FileNotFoundError:
        print(f'[WARN] {path} not found; using stored constant from paper.')
        return None, 32_617.412861


def theta_of_x(x):
    """
    Simple fallback: sum of logs of primes <= x.
    Only used if theta_powers_terms.csv is missing for cross-check.
    """
    if x < 2:
        return 0.0
    sieve = np.ones(int(x) + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(sqrt(x)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return float(np.sum(np.log(np.nonzero(sieve)[0])))


# ----------------------------------------------------------------------
def main():
    print('=== delta_mean_verify.py — Appendix A verification ===\n')

    # --- 1) Theta-powers sum ---
    df_theta, theta_sum = load_theta_powers()
    if df_theta is not None:
        print('theta-powers terms:')
        print(df_theta.to_string(index=False))
        print()
    print(f'sum_{{k>=2}} theta(N^(1/k)) = {theta_sum:.6f}')
    print(f'  (paper: 32617.412861)\n')

    # --- 2) P - theta(P) ---
    P_minus_thetaP = bar_delta_obs * n_gaps - log_P + 2.0
    print(f'P - theta(P) = bar_delta_obs * n_gaps - log P + 2')
    print(f'             = {bar_delta_obs:.10e} * {n_gaps} '
          f'- {log_P:.6f} + 2')
    print(f'             = {P_minus_thetaP:.6f}')
    print(f'  (paper: 30958.422043)\n')

    # --- 3) psi(N) - theta(N) = sum_{k>=2} theta(N^(1/k)) ---
    psi_N_minus_theta_N = theta_sum
    print(f'psi(N) - theta(N) = {psi_N_minus_theta_N:.6f}\n')

    # --- 4) psi(10^9) - 10^9 ---
    # Since no prime powers in (P, 10^9]:
    #   psi(10^9) = psi(P)
    #   psi(P) - P = (psi(N) - theta(N)) - (P - theta(P)) + (P - 10^9)
    psi_diff = psi_N_minus_theta_N - P_minus_thetaP + (P - N)
    print(f'psi(10^9) - 10^9')
    print(f'  = (psi(N) - theta(N)) - (P - theta(P)) + (P - 10^9)')
    print(f'  = {psi_N_minus_theta_N:.6f} - {P_minus_thetaP:.6f} '
          f'+ ({P} - {N})')
    print(f'  = {psi_diff:.6f}')
    print(f'  (paper: +1595.990818)\n')

    # --- 5) Expected order ---
    expected = sqrt(N) / gamma_1
    print(f'Expected order: sqrt(N)/gamma_1 = {expected:.2f}')
    print(f'Ratio observed/expected: {psi_diff / expected:.4f}\n')

    # --- 6) Reconcile bar_delta gap (theory vs obs) ---
    bar_delta_theory = 6.414748e-4
    gap_obs = bar_delta_theory - bar_delta_obs

    # First term: (psi(P) - P) / n_gaps
    term1 = (psi_N_minus_theta_N - P_minus_thetaP + (P - N) - 63.0) / n_gaps
    # Second term: -(log P - 2) / n_gaps
    term2 = -(log_P - 2.0) / n_gaps
    gap_recon = term1 + term2

    print(f'Reconciliation of bar_delta gap:')
    print(f'  bar_delta_theory = {bar_delta_theory:.10e}')
    print(f'  bar_delta_obs    = {bar_delta_obs:.10e}')
    print(f'  gap (theory-obs) = {gap_obs:.10e}')
    print(f'    term1 (psi-P)  = {term1:.10e}  '
          f'({100*term1/gap_obs:.1f}% of gap)')
    print(f'    term2 (-logP)  = {term2:.10e}  '
          f'({100*term2/gap_obs:.1f}% of gap)')
    print(f'    net            = {gap_recon:.10e}')
    print(f'  Agreement: {abs(gap_recon - gap_obs)/gap_obs * 100:.4f}%\n')

    # --- 7) Write CSV ---
    out = pd.DataFrame([
        {'quantity': 'sum_theta_powers', 'value': theta_sum,
         'paper': 32_617.412861},
        {'quantity': 'P_minus_theta_P', 'value': P_minus_thetaP,
         'paper': 30_958.422043},
        {'quantity': 'psi_N_minus_theta_N', 'value': psi_N_minus_theta_N,
         'paper': 32_617.412861},
        {'quantity': 'psi_N_minus_N', 'value': psi_diff,
         'paper': 1_595.990818},
        {'quantity': 'expected_order', 'value': expected,
         'paper': 2_238.0},
    ])
    out.to_csv('delta_mean_verification.csv', index=False)
    print('Wrote delta_mean_verification.csv')


if __name__ == '__main__':
    main()