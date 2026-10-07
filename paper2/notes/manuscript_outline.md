# Manuscript outline

## Title
Windowed Phase Conditioning and Lag-One Covariance of Consecutive Prime Gaps

## Abstract
We measure how conditioning on opening-prime residue classes modulo primorial moduli changes the lag-one covariance of consecutive prime gaps. Using logarithmic windows, local detrending, block bootstrap intervals, and an independent validation window extending to 3x10^10, we decompose the covariance into phase, cross, and within-phase terms. The phase-conditioned residual remains negative at every tested scale. The results are presented as finite-range empirical evidence; no asymptotic law or theorem-level identification with Hardy–Littlewood or Montgomery–Soundararajan is claimed.

## 1. Introduction
- Prime-gap lag-one dependence.
- Why pooled nonstationarity is dangerous.
- Distinction from residue recurrence statistics.
- Relation to primorial projection studies (Murray).

## 2. Data and definitions
- N=10^9 main dataset and [10^10,3x10^10) validation.
- Opening prime p_n and gap g_n.
- Windows and log-p detrending.
- Phase r_n=p_n mod M.

## 3. Estimators
- Var, Cov, rho_1.
- Conditional means mu_r and residual u_n.
- A+B+C+D decomposition.
- Class-size and overfitting correction.
- Bootstrap and jackknife protocol.

## 4. Primorial phase ladder
- M=6,30,210,2310,30030,510510,9699690 where class sizes permit.
- Report delta_rho2 and residual rho.
- Emphasize that the ladder is descriptive, not an asymptotic limit.

## 5. Covariance and variance decomposition
- Between/within variance.
- A,B,C,D across windows.
- Sign and scale of each component.

## 6. Scale validation at 10^10
- Frozen predictions and observed results.
- n=844,953,415 gaps.
- Var=438.670, sd=20.944, d=2.631, Cov=-11.001, rho_1=-0.025078.
- Exact ABCD closure.

## 7. Relation to theory
- HL singular series as conditional heuristic.
- LO-S as residue-bias context.
- Montgomery–Soundararajan as a conditional variance consistency check.
- Explicit statement that the extrapolation is outside the theorem's stated range.

## 8. Limitations
- Six windows and finite moduli.
- Dependence between phase classes and overlapping pairs.
- No claim of causal mechanism.
- No proof that residual tends to zero.

## Appendices
- A: estimator details.
- B: bootstrap/jackknife details.
- C: full CSV tables and checksums.
