# 8. Limitations

The present study is a finite-range empirical analysis. The results reported in Sections 3–6 are conditional on the specific window, the specific detrending, and the specific ladder of moduli chosen. The limitations below are stated so that no reader mistakes the reported decompositions for asymptotic statements.

## 8.1 Finite window

The analysis is confined to the single window \([10^{10},\,3\times 10^{10})\), containing \(844{,}953{,}415\) consecutive prime gaps. No claim is made about the behaviour of the lag-one covariance outside this range. The agreement with the pre-registered predictions (Section 5) is an out-of-sample consistency check at a fixed scale, not evidence of an asymptotic law.

## 8.2 Finite ladder of moduli

Seven primorial moduli are used:
\[
M \in \{6,\,30,\,210,\,2310,\,30030,\,510510,\,9699690\}.
\]
The monotone growth of `var_between` and the monotone decline in magnitude of \(\rho_{\text{resid}}\) across this ladder are observed facts about these seven points. No extrapolation to \(M\to\infty\) is supported. In particular, the data do not establish that \(\rho_{\text{resid}}(M)\to 0\), nor that it converges to any nonzero limit.

A note on finite-sample estimation. At \(M=9{,}699{,}690\) the number of phase classes is \(\varphi(M) = 1{,}658{,}880\), with an average of 509 gaps per class and a smallest class of 433. The phase classes form an exact partition, so the decomposition itself is not affected by empty classes. Finite class sizes nevertheless affect the stability of the estimated conditional means and hence of \(D\), especially at the largest modulus. We therefore interpret the largest-modulus result as a high-resolution diagnostic rather than as an independent validation. The observed decline in \(\rho_{\text{resid}}(M)\) is not asserted to be free of finite-class estimation error.

## 8.3 No causal claim

The decomposition
\[
\operatorname{Cov}(e_n,\,e_{n+1}) = A + X + Y + D
\]
is an algebraic identity. All four terms are computed from the same data. The decomposition partitions the covariance; it does not identify a causal mechanism. In particular, the phase component \(A\) should not be read as "the cause" of the observed covariance, and the residual term \(D\) should not be read as "the true dependence" stripped of phase effects.

## 8.4 Detrending is a modelling choice

The local detrending \(g_n \mapsto e_n = g_n - (\alpha + \beta\log p_n)\) is a specific modelling choice, frozen in the protocol and inherited from Berramdane (2026a). A different detrending would produce different residuals and hence a different decomposition. The consistency of the present results with the frozen \(10^9\) predictions is a check on this choice, not a proof that it is unique.

## 8.5 Phase conditioning is on the opening prime only

Conditioning is performed on \(r_n = p_n \bmod M\), the residue of the opening prime. The closing prime \(p_{n+1}\) also has a residue modulo \(M\), and joint conditioning on both residues would produce a finer decomposition. The present analysis does not pursue this refinement; the observable of interest is the lag-one covariance of numerical gap sizes, conditioned on the opening residue.

## 8.6 Bootstrap scope

The 95% intervals for \(\rho_{\text{resid}}\) are obtained by a moving-block bootstrap with 40 contiguous blocks and 1000 replicates. The block bootstrap is appropriate under the assumption that the residual sequence is approximately stationary within the window and that the block length is large compared to the dependence range. The intervals are conditional on these choices. They are not i.i.d. bootstrap intervals, and they do not constitute a test of stationarity.

## 8.7 No relation to HL/MS as a theorem

The present results are compatible with a sieve-phase interpretation of the lag-one covariance, but they do not prove the Hardy–Littlewood conjecture, the Montgomery–Soundararajan law, or any asymptotic statement about prime gaps. The Montgomery–Soundararajan framework provides a conditional variance-consistency benchmark for a cumulative correlation statistic; applying it to a fixed number of consecutive gaps is an extrapolation beyond its stated asymptotic regime, and is not made here.

## 8.8 No comment on Murray

This paper does not evaluate, extend, or contradict the results of Murray (SSRN 6947578, SSRN 7426882). The distinction between the two observables is recorded in Section 9. No comparison of numerical values is made, because the observables are different.

## 8.9 The IE diagnostic is optional future work

The inclusion–exclusion (IE) diagnostic discussed in Section 7.6 is a possible future-work heuristic. It is not used to define the observable, and its outcome would not modify the definitions of \(A\), \(X\), \(Y\), or \(D\). A positive or negative IE result would change the theoretical interpretation, not the measured decomposition.

## 8.10 Reproducibility

All numerical values in this paper are computed from the frozen \(10^{10}\) run. The full source, the full output log, the CSV, and the SHA256 checksums are archived on GitHub (commit `bb3548c`, `paper2/SHA256SUMS`). No value is reported without a corresponding run. The small-scale validation at \(10^6\) and \(10^7\) passes before the large run; the two-pass streaming implementation uses 64-bit prime coordinates and does not reconstruct primes from floating-point logarithms.
