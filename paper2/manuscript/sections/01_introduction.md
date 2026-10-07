# 1. Introduction

Let \(p_n\) be the \(n\)-th prime and let
\[
g_n = p_{n+1}-p_n
\]
be the \(n\)-th consecutive prime gap. The distribution and size of prime gaps have been studied extensively since Cramér (1936), with major contributions by Gallagher (1976) and Granville (1995), and more recent work on arithmetic structure and correlations by Holt (arXiv:2405.03540), the Lemke Oliver–Soundararajan work on consecutive-prime biases (LOS, arXiv:1603.03720), the Montgomery–Soundararajan framework (MS, arXiv:math/0409258), Ash–Beltis–Gross–Sinnott (2011), and Banks–Ford–Tao (arXiv:1908.08613). Much less is known, however, about the serial dependence of consecutive prime gaps after local scale effects have been removed.

In Berramdane (2026a, DOI 10.5281/zenodo.23162269), the first paper in this series, we studied locally detrended prime-gap windows and reported persistent negative lag-one dependence up to \(N=10^9\). The present paper extends that investigation to the range
\[
[10^{10},\,3\times 10^{10})
\]
and introduces a sieve-phase conditioning step. The central observable is the lag-one covariance of consecutive prime gaps after local detrending and after conditioning on the opening-prime residue modulo a primorial \(M\).

The data set analyzed here contains
\[
n = 844{,}953{,}415
\]
gaps in the interval \([10^{10},3\times 10^{10})\). For the detrended gap sequence we find
\[
\operatorname{Var}(e)=438.670,\qquad \operatorname{sd}(e)=20.944,
\]
with
\[
d = L-\operatorname{sd} = 2.631,\qquad
|\operatorname{Cov}|/(L/2)=0.933,
\]
and
\[
\rho_1=-0.025078.
\]
These values reproduce the frozen predictions made before the \(10^{10}\) run.

The main technical device is an exact covariance decomposition. For each primorial modulus
\[
M\in\{6,30,210,2310,30030,510510,9699690\},
\]
we decompose the lag-one covariance into four terms
\[
A+X+Y+D=\operatorname{Cov},
\]
where \(A\) is the phase component, \(X\) and \(Y\) are cross terms, and \(D\) is the residual term. The decomposition closes exactly to within \(10^{-6}\) for every \(M\) in the table. The residual correlation \(\rho_{\text{resid}}\) remains negative for all seven moduli, ranging from \(-0.024703\) at \(M=6\) to \(-0.010306\) at \(M=9699690\). Thus phase conditioning reduces the residual term D from −10.83 at M=6 to −4.35 at M=9699690 and redistributes the covariance among A, X, Y, D, while the total remains fixed at −11.001171. The residual correlation remains negative at all tested scales.

This paper makes five contributions. First, it reproduces the frozen \(10^{10}\) predictions and confirms the stability of the negative lag-one dependence at the larger scale. Second, it gives an exact sieve-phase decomposition of the lag-one covariance. Third, it reports the decomposition across seven primorial scales. Fourth, it studies the variance structure of the gap sequence through between-phase and within-phase components. Fifth, it states the limitations of the present analysis and distinguishes the observable studied here from related recurrence-statistics approaches.

We do not claim that phase conditioning explains all dependence, nor that the residual converges to zero, nor that the present results prove the Hardy–Littlewood or Montgomery–Soundararajan laws. The aim is narrower: to measure, decompose, and stress-test the lag-one covariance of consecutive prime gaps in a controlled windowed setting, and to report what remains after the dominant sieve-phase component has been removed.

The paper is organized as follows. Section 2 defines the data and windowing procedure. Section 3 presents the sieve-phase decomposition. Section 4 gives the covariance decomposition. Section 5 reports scale validation. Section 6 analyzes the variance structure. Section 7 discusses theory and records a finite inclusion–exclusion diagnostic as future work. Section 8 states limitations. Section 9 records the distinction from Murray (SSRN 6947578, SSRN 7426882).
