# 5. Scale validation

We repeated the frozen analysis on the disjoint opening-prime window
\[
[10^{10},\,3\times 10^{10}),
\]
containing \(844{,}953{,}415\) gaps.

Before the run, the following expected ranges were specified in the frozen protocol:

- \(\operatorname{Var}(e)\in[438,\,445]\),
- \(d = L-\operatorname{sd}(e)\in[2.5,\,2.7]\),
- \(\operatorname{Cov}(e_n,e_{n+1})\in[-11.8,\,-11.0]\),
- \(\rho_1\in[-0.0269,\,-0.0248]\).
- \(\operatorname{sd}(e)\in[20.9,\,21.1]\),
- \(|\operatorname{Cov}|/(L/2)\in[0.93,\,1.00]\).

The observed values were
\[
\operatorname{Var}(e) = 438.670412,\quad
\operatorname{sd}(e) = 20.944,
\]
\[
d = L-\operatorname{sd}(e) = 2.631,\quad
\operatorname{Cov}(e_n,e_{n+1}) = -11.001171,
\]
\[
\rho_1 = -0.025078,\quad
|\operatorname{Cov}|/(L/2) = 0.933285.
\]
All six values fall inside their pre-specified ranges. Two of them lie close to a range boundary: Cov = −11.001171 lies 0.001171 above the lower edge of [−11.8, −11.0], and |Cov|/(L/2) = 0.933285 lies 0.003285 above the lower edge of [0.93, 1.00]. This agreement is a finite-range out-of-sample consistency check, not an asymptotic validation. The pre-registered ranges were specified for the detrended quantities in the frozen protocol file `extensions/ie_test/future_work/protocol.md`, committed at source/data commit `bb3548c` before the 10^10 run.

The window is defined by the opening prime \(p_n\); the closing prime of the final gap may exceed the upper endpoint. The two-pass streaming implementation used 64-bit prime coordinates and did not reconstruct primes from floating-point logarithms. This is the correction that resolved the earlier overflow at \(p > 2^{32}\).
