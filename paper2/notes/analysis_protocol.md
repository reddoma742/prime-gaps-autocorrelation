# Frozen analysis protocol

## Windows
Main: [10^6,3x10^6), [3x10^6,10^7), [10^7,3x10^7), [3x10^7,10^8), [10^8,3x10^8), [3x10^8,10^9).
Validation: [10^10,3x10^10), opening prime in the stated interval.

## Moduli
M=2310 and 30030 for all main windows. M=510510 only where n_c>=50. M=9699690 is diagnostic only unless minimum class size is documented.

## Detrending
Within each window fit g_n = alpha + beta log(p_n), and use residuals e_n.

## Phase residual
mu_r is the sample mean of e_n among p_n mod M=r. Set u_n=e_n-mu_r.

## Report
Report raw and corrected rho_resid, n_c, minimum class count, A,B,C,D, ABCD-Cov, and bootstrap CI. Do not apply the (1-1/n_c) correction to A,B,C,D.

## Bootstrap
Use independent RNG streams for each window and modulus. Record seed, block count, block length, number of replicates, and valid replicate count.
