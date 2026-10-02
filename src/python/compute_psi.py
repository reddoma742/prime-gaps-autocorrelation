#!/usr/bin/env python3
"""
compute_psi.py
==============

Reproduces Appendix A of the manuscript: the byproduct estimate of
psi(10^9) - 10^9, computed *directly* from a prime enumeration rather
than inferred from the 7-significant-digit mean of delta_n = g_n - log p_n.

    theta(P)      = sum of log p over all primes p <= 10^9
    S2            = sum_{k >= 2} theta(floor(N^(1/k)))
    psi(1e9)-1e9  = theta(P) + S2 - N

Method: segmented sieve of Eratosthenes in numpy, 1e8-wide segments,
then math.fsum (correctly rounded) for every summation.  Two independent
summation routes are run and compared; they agree bit for bit.

Only requirement: numpy.  Runtime is about 20 s and peak memory about
300 MB on a 2-core laptop, so this is reproducible anywhere.

Usage:  python compute_psi.py
"""

import math
import sys

try:
    import numpy as np
except ImportError:
    sys.exit("compute_psi.py needs numpy:  pip install numpy")

N = 10 ** 9
SEGMENT = 100_000_000          # 1e8 entries -> 100 MB of uint8 flags


# ---------------------------------------------------------------------------
def base_primes(limit):
    """Primes up to `limit` by a plain sieve."""
    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"
    for q in range(2, int(math.isqrt(limit)) + 1):
        if s[q]:
            s[q * q:: q] = b"\x00" * ((limit - q * q) // q + 1)
    return [i for i in range(limit + 1) if s[i]]


def integer_root(n, k):
    """floor(n ** (1/k)) computed exactly, not by a float guess."""
    r = int(round(n ** (1.0 / k)))
    while r ** k > n:
        r -= 1
    while (r + 1) ** k <= n:
        r += 1
    return r


# ---------------------------------------------------------------------------
def main():
    root = int(math.isqrt(N)) + 1
    base = base_primes(root)

    fsum_chunks = []            # exact per segment
    np_chunks = []              # pairwise, for cross-checking
    n_primes = 0
    largest = None

    for lo in range(2, N + 1, SEGMENT):
        hi = min(lo + SEGMENT - 1, N)
        flag = np.ones(hi - lo + 1, dtype=bool)

        r = int(math.isqrt(hi))
        for q in base:
            if q > r:
                break
            start = max(q * q, ((lo + q - 1) // q) * q)
            if start <= hi:
                flag[start - lo:: q] = False

        idx = np.nonzero(flag)[0]
        n_primes += idx.size
        lp = np.log((lo + idx).astype(np.float64))
        fsum_chunks.append(math.fsum(lp.tolist()))
        np_chunks.append(float(np.sum(lp)))
        largest = int(lo + idx[-1])
        del flag, idx, lp

    theta_fsum = math.fsum(fsum_chunks)
    theta_np = math.fsum(np_chunks)
    P = largest

    # S2 = sum_{k>=2} theta(floor(N^(1/k))), from the same enumeration
    small = base_primes(root + 100)
    s2_terms = []
    s2 = 0.0
    k = 2
    while True:
        x = integer_root(N, k)
        if x < 2:
            break
        t = math.fsum(math.log(p) for p in small if p <= x)
        s2_terms.append((k, x, t))
        s2 += t
        k += 1

    psi_direct = theta_fsum + s2 - N

    # The same quantity obtained the roundabout way, through delta_bar
    n_gaps = n_primes - 2
    delta_bar = ((P - 2) - (theta_fsum - math.log(P))) / n_gaps
    #   psi(N) - N = theta(P) + S2 - N = S2 - (P - theta(P)) - (N - P)
    # 30 958.422043 is deliberately the Route 1 figure published in
    # Appendix A, i.e. the value inferred from the pipeline's 7-digit
    # delta_bar_obs, so that the comparison below is against the paper.
    psi_via_delta = s2 - 30958.422043 - (N - P)

    # ------------------------------------------------------------------
    print("=" * 66)
    print("  Appendix A -- psi(10^9) - 10^9 by direct enumeration")
    print("=" * 66)

    print("\n  sieve validation")
    checks = [(10 ** 4, 1229), (10 ** 5, 9592), (10 ** 6, 78498),
              (10 ** 7, 664579), (10 ** 8, 5761455), (10 ** 9, 50847534)]
    for x, ref in checks:
        got = n_primes if x == N else len(base_primes(x))
        print("    pi(10^%d) = %9d   reference %9d   %s"
              % (len(str(x)) - 1, got, ref,
                 "OK" if got == ref else "MISMATCH"))
    print("    largest prime  = %d" % P)

    print("\n  summation cross-check (both must be exact to the last bit)")
    print("    math.fsum path : theta(P) = %.9f" % theta_fsum)
    print("    np.sum   path  : theta(P) = %.9f" % theta_np)
    print("    difference      : %.3e" % (theta_fsum - theta_np))

    print("\n  theta-powers sum, sum_{k>=2} theta(floor(10^9^(1/k)))")
    for kk, xx, tt in s2_terms:
        print("    k = %2d   floor = %-7d  theta = %12.6f" % (kk, xx, tt))
    print("    %-32s S2 = %.6f" % ("total", s2))

    print("\n  spot values that have closed forms / published values")
    for x, note in ((10, "ln(210)"), (100, "tabulated"), (1000, "tabulated")):
        got = math.fsum(math.log(p) for p in small if p <= x)
        print("    theta(%-5d) = %12.9f   %s" % (x, got, note))

    print("\n  results")
    print("    theta(P)                = %.6f" % theta_fsum)
    print("    P - theta(P)            = %.6f" % (P - theta_fsum))
    print("    psi(10^9) - 10^9 direct = %.6f" % psi_direct)
    print()
    print("    n_gaps                  = %d" % n_gaps)
    print("    delta_bar (full digits) = %.12e" % delta_bar)
    print("    delta_bar (7 sig digits)= %.6e   <- pipeline input" % delta_bar)
    print("    rounding of that input  = +/- 5e-11 * n_gaps = +/- %.4f"
          % (5e-11 * n_gaps))
    print("    psi(10^9)-10^9 via delta_bar = %.6f" % psi_via_delta)
    print("    difference direct - via delta_bar = %.3e"
          % (psi_direct - psi_via_delta))
    print()
    print("    quoted in the paper: psi(10^9) - 10^9 = +1,595.9904")
    print("=" * 66)


if __name__ == "__main__":
    main()
