# IE Diagnostic — Future Work Protocol (Frozen, Not Executed)

Status: FROZEN, NOT EXECUTED.

This protocol is written in advance so that, if the IE diagnostic is ever
run, it is pre-registered. Running it is optional. Its outcome does not
modify any empirical claim of Paper 2.

## Purpose

Estimate the order of magnitude of d(L) using a finite truncated
Hardy–Littlewood inclusion–exclusion calculation for consecutive gaps.

## Not a hypothesis test

No p-value, confidence level, or theorem-level prediction is claimed.
The output is a descriptive magnitude, to be compared against
pre-registered descriptive bands.

## Input

The same logarithmic windows used in Paper 2:
  - [10^10, 3 x 10^10)

## Pair density

Specify the exact Hardy–Littlewood pair-density formula and normalization.
Choice to be recorded here before execution. Default proposal:
  Use the standard HL pair singular series with the full Euler product
  truncated at primes p <= 10^6.

## Intermediate-prime exclusion

List the allowed intermediate offsets, truncation depth, and treatment of
overlapping tuples. Default proposal:
  - Allowed intermediate offsets: even offsets in {2, 4, ..., 2*k_max}
  - k_max = 20
  - Overlapping tuples: resolved by standard Möbius-style inclusion–exclusion

## Tail

Specify h_max, tail treatment, and whether renormalization is applied.
Default proposal:
  - h_max = 500
  - Tail: absorbed by renormalization to sum P_IE(h) = 1
  - Renormalization: applied once, at the end

## Output

P_IE(h), E_IE[h], E_IE[h^2], Var_IE(g), sd_IE(g), d_IE(L).

## Pre-registered descriptive bands

  2.3 <= d_IE <= 2.7       : agreement in magnitude
  2.0 <= d_IE < 2.3        : partial agreement
  2.7 <  d_IE <= 3.0       : partial agreement
  d_IE < 2.0 or d_IE > 3.0 : does not account for the observed level

These are descriptive bands, not statistical-significance thresholds.

## Interpretation

The calculation is a diagnostic heuristic only. It neither defines the
empirical estimator nor proves an asymptotic law. It is reported as a
finite truncated inclusion–exclusion diagnostic, not as a test.

## Execution note

Executing this protocol requires a separate task with its own ID.
The current task paper2_build_clean_layout_v1 does not execute it.
