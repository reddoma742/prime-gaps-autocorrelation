// ================================================================
// theta_powers.cpp
//
// Computes the exact theoretical value of
//     δ_mean_theory = (1 / n_gaps) · Σ_{k ≥ 2} θ(N^{1/k})
// where θ(x) = Σ_{p ≤ x} log p.
//
// Uses the telescoping identity:
//     Σ_i (g_i − log p_i) = (P − 2) − θ(P) + log P
// with P the largest prime ≤ N.
//
// For N = 10^9: P = 999,999,937.
//
// Compile:
//   g++ -O3 -std=c++17 -o theta_powers theta_powers.cpp
//
// Run:
//   ./theta_powers 1000000000
//
// Output:
//   - theta_powers_terms.csv: k, N^{1/k}, θ(N^{1/k}) for k ≥ 2
//   - stdout: sum, δ_mean_theory, comparison with observed
// ================================================================

#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using ld = long double;

int main(int argc, char** argv) {
    ll N = (argc > 1) ? stoll(argv[1]) : 1000000000LL;
    auto t0 = chrono::high_resolution_clock::now();

    // ---- 1. Sieve of Eratosthenes up to sqrt(N) ----
    ll limit = (ll)sqrt((ld)N) + 2;
    vector<bool> sieve(limit + 1, true);
    if (limit >= 0) sieve[0] = false;
    if (limit >= 1) sieve[1] = false;
    for (ll i = 2; i * i <= limit; i++)
        if (sieve[i])
            for (ll j = i * i; j <= limit; j += i)
                sieve[(size_t)j] = false;

    vector<int> primes;
    for (ll i = 2; i <= limit; i++)
        if (sieve[i]) primes.push_back((int)i);

    cerr << "Sieved " << primes.size()
         << " primes up to sqrt(N) = " << limit << "\n";

    // ---- 2. θ(x) via binary search on the prime vector ----
    auto theta = [&](ll x) -> ld {
        if (x < 2) return 0.0L;
        int hi = (int)(upper_bound(primes.begin(), primes.end(), (int)x)
                       - primes.begin());
        ld s = 0.0L;
        for (int i = 0; i < hi; i++) s += logl((ld)primes[i]);
        return s;
    };

    // ---- 3. Integer k-th root (largest r with r^k ≤ N) ----
    auto kth_root = [&](int k) -> ll {
        if (k == 1) return N;
        ll r = (ll)pow((ld)N, 1.0L / (ld)k);
        // correct downward
        while (r > 1) {
            ll rk = 1;
            bool over = false;
            for (int i = 0; i < k; i++) {
                if (rk > N / r) { over = true; break; }
                rk *= r;
            }
            if (!over && rk <= N) break;
            r--;
        }
        // correct upward
        while (true) {
            ll rk = 1;
            bool over = false;
            for (int i = 0; i < k; i++) {
                if (rk > N / (r + 1)) { over = true; break; }
                rk *= (r + 1);
            }
            if (over || rk > N) break;
            r++;
        }
        return r;
    };

    // ---- 4. Sum over k ≥ 2 of θ(N^{1/k}) ----
    ld sum_theta = 0.0L;
    int n_terms = 0;
    int kmax = (int)floorl(log2l((ld)N)) + 2; // safety margin

    ofstream fcsv("theta_powers_terms.csv");
    fcsv << "k,N_root,theta\n";
    fcsv << fixed << setprecision(10);

    cout << fixed << setprecision(6);
    cout << "  Individual terms θ(N^(1/k)):\n";
    cout << "  " << setw(4) << "k" << "  "
         << setw(14) << "N^(1/k)" << "  "
         << setw(18) << "theta" << "\n";
    cout << "  " << string(40, '-') << "\n";

    for (int k = 2; k <= kmax; k++) {
        ll root = kth_root(k);
        if (root < 2) break;
        ld th = theta(root);
        sum_theta += th;
        n_terms++;

        cout << "  " << setw(4) << k << "  "
             << setw(14) << root << "  "
             << setw(18) << (double)th << "\n";
        fcsv << k << "," << root << "," << (double)th << "\n";
    }
    fcsv.close();

    // ---- 5. Known constants for N = 10^9 ----
    ll n_primes_N = (N == 1000000000LL) ? 50847534LL : 0;
    ll n_gaps_N   = (N == 1000000000LL) ? 50847532LL : 0;
    ll P          = (N == 1000000000LL) ? 999999937LL : 0;
    ld delta_obs  = 0.0006092163L;

    cout << "\n========================================\n";
    cout << "  theta-powers calculation\n";
    cout << "========================================\n";
    cout << "  N             = " << N << "\n";
    cout << "  n_terms       = " << n_terms << "\n";
    cout << "  sum_theta     = " << (double)sum_theta << "\n";

    if (n_gaps_N > 0) {
        ld delta_theory = sum_theta / (ld)n_gaps_N;
        cout << "\n  n_primes(N)   = " << n_primes_N << "\n";
        cout << "  n_gaps        = " << n_gaps_N << "\n";
        cout << "  P             = " << P << "\n";
        cout << "\n  delta_mean (theory)   = " << (double)delta_theory << "\n";
        cout << "  delta_mean (observed) = " << (double)delta_obs << "\n";
        cout << "  ratio                 = "
             << (double)(delta_theory / delta_obs) << "\n";
        cout << "  difference            = "
             << (double)(delta_theory - delta_obs) << "\n";

        // Extract ψ(N) − N (assuming identity)
        //   sum_theta = ψ(N) − θ(N)
        //   P − θ(P) = n_gaps · δ_obs − log P + 2
        ld log_P = logl((ld)P);
        ld P_minus_thetaP = (ld)n_gaps_N * delta_obs - log_P + 2.0L;
        ld psiP_minus_P = sum_theta - P_minus_thetaP;
        ld psiN_minus_N = psiP_minus_P + ((ld)P - (ld)N);

        cout << "\n  --- Extraction of ψ(N) − N ---\n";
        cout << "  log P                 = " << (double)log_P << "\n";
        cout << "  P − θ(P)              = " << (double)P_minus_thetaP << "\n";
        cout << "  ψ(P) − P              = " << (double)psiP_minus_P << "\n";
        cout << "  ψ(N) − N              = " << (double)psiN_minus_N << "\n";
        cout << "  (expected ~ +1596 based on sqrt(N)/γ₁ ≈ 2238)\n";
    }

    auto t1 = chrono::high_resolution_clock::now();
    double sec = chrono::duration<double>(t1 - t0).count();
    cout << "\n  Elapsed: " << sec << " s\n";
    cout << "  Saved: theta_powers_terms.csv\n";

    return 0;
}
