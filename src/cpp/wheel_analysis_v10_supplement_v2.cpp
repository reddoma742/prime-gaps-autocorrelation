// ================================================================
// wheel_analysis_v10_supplement_v2.cpp
//
// v10 supplement, CORRECTED:
//   [C1] task_global_correlogram now computes delta_n = g_n - log(p_n)
//        (global single subtraction, no windowed regression), then
//        passes delta to the correlogram. This matches the delta_acf_1b
//        pipeline that produced the k=1024 result reported by Kimi.
//   [C2] task_consistency: expected_sum fixed from P-2 to P-3
//        (all_gaps excludes the first gap 2→3 = 1).
//
// Unchanged from v10_supplement:
//   [G2] LO-S stay probabilities
//   [G3] rho_k / rho_1 per window
//   [G4] Consistency checks (with C2 fix)
//
// Compile:
//   g++ -O3 -std=c++17 -march=native -pthread -o wheel_v10_supp_v2 \
//       wheel_analysis_v10_supplement_v2.cpp
//
// Run:
//   ./wheel_v10_supp_v2 1000000000
//
// Outputs:
//   - correlogram_k1024_v10_v2.csv   (delta-based, correct)
//   - correlogram_k1024_raw_v10_v2.csv (raw gaps, for §4.3 extension)
//   - los_stay_v10.csv                (unchanged)
//   - los_joint_v10.csv               (unchanged)
//   - rho_over_k_v10.csv              (unchanged)
//   - results_v10_supplement_v2.txt
// ================================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <numeric>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>
#include <cstdint>
#include <random>
#include <array>
#include <climits>
#include <thread>
#include <atomic>
#include <mutex>
#include <stdexcept>

using namespace std;
using ll = long long;
using ld = long double;

// ================================================================
// Configuration
// ================================================================

const ll DEFAULT_N = 1000000000LL;
const int BLOCK_SIZE = 1 << 20;
const int K_GLOBAL_MAX = 1024;
const int K_RATIO = 20;

const vector<double> WINDOW_EDGES = {
    1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9
};

// ================================================================
// Tee streambuf
// ================================================================

class TeeBuf : public streambuf {
    streambuf* buf1; streambuf* buf2;
public:
    TeeBuf(streambuf* b1, streambuf* b2) : buf1(b1), buf2(b2) {}
protected:
    int overflow(int c) override {
        if (c == EOF) return !EOF;
        int r1 = buf1->sputc((char)c);
        int r2 = buf2->sputc((char)c);
        return (r1 == EOF || r2 == EOF) ? EOF : c;
    }
    int sync() override {
        int r1 = buf1->pubsync();
        int r2 = buf2->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }
};

// ================================================================
// Sieve
// ================================================================

struct PrimeGenerator {
    vector<int> base_primes;
    ll N;
    explicit PrimeGenerator(ll n) : N(n) {
        int limit = (int)sqrt((double)n) + 2;
        vector<bool> sieve(limit + 1, true);
        if (limit >= 0) sieve[0] = false;
        if (limit >= 1) sieve[1] = false;
        for (int i = 2; (ll)i * i <= limit; i++)
            if (sieve[i])
                for (ll j = (ll)i * i; j <= limit; j += i)
                    sieve[(size_t)j] = false;
        for (int i = 2; i <= limit; i++)
            if (sieve[i]) base_primes.push_back(i);
    }
    template <typename Func>
    void generate_all(Func callback) const {
        for (int p : base_primes)
            if ((ll)p < N) callback((ll)p);
        if (base_primes.empty()) return;
        ll start = (ll)base_primes.back() + 1;
        for (ll lo = start; lo < N; lo += BLOCK_SIZE) {
            ll hi = min(lo + (ll)BLOCK_SIZE, N);
            int sz = (int)(hi - lo);
            vector<bool> is_prime(sz, true);
            for (int p : base_primes) {
                ll s = ((lo + p - 1) / p) * p;
                if (s < (ll)p * p) s = (ll)p * p;
                if (s >= hi) continue;
                for (ll j = s - lo; j < sz; j += p)
                    is_prime[(size_t)j] = false;
            }
            for (int i = 0; i < sz; i++)
                if (is_prime[i]) callback(lo + i);
        }
    }
};

// ================================================================
// [1] Global correlogram, K = 1024 — TEMPLATED
// ================================================================

template <typename T>
struct GlobalCorrelogram {
    vector<double> rho;
    double mean = 0, var = 0;
    ll n = 0;
};

template <typename T>
GlobalCorrelogram<T> compute_global_correlogram(const vector<T>& x, int K) {
    GlobalCorrelogram<T> G;
    int n = (int)x.size();
    G.n = n;
    G.rho.assign(K + 1, 0.0);
    if (n < 2) return G;

    ld sum = 0;
    for (const auto& v : x) sum += (ld)v;
    G.mean = (double)(sum / n);

    ld var_sum = 0;
    for (const auto& v : x) {
        ld d = (ld)v - G.mean;
        var_sum += d * d;
    }
    G.var = (double)(var_sum / n);
    if (G.var < 1e-12) return G;

    vector<ld> cov(K + 1, 0.0);
    vector<ll> cnt(K + 1, 0);
    for (int i = 0; i < n; i++) {
        ld xi = (ld)x[i] - G.mean;
        int max_k = min(K, i);
        for (int k = 1; k <= max_k; k++) {
            cov[k] += xi * ((ld)x[i - k] - G.mean);
            cnt[k]++;
        }
    }
    for (int k = 1; k <= K; k++)
        if (cnt[k] > 0) G.rho[k] = (double)((cov[k] / cnt[k]) / G.var);
    return G;
}

// ----------------------------------------------------------------
// Helper: compute statistics and save CSV
// ----------------------------------------------------------------

template <typename T>
void report_and_save_correlogram(
    const GlobalCorrelogram<T>& G,
    int K,
    ll N,
    const string& csv_path,
    const string& label
) {
    cout << "  --- " << label << " ---\n";
    cout << "  n = " << G.n << "\n";

    ll n_pos = 0;
    int first_pos_k = -1;
    double max_pos = -1e-30;
    int max_pos_k = -1;
    double min_rho = +1e30, max_rho = -1e30;

    for (int k = 1; k <= K; k++) {
        if (G.rho[k] > 0) {
            n_pos++;
            if (first_pos_k < 0) first_pos_k = k;
            if (G.rho[k] > max_pos) {
                max_pos = G.rho[k];
                max_pos_k = k;
            }
        }
        min_rho = min(min_rho, G.rho[k]);
        max_rho = max(max_rho, G.rho[k]);
    }

    const int TAIL_LO = 160;
    double mean_tail = 0.0;
    ll n_tail = 0;
    for (int k = TAIL_LO; k <= K; k++) {
        mean_tail += G.rho[k];
        n_tail++;
    }
    mean_tail /= n_tail;

    double SE = 1.0 / sqrt((double)G.n);
    double z_bonf = 4.0612;
    double threshold = z_bonf * SE;

    cout << "    n_positive: " << n_pos << " / " << K << "\n";
    cout << "    first_positive_k: " << first_pos_k << "\n";
    cout << "    max_positive: +" << setprecision(8) << max_pos
         << " at k = " << max_pos_k << "\n";
    cout << "    min_rho: " << setprecision(8) << min_rho << "\n";
    cout << "    max_rho: " << setprecision(8) << max_rho << "\n";
    cout << "    SE (1/sqrt(n)) = " << setprecision(8) << SE << "\n";
    cout << "    Bonferroni threshold (z=4.0612) = " << threshold << "\n";
    cout << "    max_rho/threshold = " << setprecision(4)
         << max_rho / threshold << "\n";
    cout << "    tail_mean (k ∈ [" << TAIL_LO << ", " << K << "]) = "
         << setprecision(8) << mean_tail << "\n";
    cout << "    tail_n = " << n_tail << "\n\n";

    ofstream fout(csv_path);
    fout << "# N = " << N << "\n";
    fout << "# label = " << label << "\n";
    fout << "# n_gaps = " << G.n << "\n";
    fout << "# mean = " << setprecision(10) << G.mean << "\n";
    fout << "# var = " << setprecision(10) << G.var << "\n";
    fout << "# SE = " << setprecision(10) << SE << "\n";
    fout << "# z_bonf = " << z_bonf << "\n";
    fout << "# threshold = " << setprecision(10) << threshold << "\n";
    fout << "# n_positive = " << n_pos << "\n";
    fout << "# first_positive_k = " << first_pos_k << "\n";
    fout << "# max_rho = " << setprecision(10) << max_rho << "\n";
    fout << "# max_rho_k = " << max_pos_k << "\n";
    fout << "# tail_mean = " << setprecision(10) << mean_tail << "\n";
    fout << "k,rho_k\n";
    for (int k = 1; k <= K; k++)
        fout << k << "," << setprecision(12) << G.rho[k] << "\n";
    fout.close();
    cout << "  Saved: " << csv_path << "\n";
}

// ----------------------------------------------------------------
// Task G1 — global correlogram on DELTA (main) and RAW (comparison)
// ----------------------------------------------------------------

void task_global_correlogram_delta(
    const vector<int>& all_gaps,
    const vector<ll>& primes,
    ll N
) {
    cout << "\n[TASK G1] Global correlogram on delta_n = g_n - log(p_n)\n";
    cout << "  Direct O(n*K) loop, no FFT\n";
    cout << "  Global mean subtraction only (no windowed regression)\n";

    // ---- Build delta_n ----
    // delta_n is indexed by the gap position. The gap g_n = p_{n+1}-p_n
    // corresponds to p_n (opening prime) for n from 1 to (primes.size()-1).
    // We use the v9.1 convention: exclude the initial gap 2→3 = 1.
    vector<double> delta;
    delta.reserve(all_gaps.size());

    size_t gap_idx = 0;
    for (size_t i = 0; i + 1 < primes.size(); i++) {
        ll p = primes[i];
        int gap = (int)(primes[i + 1] - primes[i]);
        if (gap <= 1) continue;              // exclude 2→3
        delta.push_back((double)gap - log((double)p));
        gap_idx++;
    }

    cout << "  |delta| = " << delta.size() << " (should equal all_gaps)\n";
    cout << "  (all_gaps size = " << all_gaps.size() << ")\n";

    auto t0 = chrono::high_resolution_clock::now();
    auto G = compute_global_correlogram(delta, K_GLOBAL_MAX);
    double dt = chrono::duration<double>(
        chrono::high_resolution_clock::now() - t0).count();
    cout << "  Elapsed: " << setprecision(2) << dt << " s\n\n";

    report_and_save_correlogram(G, K_GLOBAL_MAX, N,
                                 "correlogram_k1024_v10_v2.csv",
                                 "delta_n = g_n - log(p_n)");
}

void task_global_correlogram_raw(
    const vector<int>& all_gaps,
    ll N
) {
    cout << "\n[TASK G1b] Global correlogram on raw gaps g_n (comparison)\n";

    auto t0 = chrono::high_resolution_clock::now();
    auto G = compute_global_correlogram(all_gaps, K_GLOBAL_MAX);
    double dt = chrono::duration<double>(
        chrono::high_resolution_clock::now() - t0).count();
    cout << "  Elapsed: " << setprecision(2) << dt << " s\n\n";

    report_and_save_correlogram(G, K_GLOBAL_MAX, N,
                                 "correlogram_k1024_raw_v10_v2.csv",
                                 "raw gaps g_n");
}

// ================================================================
// [2] LO-S stay probabilities mod 3 (unchanged)
// ================================================================

void task_los_stay(const vector<ll>& primes, ll N) {
    cout << "\n[TASK G2] LO-S stay probabilities mod 3\n";

    ll cnt_pair_1 = 0, cnt_stay_1 = 0;
    ll cnt_pair_2 = 0, cnt_stay_2 = 0;
    ll Tab[2][2] = {};

    for (size_t i = 0; i + 1 < primes.size(); i++) {
        ll p = primes[i];
        ll q = primes[i + 1];
        if (p <= 3) continue;
        int a = (int)(p % 3);
        int b = (int)(q % 3);
        if (a == 0 || b == 0) continue;
        int ai = (a == 1) ? 0 : 1;
        int bi = (b == 1) ? 0 : 1;
        Tab[ai][bi]++;
        if (a == 1) {
            cnt_pair_1++;
            if (b == 1) cnt_stay_1++;
        } else {
            cnt_pair_2++;
            if (b == 2) cnt_stay_2++;
        }
    }

    double P1 = (double)cnt_stay_1 / cnt_pair_1;
    double P2 = (double)cnt_stay_2 / cnt_pair_2;
    double SE1 = sqrt(0.25 / cnt_pair_1);
    double SE2 = sqrt(0.25 / cnt_pair_2);
    double dev1 = P1 - 0.5;
    double dev2 = P2 - 0.5;
    double z1 = dev1 / SE1;
    double z2 = dev2 / SE2;

    cout << setprecision(8);
    cout << "  Class 1 (p_n ≡ 1 mod 3):\n";
    cout << "    n_pairs  = " << cnt_pair_1 << "\n";
    cout << "    n_stay   = " << cnt_stay_1 << "\n";
    cout << "    P(stay|1) = " << P1 << "\n";
    cout << "    dev      = " << dev1 << "\n";
    cout << "    z (H0: p=1/2) = " << z1 << "\n";
    cout << "  Class 2 (p_n ≡ 2 mod 3):\n";
    cout << "    n_pairs  = " << cnt_pair_2 << "\n";
    cout << "    n_stay   = " << cnt_stay_2 << "\n";
    cout << "    P(stay|2) = " << P2 << "\n";
    cout << "    dev      = " << dev2 << "\n";
    cout << "    z (H0: p=1/2) = " << z2 << "\n";

    ofstream fout("los_stay_v10.csv");
    fout << "# N = " << N << "\n";
    fout << "# full range [2, N]\n";
    fout << "# conditioning on p_n (opening prime), p_n > 3\n";
    fout << "class,n_pairs,n_stay,P_stay,dev_from_0.5,SE,z\n";
    fout << "1," << cnt_pair_1 << "," << cnt_stay_1 << ","
         << setprecision(10) << P1 << ","
         << dev1 << "," << SE1 << "," << z1 << "\n";
    fout << "2," << cnt_pair_2 << "," << cnt_stay_2 << ","
         << setprecision(10) << P2 << ","
         << dev2 << "," << SE2 << "," << z2 << "\n";
    fout.close();
    cout << "  Saved: los_stay_v10.csv\n";

    ofstream fj("los_joint_v10.csv");
    fj << "# joint table T[a][b] = # of pairs (p_n mod 3 = a, p_{n+1} mod 3 = b)\n";
    fj << "# for a, b in {1, 2}\n";
    fj << "a,b,count\n";
    for (int a : {1, 2}) {
        for (int b : {1, 2}) {
            int ai = (a == 1) ? 0 : 1;
            int bi = (b == 1) ? 0 : 1;
            fj << a << "," << b << "," << Tab[ai][bi] << "\n";
        }
    }
    fj.close();
    cout << "  Saved: los_joint_v10.csv\n";

    ll total_from_Tab = 0;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            total_from_Tab += Tab[i][j];
    ll total_pairs = cnt_pair_1 + cnt_pair_2;
    cout << "  Consistency: total pairs = " << total_pairs
         << " (from joint table: " << total_from_Tab << ")";
    if (total_pairs == total_from_Tab)
        cout << "  PASS\n";
    else
        cout << "  FAIL\n";
}

// ================================================================
// [3] Ratio rho_k / rho_1 per window (unchanged)
// ================================================================

template <typename T>
vector<double> compute_rho_vector_fast(const vector<T>& x, int K) {
    int n = (int)x.size();
    vector<double> rho(K + 1, 0.0);
    if (n < 2) return rho;
    double mean = 0;
    for (const auto& v : x) mean += (double)v;
    mean /= n;
    double var = 0;
    for (const auto& v : x) { double d = (double)v - mean; var += d * d; }
    var /= n;
    if (var < 1e-12) return rho;
    vector<double> cov(K + 1, 0.0);
    vector<int> cnt(K + 1, 0);
    for (int i = 0; i < n; i++) {
        double xi = (double)x[i] - mean;
        int max_k = min(K, i);
        for (int k = 1; k <= max_k; k++) {
            cov[k] += xi * ((double)x[i - k] - mean);
            cnt[k]++;
        }
    }
    for (int k = 1; k <= K; k++)
        if (cnt[k] > 0) rho[k] = (cov[k] / cnt[k]) / var;
    return rho;
}

void task_rho_over_k(const vector<ll>& primes, ll N) {
    cout << "\n[TASK G3] Ratio rho_k / rho_1 per window (k=2.." << K_RATIO << ")\n";

    struct Window { double lo, hi; string name; };
    vector<Window> windows;
    for (size_t i = 0; i + 1 < WINDOW_EDGES.size(); i++) {
        if ((ll)WINDOW_EDGES[i] >= N) break;
        ostringstream oss;
        oss << "[" << (ll)WINDOW_EDGES[i] << ","
            << (ll)min(WINDOW_EDGES[i+1], (double)N) << ")";
        windows.push_back({WINDOW_EDGES[i],
                           min(WINDOW_EDGES[i+1], (double)N),
                           oss.str()});
    }

    ofstream fout("rho_over_k_v10.csv");
    fout << "# N = " << N << "\n";
    fout << "# detrended by linear regression in log p within each window\n";
    fout << "# K_RATIO = " << K_RATIO << "\n";
    fout << "window_lo,window_hi,n_gaps,mean_gap,rho_1";
    for (int k = 2; k <= K_RATIO; k++) fout << ",rho_" << k;
    for (int k = 2; k <= K_RATIO; k++) fout << ",ratio_" << k << "_over_1";
    fout << "\n";

    cout << setw(28) << "Window" << " | "
         << setw(14) << "rho_1" << " | ";
    cout << "rho_k / rho_1 for k=2..8\n";

    for (auto& w : windows) {
        vector<int> wgaps;
        vector<double> lnp;
        for (size_t i = 0; i + 1 < primes.size(); i++) {
            if ((double)primes[i] >= w.lo && (double)primes[i] < w.hi) {
                int gap = (int)(primes[i+1] - primes[i]);
                if (gap > 1) {
                    wgaps.push_back(gap);
                    lnp.push_back(log((double)primes[i]));
                }
            }
        }
        if (wgaps.size() < 1000) continue;

        double mean_lnp = 0, mean_g = 0;
        for (size_t i = 0; i < wgaps.size(); i++) {
            mean_lnp += lnp[i]; mean_g += wgaps[i];
        }
        mean_lnp /= wgaps.size();
        mean_g /= wgaps.size();
        double num = 0, den = 0;
        for (size_t i = 0; i < wgaps.size(); i++) {
            num += (lnp[i] - mean_lnp) * (wgaps[i] - mean_g);
            den += (lnp[i] - mean_lnp) * (lnp[i] - mean_lnp);
        }
        double slope = (den > 1e-15) ? num / den : 0.0;
        vector<double> dt(wgaps.size());
        for (size_t i = 0; i < wgaps.size(); i++)
            dt[i] = wgaps[i] - (mean_g + slope * (lnp[i] - mean_lnp));

        auto rho = compute_rho_vector_fast(dt, K_RATIO);

        fout << (ll)w.lo << "," << (ll)w.hi << "," << wgaps.size() << ","
             << setprecision(6) << mean_g << ","
             << setprecision(10) << rho[1];
        for (int k = 2; k <= K_RATIO; k++)
            fout << "," << setprecision(10) << rho[k];
        for (int k = 2; k <= K_RATIO; k++) {
            double ratio = (rho[1] != 0) ? (rho[k] / rho[1]) : 0.0;
            fout << "," << setprecision(10) << ratio;
        }
        fout << "\n";

        cout << "  " << setw(28) << w.name << " | "
             << setw(14) << setprecision(6) << rho[1] << " | ";
        for (int k = 2; k <= min(K_RATIO, 8); k++)
            cout << setw(7) << setprecision(3) << (rho[k] / rho[1]) << " ";
        cout << "\n";
    }
    fout.close();
    cout << "  Saved: rho_over_k_v10.csv\n";
}

// ================================================================
// [4] Consistency checks (with C2 fix: P-3)
// ================================================================

void task_consistency(
    const vector<int>& all_gaps,
    const vector<ll>& primes,
    ll N
) {
    cout << "\n[TASK G4] Consistency checks\n";

    ll prime_count = (ll)primes.size();
    ll gap_count = (ll)all_gaps.size();

    ll exp_primes = 50847534LL;
    ll exp_gaps = 50847532LL;

    cout << "  primes = " << prime_count
         << " (expected " << exp_primes << ")";
    if (N == 1000000000LL)
        cout << (prime_count == exp_primes ? "  PASS" : "  FAIL");
    cout << "\n";

    cout << "  gaps = " << gap_count
         << " (expected " << exp_gaps << ")";
    if (N == 1000000000LL)
        cout << (gap_count == exp_gaps ? "  PASS" : "  FAIL");
    cout << "\n";

    // Sum of all_gaps = P - 3, because:
    //   sum of ALL gaps from 2 to P is P - 2
    //   all_gaps excludes the first gap (2→3 = 1)
    //   So sum(all_gaps) = (P - 2) - 1 = P - 3.
    ll P = primes.empty() ? 0 : primes.back();
    ll sum_gaps = 0;
    for (int g : all_gaps) sum_gaps += g;
    ll expected_sum = P - 3;   // ★ C2 fix (was P - 2)
    cout << "  sum(gaps) = " << sum_gaps
         << " (expected P-3 = " << expected_sum << ")";
    if (sum_gaps == expected_sum)
        cout << "  PASS  [note: all_gaps excludes the first gap 2→3=1]\n";
    else
        cout << "  FAIL\n";
}

// ================================================================
// MAIN
// ================================================================

int main(int argc, char* argv[]) {
    ll N = DEFAULT_N;
    if (argc > 1) {
        try { N = stoll(argv[1]); }
        catch (...) { cerr << "Invalid N. Using default.\n"; }
    }
    if (N < 1000) { cerr << "N too small.\n"; return 1; }

    ofstream log_file("results_v10_supplement_v2.txt");
    streambuf* original_cout = cout.rdbuf();
    TeeBuf tee(original_cout, log_file.rdbuf());
    cout.rdbuf(&tee);

    auto t_start = chrono::high_resolution_clock::now();
    cout << fixed << setprecision(6);

    cout << "==========================================================================================\n";
    cout << "  PRIME WHEELS — v10 SUPPLEMENT v2 (CORRECTED)\n";
    cout << "  N = " << N << "\n";
    cout << "  K_GLOBAL_MAX = " << K_GLOBAL_MAX << "\n";
    cout << "  [C1] task_global_correlogram uses delta_n = g_n - log(p_n)\n";
    cout << "  [C2] task_consistency uses P-3\n";
    cout << "  Purpose: [1] K=1024 correlogram on delta, [2] LO-S stay,\n";
    cout << "           [3] rho_k/rho_1, [4] consistency\n";
    cout << "  Does NOT re-run bootstrap or permutation (v9.1 provides those).\n";
    cout << "==========================================================================================\n";

    // ---------------- Sieve ----------------
    cout << "\n[Phase 0] Sieving...\n";
    vector<int> all_gaps;
    vector<ll> primes;
    all_gaps.reserve((size_t)(N / 15));
    primes.reserve((size_t)(N / 15));

    PrimeGenerator gen(N);
    ll prev = -1, prime_count = 0, last_rep = 0;
    ld sum_g = 0;

    gen.generate_all([&](ll p) {
        prime_count++;
        primes.push_back(p);
        if (prev > 0) {
            int gap = (int)(p - prev);
            if (gap > 1) {
                all_gaps.push_back(gap);
                sum_g += gap;
            }
        }
        prev = p;
        if (p - last_rep > N / 20) {
            cerr << "\r  Progress: " << setprecision(1)
                 << (100.0 * (double)p / (double)N) << "%"
                 << "  primes: " << prime_count << "     " << flush;
            last_rep = p;
        }
    });
    cerr << "\r  Progress: 100%  primes: " << prime_count << "     \n";

    cout << "  Primes: " << prime_count << "\n";
    cout << "  Gaps:   " << all_gaps.size() << "\n";
    cout << "  Time:   " << setprecision(2)
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_start).count()
         << " s\n";

    // ---------------- Task G1: delta-based correlogram (MAIN) ----------------
    task_global_correlogram_delta(all_gaps, primes, N);

    // ---------------- Task G1b: raw-gap correlogram (comparison) ----------------
    task_global_correlogram_raw(all_gaps, N);

    // ---------------- Task G2 ----------------
    task_los_stay(primes, N);

    // ---------------- Task G3 ----------------
    task_rho_over_k(primes, N);

    // ---------------- Task G4 ----------------
    task_consistency(all_gaps, primes, N);

    // ---------------- Summary ----------------
    cout << "\n==========================================================================================\n";
    cout << "[FINAL] Summary\n";
    cout << "  N = " << N << "\n";
    cout << "  primes = " << prime_count << "\n";
    cout << "  gaps = " << all_gaps.size() << "\n";
    cout << "  Total time: " << setprecision(2)
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_start).count()
         << " s\n";
    cout << "\n  Files saved:\n";
    cout << "    - correlogram_k1024_v10_v2.csv     (delta-based, main)\n";
    cout << "    - correlogram_k1024_raw_v10_v2.csv (raw gaps, comparison)\n";
    cout << "    - los_stay_v10.csv\n";
    cout << "    - los_joint_v10.csv\n";
    cout << "    - rho_over_k_v10.csv\n";
    cout << "    - results_v10_supplement_v2.txt\n";
    cout << "==========================================================================================\n";

    cout.flush();
    log_file.flush();
    cout.rdbuf(original_cout);
    return 0;
}
