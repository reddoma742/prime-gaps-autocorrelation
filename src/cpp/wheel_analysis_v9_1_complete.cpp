// ================================================================
// wheel_analysis_v9_1_complete.cpp
//
// Single-file, end-to-end analysis of prime gaps up to N = 10^9.
//
// This is v9.1 with corrections from kimi's review:
//   [C1] gap_pairs built directly in sieve loop (no rebuild)
//   [C2] SingularSeries cutoff renamed (finite product, not analytic tail)
//   [C3] S() throws on pattern too large (no silent 0.0)
//   [C4] S() cached via unordered_map
//   [C5] Consistency checks (gaps, pairs, marginals, known values)
//   [C6] CSV metadata renamed correctly
//
// Compile:
//   g++ -O3 -std=c++17 -march=native -pthread -o wheel_v9_1 wheel_analysis_v9_1_complete.cpp
//
// Run (small test first):
//   ./wheel_v9_1 1000000
//   ./wheel_v9_1 10000000
//   ./wheel_v9_1                       # default N = 10^9
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
#include <unordered_map>
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
// Global configuration
// ================================================================

const ll DEFAULT_N = 1000000000LL;
const int BLOCK_SIZE = 1 << 20;
const int K_SHORT = 10;
const int K_PERM = 10;
const int K_EXTENDED = 64;
const int N_PERM = 2000;
const int N_BOOT = 5000;
const int N_THREADS = min(16u, max(1u, thread::hardware_concurrency()));

// Task 9: inclusion-exclusion model
const int IE_R_MAX = 4;
const int IE_PRIME_LIMIT = 1000000;   // finite product cutoff
const int K_MAX_IE = 20;

// Seeds
const uint64_t SEED_PERM = 42;
const uint64_t SEED_BOOT = 7919;
const uint64_t SEED_WIN_BOOT = 1234;

// 17 target pairs
const vector<pair<int,int>> PAIRS_17 = {
    {2,10}, {10,2}, {2,28}, {28,2},
    {14,16}, {16,14}, {26,4}, {4,26},
    {6,6}, {6,12}, {12,6}, {6,4}, {4,6},
    {2,4}, {4,2}, {2,6}, {6,2}
};

// ================================================================
// Utility
// ================================================================

void print_sep(int width = 90, char c = '=') {
    for (int i = 0; i < width; i++) cout << c;
    cout << '\n';
}

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
// Part 1: Segmented sieve of Eratosthenes
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
// Part 2: Rho computation (single-pass for all lags)
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

template <typename T>
double compute_rho_single_fast(const vector<T>& x, int k) {
    if (k < 1 || (int)x.size() <= k + 1) return 0.0;
    int n = (int)x.size();
    double mean = 0;
    for (const auto& v : x) mean += (double)v;
    mean /= n;
    double var = 0;
    for (const auto& v : x) { double d = (double)v - mean; var += d * d; }
    var /= n;
    if (var < 1e-12) return 0.0;
    double cov = 0;
    int m = n - k;
    for (int i = k; i < n; i++)
        cov += ((double)x[i] - mean) * ((double)x[i - k] - mean);
    cov /= m;
    return cov / var;
}

// ================================================================
// Part 3: Task 1 — Bootstrap
// ================================================================

struct BootstrapCI {
    double mean = 0, std_dev = 0;
    double ci_lo_95 = 0, ci_hi_95 = 0;
    double ci_lo_bonf = 0, ci_hi_bonf = 0;
    bool sig_bonf = false;
    double rho_hat = 0;
};

vector<BootstrapCI> parallel_bootstrap_multi_fast(
    const vector<int>& gaps, int K,
    int block_size, int B, int n_threads
) {
    int n = (int)gaps.size();
    vector<BootstrapCI> results(K + 1);
    if (n < block_size + 10) return results;
    int n_blocks = n - block_size + 1;

    auto obs = compute_rho_vector_fast(gaps, K);

    vector<vector<vector<double>>> per_thread(n_threads);
    vector<thread> threads;

    for (int t = 0; t < n_threads; t++) {
        threads.emplace_back([&, t]() {
            mt19937_64 rng(SEED_BOOT + (uint64_t)t * 7919ULL);
            uniform_int_distribution<int> bd(0, n_blocks - 1);
            vector<int> resampled;
            resampled.reserve(n);
            int per = B / n_threads;
            int extra = (t < (B % n_threads)) ? 1 : 0;
            int my_B = per + extra;
            per_thread[t].resize(K + 1);
            for (int k = 1; k <= K; k++) per_thread[t][k].reserve(my_B);
            for (int b = 0; b < my_B; b++) {
                resampled.clear();
                while ((int)resampled.size() < n) {
                    int s = bd(rng);
                    for (int j = 0; j < block_size && (int)resampled.size() < n; j++)
                        resampled.push_back(gaps[s + j]);
                }
                auto rhos = compute_rho_vector_fast(resampled, K);
                for (int k = 1; k <= K; k++)
                    per_thread[t][k].push_back(rhos[k]);
            }
        });
    }
    for (auto& th : threads) th.join();

    double alpha_b = 0.05 / (2.0 * K);
    for (int k = 1; k <= K; k++) {
        vector<double> s;
        s.reserve(B);
        for (int t = 0; t < n_threads; t++)
            for (double v : per_thread[t][k]) s.push_back(v);
        sort(s.begin(), s.end());
        int m = (int)s.size();
        double sum = 0; for (double v : s) sum += v;
        results[k].mean = sum / m;
        double var = 0;
        for (double v : s) var += (v - results[k].mean) * (v - results[k].mean);
        results[k].std_dev = sqrt(var / m);
        int lo95 = max(0, (int)(0.025 * m));
        int hi95 = min(m - 1, (int)(0.975 * m));
        int lo_b = max(0, (int)(alpha_b * m));
        int hi_b = min(m - 1, (int)((1 - alpha_b) * m));
        results[k].ci_lo_95 = s[lo95];
        results[k].ci_hi_95 = s[hi95];
        results[k].ci_lo_bonf = s[lo_b];
        results[k].ci_hi_bonf = s[hi_b];
        results[k].sig_bonf = (results[k].ci_lo_bonf > 0) ||
                              (results[k].ci_hi_bonf < 0);
        results[k].rho_hat = obs[k];
    }
    return results;
}

void task1_bootstrap(const vector<int>& gaps) {
    print_sep();
    cout << "\n[TASK 1] Bootstrap rho_k (k=1.." << K_EXTENDED
         << ", B=" << N_BOOT << ")\n\n";
    int block_size = max(100, (int)pow((double)gaps.size(), 1.0/3.0));
    cout << "  Block size = " << block_size
         << ", threads = " << N_THREADS << "\n";

    auto t0 = chrono::high_resolution_clock::now();
    auto boot = parallel_bootstrap_multi_fast(gaps, K_EXTENDED,
                                               block_size, N_BOOT, N_THREADS);
    double dt = chrono::duration<double>(
        chrono::high_resolution_clock::now() - t0).count();
    cout << "  Elapsed: " << setprecision(2) << dt << " s\n\n";

    ofstream fout("bootstrap_v9_1.csv");
    fout << "# seeds: perm=" << SEED_PERM << " boot=" << SEED_BOOT << "\n";
    fout << "k,rho_hat,boot_mean,boot_std,ci_lo_95,ci_hi_95,"
            "ci_lo_bonf,ci_hi_bonf,sig_bonf\n";
    for (int k = 1; k <= K_EXTENDED; k++) {
        fout << k << "," << setprecision(10) << boot[k].rho_hat << ","
             << boot[k].mean << "," << boot[k].std_dev << ","
             << boot[k].ci_lo_95 << "," << boot[k].ci_hi_95 << ","
             << boot[k].ci_lo_bonf << "," << boot[k].ci_hi_bonf << ","
             << (boot[k].sig_bonf ? 1 : 0) << "\n";
    }
    fout.close();

    cout << setw(5) << "k" << "|"
         << setw(14) << "rho_hat" << "|"
         << setw(24) << "Bonferroni CI" << "|"
         << setw(8) << "sig" << "\n";
    print_sep(60, '-');
    for (int k = 1; k <= min(K_EXTENDED, 20); k++) {
        cout << setw(5) << k << "|"
             << setw(14) << setprecision(6) << boot[k].rho_hat << "|"
             << "[" << setw(9) << boot[k].ci_lo_bonf
             << "," << setw(9) << boot[k].ci_hi_bonf << "]  |"
             << setw(8) << (boot[k].sig_bonf ? "YES" : "no") << "\n";
    }
    cout << "  ... (see bootstrap_v9_1.csv for all 64 lags)\n";
    cout << "  Saved: bootstrap_v9_1.csv\n";
}

// ================================================================
// Part 4: Task 3 — Permutation
// ================================================================

vector<vector<double>> parallel_permutation_fast(
    const vector<int>& gaps, int K, int n_perm, int n_threads
) {
    vector<vector<vector<double>>> per_thread(n_threads);
    atomic<int> counter{0};
    vector<thread> threads;

    for (int t = 0; t < n_threads; t++) {
        threads.emplace_back([&, t]() {
            mt19937_64 rng(SEED_PERM + (uint64_t)t * 1000003ULL);
            vector<int> shuffled = gaps;
            int per = n_perm / n_threads;
            int extra = (t < (n_perm % n_threads)) ? 1 : 0;
            int my_n = per + extra;
            per_thread[t].reserve(my_n);
            int sz = (int)shuffled.size();
            for (int i = 0; i < my_n; i++) {
                for (int j = sz - 1; j > 0; j--) {
                    uniform_int_distribution<int> d(0, j);
                    swap(shuffled[j], shuffled[d(rng)]);
                }
                per_thread[t].push_back(compute_rho_vector_fast(shuffled, K));
                counter.fetch_add(1, memory_order_relaxed);
            }
        });
    }
    int last = 0;
    while (counter.load() < n_perm) {
        this_thread::sleep_for(chrono::seconds(5));
        int cur = counter.load();
        if (cur != last) {
            cerr << "\r  Permutations: " << cur << "/" << n_perm
                 << " (" << (100.0 * cur / n_perm) << "%)    " << flush;
            last = cur;
        }
    }
    for (auto& th : threads) th.join();
    cerr << "\r  Permutations: " << n_perm << "/" << n_perm << " (100%)    \n";

    vector<vector<double>> null_dist(K + 1);
    for (int t = 0; t < n_threads; t++)
        for (auto& rhos : per_thread[t])
            for (int k = 1; k <= K; k++)
                null_dist[k].push_back(rhos[k]);
    return null_dist;
}

void task3_permutation(const vector<int>& gaps) {
    print_sep();
    cout << "\n[TASK 3] Permutation test (k=1.." << K_PERM
         << ", P=" << N_PERM << ")\n\n";
    auto obs = compute_rho_vector_fast(gaps, K_PERM);
    auto t0 = chrono::high_resolution_clock::now();
    auto null_dist = parallel_permutation_fast(gaps, K_PERM, N_PERM, N_THREADS);
    double dt = chrono::duration<double>(
        chrono::high_resolution_clock::now() - t0).count();
    cout << "  Elapsed: " << setprecision(2) << dt << " s\n\n";

    double alpha_bonf = 0.05 / K_PERM;
    vector<double> p_vals(K_PERM + 1);
    for (int k = 1; k <= K_PERM; k++) {
        auto& null = null_dist[k];
        sort(null.begin(), null.end());
        int extreme = 0;
        for (double v : null)
            if (fabs(v) >= fabs(obs[k])) extreme++;
        p_vals[k] = (double)(extreme + 1) / (N_PERM + 1);
    }

    vector<pair<double,int>> sorted_p;
    for (int k = 1; k <= K_PERM; k++) sorted_p.push_back({p_vals[k], k});
    sort(sorted_p.begin(), sorted_p.end());
    vector<bool> bh_sig(K_PERM + 1, false);
    for (int rank = K_PERM; rank >= 1; rank--) {
        int k = sorted_p[rank - 1].second;
        if (p_vals[k] <= (double)rank / K_PERM * 0.05) {
            for (int r = 0; r < rank; r++) bh_sig[sorted_p[r].second] = true;
            break;
        }
    }

    ofstream fout("permutation_v9_1.csv");
    fout << "# seed: " << SEED_PERM << ", P=" << N_PERM << "\n";
    fout << "k,rho_obs,null_mean,null_std,p_value,sig_bonf,sig_fdr\n";
    for (int k = 1; k <= K_PERM; k++) {
        auto& null = null_dist[k];
        sort(null.begin(), null.end());
        double sum = 0; for (double v : null) sum += v;
        double nm = sum / N_PERM;
        double var = 0;
        for (double v : null) var += (v - nm) * (v - nm);
        double ns = sqrt(var / N_PERM);
        fout << k << "," << setprecision(10) << obs[k] << ","
             << nm << "," << ns << "," << p_vals[k] << ","
             << (p_vals[k] < alpha_bonf ? 1 : 0) << ","
             << (bh_sig[k] ? 1 : 0) << "\n";
    }
    fout.close();

    cout << setw(5) << "k" << "|"
         << setw(12) << "rho_obs" << "|"
         << setw(12) << "p_value" << "|"
         << setw(6) << "Bonf" << "|"
         << setw(6) << "FDR" << "\n";
    print_sep(55, '-');
    for (int k = 1; k <= K_PERM; k++) {
        cout << setw(5) << k << "|"
             << setw(12) << setprecision(6) << obs[k] << "|"
             << setw(12) << setprecision(6) << p_vals[k] << "|"
             << setw(6) << (p_vals[k] < alpha_bonf ? "YES" : "no") << "|"
             << setw(6) << (bh_sig[k] ? "YES" : "no") << "\n";
    }
    cout << "  Saved: permutation_v9_1.csv\n";
}

// ================================================================
// Part 5: Task 7 — mod-3
// ================================================================

void task7_mod3(const vector<ll>& primes) {
    print_sep();
    cout << "\n[TASK 7] mod-3 conditional rho_k (closing prime)\n\n";
    vector<int> gaps_c1, gaps_c2;
    for (size_t i = 0; i + 1 < primes.size(); i++) {
        ll p_closing = primes[i + 1];
        if (p_closing < 5) continue;
        int gap = (int)(primes[i + 1] - primes[i]);
        if (gap <= 1) continue;
        if (p_closing % 3 == 1) gaps_c1.push_back(gap);
        else if (p_closing % 3 == 2) gaps_c2.push_back(gap);
    }
    cout << "  Class 1 (p_{n+1} ≡ 1 mod 3): n = " << gaps_c1.size() << "\n";
    cout << "  Class 2 (p_{n+1} ≡ 2 mod 3): n = " << gaps_c2.size() << "\n";

    auto rho1 = compute_rho_vector_fast(gaps_c1, K_SHORT);
    auto rho2 = compute_rho_vector_fast(gaps_c2, K_SHORT);

    ofstream fout("mod3_v9_1.csv");
    fout << "k,rho_c1,rho_c2,diff\n";
    for (int k = 1; k <= K_SHORT; k++)
        fout << k << "," << rho1[k] << "," << rho2[k] << ","
             << (rho1[k] - rho2[k]) << "\n";
    fout.close();

    cout << setw(5) << "k" << "|"
         << setw(14) << "rho_c1" << "|"
         << setw(14) << "rho_c2" << "|"
         << setw(12) << "diff" << "\n";
    print_sep(50, '-');
    for (int k = 1; k <= K_SHORT; k++) {
        cout << setw(5) << k << "|"
             << setw(14) << setprecision(6) << rho1[k] << "|"
             << setw(14) << rho2[k] << "|"
             << setw(12) << (rho1[k] - rho2[k]) << "\n";
    }
    cout << "  Saved: mod3_v9_1.csv\n";
}

// ================================================================
// Part 6: Task 8 — windowed correlogram
// ================================================================

pair<double,double> window_bootstrap_rho1(
    const vector<double>& detrended, int block_size, int B
) {
    int n = (int)detrended.size();
    if (n < block_size + 10) return {0.0, 0.0};
    int n_blocks = n - block_size + 1;
    mt19937_64 rng(SEED_WIN_BOOT);
    uniform_int_distribution<int> bd(0, n_blocks - 1);
    vector<double> resampled;
    resampled.reserve(n);
    vector<double> samples;
    samples.reserve(B);
    for (int b = 0; b < B; b++) {
        resampled.clear();
        while ((int)resampled.size() < n) {
            int s = bd(rng);
            for (int j = 0; j < block_size && (int)resampled.size() < n; j++)
                resampled.push_back(detrended[s + j]);
        }
        samples.push_back(compute_rho_single_fast(resampled, 1));
    }
    sort(samples.begin(), samples.end());
    return {samples[(int)(0.025 * B)], samples[(int)(0.975 * B)]};
}

void task8_correlogram(
    const vector<ll>& primes,
    const vector<pair<ll,ll>>& windows
) {
    print_sep();
    cout << "\n[TASK 8] Per-window correlogram (k=1.." << K_EXTENDED << ")\n\n";
    ofstream fout("correlogram_v9_1.csv");
    fout << "# window_lo,window_hi,n_gaps,mean_gap,rho1_ci_lo,rho1_ci_hi";
    for (int k = 1; k <= K_EXTENDED; k++) fout << ",rho_" << k;
    fout << "\n";

    for (auto& [lo, hi] : windows) {
        vector<int> wgaps;
        vector<double> lnp;
        for (size_t i = 0; i + 1 < primes.size(); i++) {
            if (primes[i] >= lo && primes[i] < hi) {
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

        auto rho = compute_rho_vector_fast(dt, K_EXTENDED);
        auto ci = window_bootstrap_rho1(dt, 179, 1000);
        bool sig = (ci.first > 0) || (ci.second < 0);

        fout << lo << "," << hi << "," << wgaps.size() << ","
             << mean_g << "," << ci.first << "," << ci.second;
        cout << "[" << lo << "," << hi << "): n=" << wgaps.size()
             << "  rho_1=" << setprecision(5) << rho[1]
             << "  CI=[" << setprecision(5) << ci.first << ","
             << ci.second << "]" << (sig ? " *" : "") << "\n";
        for (int k = 1; k <= K_EXTENDED; k++)
            fout << "," << setprecision(8) << rho[k];
        fout << "\n";
    }
    fout.close();
    cout << "  Saved: correlogram_v9_1.csv\n";
}

// ================================================================
// Part 7: Task 9 — Inclusion-exclusion model (with cache)
// ================================================================

struct SingularSeries {
    vector<int> primes;
    vector<double> log_tail;      // finite product, p in (30, cutoff]
    int cutoff;

    explicit SingularSeries(int prime_cutoff) : cutoff(prime_cutoff) {
        vector<bool> sv(prime_cutoff + 1, true);
        sv[0] = sv[1] = false;
        for (int i = 2; (ll)i * i <= prime_cutoff; i++)
            if (sv[i])
                for (ll j = (ll)i * i; j <= prime_cutoff; j += i)
                    sv[(size_t)j] = false;
        for (int i = 2; i <= prime_cutoff; i++)
            if (sv[i]) primes.push_back(i);

        log_tail.assign(K_MAX_IE + 1, 0.0);
        for (int k = 1; k <= K_MAX_IE; k++) {
            double s = 0.0;
            for (int p : primes) {
                if (p <= 30) continue;
                s += log1p(-(double)k / p) - k * log1p(-1.0 / p);
            }
            log_tail[k] = s;
        }
    }

    // S() is the direct (uncached) computation.
    // Cached version is S_cached() below.
    double S(const vector<int>& pattern) const {
        int k = (int)pattern.size();
        if (k > K_MAX_IE) {
            throw runtime_error(
                "Pattern size " + to_string(k) +
                " exceeds K_MAX_IE=" + to_string(K_MAX_IE)
            );
        }
        double r = exp(log_tail[k]);
        for (int p : primes) {
            if (p > 30) break;
            set<int> residues;
            for (int h : pattern)
                residues.insert(((h % p) + p) % p);
            int nu = (int)residues.size();
            if (nu == p) return 0.0;  // inadmissible
            r *= (1.0 - (double)nu / p) / pow(1.0 - 1.0 / p, k);
        }
        return r;
    }
};

// ---------- Cache for S() ----------
struct PatternHash {
    size_t operator()(const vector<int>& v) const {
        size_t h = 1469598103934665603ULL;
        for (int x : v) {
            h ^= (size_t)(x + 0x9e3779b9);
            h *= 1099511628211ULL;
        }
        return h;
    }
};
static unordered_map<vector<int>, double, PatternHash> g_S_cache;

double S_cached(const SingularSeries& ss, const vector<int>& pattern) {
    auto it = g_S_cache.find(pattern);
    if (it != g_S_cache.end()) return it->second;
    double val = ss.S(pattern);
    g_S_cache[pattern] = val;
    return val;
}

// Inclusion-exclusion coefficients
vector<double> ie_coeffs(const SingularSeries& ss,
                         const vector<int>& base, int span,
                         int excl, int R_max)
{
    vector<int> E;
    for (int h = 2; h < span; h += 2)
        if (h != excl) E.push_back(h);

    int r_top = min(R_max, (int)E.size());
    vector<double> e(r_top + 1, 0.0);

    for (int r = 0; r <= r_top; r++) {
        if (r == 0) {
            e[0] = S_cached(ss, base);
            continue;
        }
        vector<int> idx(r);
        for (int i = 0; i < r; i++) idx[i] = i;
        double s = 0.0;
        while (true) {
            vector<int> pat = base;
            for (int i = 0; i < r; i++) pat.push_back(E[idx[i]]);
            s += S_cached(ss, pat);
            int i = r - 1;
            while (i >= 0 && idx[i] == (int)E.size() - r + i) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < r; j++) idx[j] = idx[j-1] + 1;
        }
        e[r] = s;
    }
    return e;
}

double ie_D(const SingularSeries& ss,
            const vector<int>& base, int span, int excl,
            double log_x, int R_max)
{
    auto e = ie_coeffs(ss, base, span, excl, R_max);
    int k = (int)base.size();
    double result = 0.0;
    for (int r = 0; r < (int)e.size(); r++) {
        double term = e[r] / pow(log_x, k + r);
        result += (r % 2 == 0 ? term : -term);
    }
    return result;
}

double ie_R_model(const SingularSeries& ss, int a, int b,
                  double log_x, int R_max)
{
    vector<int> base_a = {0, a};
    vector<int> base_b = {0, b};
    vector<int> base_ab = {0, a, a + b};
    double Da = ie_D(ss, base_a, a, -1, log_x, R_max);
    double Db = ie_D(ss, base_b, b, -1, log_x, R_max);
    double Dab = ie_D(ss, base_ab, a + b, a, log_x, R_max);
    if (Da * Db == 0.0) return 0.0;
    return Dab / (log_x * Da * Db);
}

double ie_R_HL(const SingularSeries& ss, int a, int b) {
    vector<int> base_a = {0, a};
    vector<int> base_b = {0, b};
    vector<int> base_ab = {0, a, a + b};
    double Sa = S_cached(ss, base_a);
    double Sb = S_cached(ss, base_b);
    double Sab = S_cached(ss, base_ab);
    if (Sa * Sb < 1e-15) return 0.0;
    return Sab / (Sa * Sb);
}

// ================================================================
// Part 8: Task 6 — save marginals and pair counts
// ================================================================

void task6_save_marginals(
    const map<int, ll>& gap_dist,
    const map<pair<int,int>, ll>& gap_pairs,
    ll total_gaps
) {
    print_sep();
    cout << "\n[TASK 6] HL marginals and pair counts\n\n";

    ofstream fm("marginal_v9_1.csv");
    fm << "gap,count\n";
    for (auto& kv : gap_dist)
        fm << kv.first << "," << kv.second << "\n";
    fm.close();

    ofstream fp("pair_counts_v9_1.csv");
    fp << "a,b,count\n";
    for (auto& kv : gap_pairs)
        fp << kv.first.first << "," << kv.first.second << ","
           << kv.second << "\n";
    fp.close();

    cout << "  Marginals saved: marginal_v9_1.csv (" << gap_dist.size()
         << " entries)\n";
    cout << "  Pair counts saved: pair_counts_v9_1.csv ("
         << gap_pairs.size() << " entries)\n";
    cout << "  Total gaps = " << total_gaps << "\n";

    // R_obs recomputed from counts
    ll N_pairs = total_gaps - 1;
    ofstream fr("R_obs_v9_1.csv");
    fr << "a,b,N_ab,R_obs\n";
    for (auto& pr : PAIRS_17) {
        int a = pr.first, b = pr.second;
        auto it = gap_pairs.find({a, b});
        if (it == gap_pairs.end()) continue;
        auto ia = gap_dist.find(a);
        auto ib = gap_dist.find(b);
        if (ia == gap_dist.end() || ib == gap_dist.end()) continue;
        ll cnt = it->second;
        double p_ab = (double)cnt / N_pairs;
        double p_a = (double)ia->second / total_gaps;
        double p_b = (double)ib->second / total_gaps;
        double R_obs = p_ab / (p_a * p_b);
        fr << a << "," << b << "," << cnt << ","
           << setprecision(10) << R_obs << "\n";
    }
    fr.close();
    cout << "  R_obs for 17 pairs saved: R_obs_v9_1.csv\n";
}

// ================================================================
// Part 9: Task 9 — inclusion-exclusion on 17 pairs
// ================================================================

void task9_ie_model(
    const map<int, ll>& gap_dist,
    const map<pair<int,int>, ll>& gap_pairs,
    ll total_gaps, ll N
) {
    print_sep();
    cout << "\n[TASK 9] Inclusion-exclusion model (Kimi)\n";
    cout << "  R_max = " << IE_R_MAX
         << ", prime cutoff = " << IE_PRIME_LIMIT << "\n\n";

    cout << "  Building singular series...\n";
    auto t0 = chrono::high_resolution_clock::now();
    SingularSeries ss(IE_PRIME_LIMIT);
    double dt_build = chrono::duration<double>(
        chrono::high_resolution_clock::now() - t0).count();
    cout << "  Built in " << setprecision(2) << dt_build << " s ("
         << ss.primes.size() << " primes)\n\n";

    double log_x = log((double)N);
    ll N_pairs = total_gaps - 1;

    ofstream fout("hl_pairs_v9_1.csv");
    fout << "# N = " << N << "\n";
    fout << "# log_x = " << setprecision(10) << log_x << "\n";
    fout << "# R_max = " << IE_R_MAX << "\n";
    fout << "# singular_series_finite_product_cutoff = "
         << IE_PRIME_LIMIT << "\n";
    fout << "# primes_above_cutoff_not_included = true\n";
    fout << "# analytic_tail_correction = none\n";
    fout << "a,b,N_ab,R_obs,R_HL,R_model,rel_err_HL,rel_err_model\n";

    cout << setw(10) << "Pair" << "|"
         << setw(10) << "R_obs" << "|"
         << setw(10) << "R_HL" << "|"
         << setw(10) << "R_model" << "|"
         << setw(10) << "rel_HL" << "|"
         << setw(10) << "rel_md" << "\n";
    print_sep(70, '-');

    double sum_abs_HL = 0, sum_abs_md = 0;
    int cnt = 0;

    for (auto& pr : PAIRS_17) {
        int a = pr.first, b = pr.second;
        auto it = gap_pairs.find({a, b});
        if (it == gap_pairs.end()) continue;
        auto ia = gap_dist.find(a);
        auto ib = gap_dist.find(b);
        if (ia == gap_dist.end() || ib == gap_dist.end()) continue;

        ll cnt_ab = it->second;
        double p_ab = (double)cnt_ab / N_pairs;
        double p_a = (double)ia->second / total_gaps;
        double p_b = (double)ib->second / total_gaps;
        double R_obs = p_ab / (p_a * p_b);

        double R_HL_ = ie_R_HL(ss, a, b);
        double R_md = ie_R_model(ss, a, b, log_x, IE_R_MAX);

        double rel_HL = (R_HL_ - R_obs) / R_obs;
        double rel_md = (R_md - R_obs) / R_obs;

        cout << "(" << setw(3) << a << "," << setw(3) << b << ")  |"
             << setw(10) << setprecision(5) << R_obs << "|"
             << setw(10) << R_HL_ << "|"
             << setw(10) << R_md << "|"
             << setw(10) << setprecision(4) << rel_HL << "|"
             << setw(10) << rel_md << "\n";

        fout << a << "," << b << "," << cnt_ab << ","
             << setprecision(10) << R_obs << ","
             << R_HL_ << ","
             << R_md << ","
             << scientific << rel_HL << ","
             << rel_md << "\n";

        sum_abs_HL += fabs(rel_HL);
        sum_abs_md += fabs(rel_md);
        cnt++;
    }
    fout.close();

    cout << "\n  Summary (mean |relative error|):\n";
    cout << "    R_HL    : " << setprecision(4)
         << (sum_abs_HL / cnt) * 100 << "%\n";
    cout << "    R_model : " << (sum_abs_md / cnt) * 100 << "%\n";
    cout << "  Cache size = " << g_S_cache.size() << " patterns\n";
    cout << "  Saved: hl_pairs_v9_1.csv\n";
}

// ================================================================
// Part 10: Consistency checks
// ================================================================

void run_consistency_checks(
    const vector<int>& all_gaps,
    const vector<ll>& all_primes,
    const map<int, ll>& gap_dist,
    const map<pair<int,int>, ll>& gap_pairs,
    ll N
) {
    print_sep();
    cout << "\n[CHECKS] Consistency verification\n\n";

    int n_fail = 0;

    // Check 1: gaps = primes - 1 (excluding initial gap=1)
    ll prime_count = (ll)all_primes.size();
    ll expected_gaps = prime_count - 2;  // exclude first gap (2->3 = 1)
    cout << "  Check 1: gaps = " << all_gaps.size()
         << " (expected " << expected_gaps << ")";
    if ((ll)all_gaps.size() == expected_gaps) cout << "  PASS\n";
    else { cout << "  FAIL\n"; n_fail++; }

    // Check 2: total pair count
    ll total_pairs = 0;
    for (auto& kv : gap_pairs) total_pairs += kv.second;
    ll expected_pairs = (ll)all_gaps.size() - 1;
    cout << "  Check 2: pairs = " << total_pairs
         << " (expected " << expected_pairs << ")";
    if (total_pairs == expected_pairs) cout << "  PASS\n";
    else { cout << "  FAIL\n"; n_fail++; }

    // Check 3: marginal consistency
    map<int, ll> marginal_from_pairs;
    for (auto& kv : gap_pairs)
        marginal_from_pairs[kv.first.first] += kv.second;
    int mismatches = 0;
    for (auto& kv : marginal_from_pairs) {
        int a = kv.first;
        ll sum_b = kv.second;
        ll N_a = gap_dist.at(a);
        // N_a - sum_b = 1 if a appears as last gap, 0 otherwise
        if (sum_b > N_a || N_a - sum_b > 1) mismatches++;
    }
    cout << "  Check 3: marginal consistency (sum_b N_ab = N_a ± 1)";
    if (mismatches == 0) cout << "  PASS\n";
    else { cout << "  FAIL (" << mismatches << " mismatches)\n"; n_fail++; }

    // Check 4: known prime count at N=10^9
    if (N == 10000000LL) {
        cout << "  Check 4: primes = " << prime_count
             << " (expected 50847534)";
        if (prime_count == 50847534LL) cout << "  PASS\n";
        else { cout << "  FAIL\n"; n_fail++; }

        cout << "  Check 5: gaps = " << all_gaps.size()
             << " (expected 50847532)";
        if ((ll)all_gaps.size() == 50847532LL) cout << "  PASS\n";
        else { cout << "  FAIL\n"; n_fail++; }
    }

    cout << "\n  Total failures: " << n_fail << "\n";
    if (n_fail == 0) cout << "  All checks PASSED.\n";
    else cout << "  *** Some checks FAILED — investigate before trusting results. ***\n";
}

// ================================================================
// Main
// ================================================================

int main(int argc, char* argv[]) {
    ll N = DEFAULT_N;
    if (argc > 1) {
        try { N = stoll(argv[1]); }
        catch (...) { cerr << "Invalid N. Using default.\n"; }
    }
    if (N < 1000) { cerr << "N too small.\n"; return 1; }

    ofstream log_file("results_v9_1.txt");
    streambuf* original_cout = cout.rdbuf();
    TeeBuf tee(original_cout, log_file.rdbuf());
    cout.rdbuf(&tee);

    auto t_start = chrono::high_resolution_clock::now();
    cout << fixed << setprecision(6);

    print_sep();
    cout << "  PRIME WHEELS & PRIME GAPS — v9.1 (complete, corrected)\n";
    cout << "  N = " << N << "\n";
    cout << "  Threads = " << N_THREADS << "\n";
    cout << "  K_PERM = " << K_PERM << ", K_EXTENDED = " << K_EXTENDED << "\n";
    cout << "  N_PERM = " << N_PERM << ", N_BOOT = " << N_BOOT << "\n";
    cout << "  IE: R_max = " << IE_R_MAX
         << ", prime cutoff = " << IE_PRIME_LIMIT << "\n";
    print_sep();

    // ---------------- Phase 0: Sieve ----------------
    cout << "\n[Phase 0] Sieving...\n";
    vector<double> edges = {1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6,
                            3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10};

    vector<int> all_gaps;
    vector<ll> all_primes;
    all_gaps.reserve((size_t)(N / 15));
    all_primes.reserve((size_t)(N / 15));

    map<int, ll> gap_dist;
    map<pair<int,int>, ll> gap_pairs;
    ll prev = -1, prev_gap = -1, prime_count = 0, last_rep = 0;
    ld sum_g = 0;
    bool first_gap_seen = false;

    PrimeGenerator gen(N);
    gen.generate_all([&](ll p) {
        prime_count++;
        all_primes.push_back(p);
        if (prev > 0) {
            int gap = (int)(p - prev);
            // استبعاد الفجوة الأولى (2→3 = 1) صراحةً
            if (gap > 1) {
                all_gaps.push_back(gap);
                gap_dist[gap]++;
                sum_g += gap;
                // بناء gap_pairs مباشرة
                if (first_gap_seen) {
                    gap_pairs[{prev_gap, gap}]++;
                }
                prev_gap = gap;
                first_gap_seen = true;
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
    cout << "  g_mean = " << setprecision(4)
         << (double)(sum_g / all_gaps.size()) << "\n";
    cout << "  Time:   " << setprecision(2)
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_start).count()
         << " s\n";

    // ---------------- Tasks ----------------
    auto t_task = chrono::high_resolution_clock::now();
    task1_bootstrap(all_gaps);
    cerr << "  [Task 1 done in "
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_task).count()
         << " s]\n";

    t_task = chrono::high_resolution_clock::now();
    task3_permutation(all_gaps);
    cerr << "  [Task 3 done in "
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_task).count()
         << " s]\n";

    t_task = chrono::high_resolution_clock::now();
    task6_save_marginals(gap_dist, gap_pairs, (ll)all_gaps.size());
    cerr << "  [Task 6 done in "
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_task).count()
         << " s]\n";

    t_task = chrono::high_resolution_clock::now();
    task7_mod3(all_primes);
    cerr << "  [Task 7 done in "
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_task).count()
         << " s]\n";

    t_task = chrono::high_resolution_clock::now();
    vector<pair<ll,ll>> log_windows;
    for (size_t i = 0; i + 1 < edges.size(); i++) {
        if ((ll)edges[i] >= N) break;
        log_windows.push_back({(ll)edges[i], min((ll)edges[i+1], N)});
    }
    task8_correlogram(all_primes, log_windows);
    cerr << "  [Task 8 done in "
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_task).count()
         << " s]\n";

    t_task = chrono::high_resolution_clock::now();
    task9_ie_model(gap_dist, gap_pairs, (ll)all_gaps.size(), N);
    cerr << "  [Task 9 done in "
         << chrono::duration<double>(
                chrono::high_resolution_clock::now() - t_task).count()
         << " s]\n";

    // ---------------- Consistency checks ----------------
    run_consistency_checks(all_gaps, all_primes, gap_dist, gap_pairs, N);

    // ---------------- Summary ----------------
    print_sep();
    cout << "\n[FINAL] Summary\n\n";
    cout << "  N = " << N << "\n";
    cout << "  primes = " << prime_count << "\n";
    cout << "  gaps = " << all_gaps.size() << "\n";
    cout << "  g_mean = " << setprecision(4)
         << (double)(sum_g / all_gaps.size()) << "\n";
    auto t_end = chrono::high_resolution_clock::now();
    cout << "  Total time: " << setprecision(2)
         << chrono::duration<double>(t_end - t_start).count() << " s\n";
    cout << "\n  Files saved:\n";
    cout << "    - results_v9_1.txt\n";
    cout << "    - bootstrap_v9_1.csv\n";
    cout << "    - permutation_v9_1.csv\n";
    cout << "    - marginal_v9_1.csv\n";
    cout << "    - pair_counts_v9_1.csv\n";
    cout << "    - R_obs_v9_1.csv\n";
    cout << "    - mod3_v9_1.csv\n";
    cout << "    - correlogram_v9_1.csv\n";
    cout << "    - hl_pairs_v9_1.csv\n";

    cout.flush();
    log_file.flush();
    cout.rdbuf(original_cout);

    return 0;
}
