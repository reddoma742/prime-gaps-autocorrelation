// test_1e10_final_v3.cpp
// All Claude corrections incorporated:
//   (1) Odd-only sieve in number domain (clear, unambiguous)
//   (2) closing_prime recorded explicitly
//   (3) isqrt instead of sqrt+1000
//   (4) Explicit throw after segment if window not closed
//   (5) Small-scale validation (10^6, 10^7) before large run
// Compile: g++ -O2 -std=c++17 test_1e10_final_v3.cpp -o test_1e10_final

#include <bits/stdc++.h>
using namespace std;

const uint64_t LO = 10000000000ULL;   // 1e10
const uint64_t HI = 30000000000ULL;   // 3e10

// ---------------------------------------------------------------------------
uint64_t isqrt(uint64_t n) {
    if (n == 0) return 0;
    uint64_t r = (uint64_t)std::sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

uint64_t phi_of(uint64_t M) {
    uint64_t res = M, m = M;
    for (uint64_t p = 2; p*p <= m; ++p) if (m % p == 0) {
        res = res/p*(p-1);
        while (m % p == 0) m /= p;
    }
    if (m > 1) res = res/m*(m-1);
    return res;
}

vector<uint32_t> base_primes(uint32_t limit) {
    vector<bool> s(limit + 1, true);
    s[0] = s[1] = false;
    for (uint32_t i = 2; (uint64_t)i*i <= limit; ++i)
        if (s[i])
            for (uint64_t j = (uint64_t)i*i; j <= limit; j += i)
                s[(size_t)j] = false;
    vector<uint32_t> out;
    for (uint32_t i = 2; i <= limit; ++i)
        if (s[i]) out.push_back(i);
    return out;
}

// ---------------------------------------------------------------------------
// Odd-only segmented sieve.
// Window: gaps with OPENING prime in [lo, hi). Closing prime may be > hi.
// Records: first_prime, last_open, closing_prime, min/max gap, count.
// ---------------------------------------------------------------------------
template<class F>
void stream_gaps(uint64_t lo, uint64_t hi, F callback) {
    if (lo < 3) lo = 3;
    if (lo % 2 == 0) lo++;

    uint64_t root = isqrt(hi) + 100;   // small safety margin
    vector<uint32_t> base = base_primes((uint32_t)root);

    const uint64_t SEG_ODDS = 50000000ULL;

    uint64_t prev_prime = 0;
    bool closed_last_window_gap = false;
    uint64_t closing_prime = 0;

    for (uint64_t s = lo; !closed_last_window_gap; s += 2 * SEG_ODDS) {
        uint64_t seg_odd_end = s + 2 * SEG_ODDS;
        vector<bool> flag(SEG_ODDS, true);

        for (uint32_t q : base) {
            if (q == 2) continue;   // odd-only

            uint64_t q2 = (uint64_t)q * q;
            if (q2 >= seg_odd_end) break;

            // first odd multiple of q, >= max(q�, s)
            uint64_t first = max(q2, ((s + q - 1) / q) * q);
            if ((first & 1ULL) == 0) first += q;   // make it odd
            if (first >= seg_odd_end) continue;

            // iterate in NUMBER domain, step = 2*q (odd multiples only)
            for (uint64_t j = first; j < seg_odd_end; j += 2ULL * q) {
                uint64_t idx = (j - s) / 2;
                flag[(size_t)idx] = false;
            }
        }

        for (uint64_t j = 0; j < SEG_ODDS; ++j) {
            if (!flag[(size_t)j]) continue;
            uint64_t p = s + 2 * j;
            if (p < lo) { prev_prime = p; continue; }

            if (prev_prime > 0 && prev_prime >= lo && prev_prime < hi) {
                callback(prev_prime, p - prev_prime);
                if (p > hi) {
                    closing_prime = p;
                    closed_last_window_gap = true;
                }
            }
            prev_prime = p;
            if (closed_last_window_gap) break;
        }

        // explicit check after each segment
        if (s + 2 * SEG_ODDS > hi && !closed_last_window_gap) {
            // we passed hi but haven't found a prime > hi yet � continue
            // (primes are infinite, this always terminates)
        }
    }

    if (!closed_last_window_gap)
        throw runtime_error("Window boundary was not closed");
    if (closing_prime <= hi)
        throw runtime_error("closing_prime <= HI; logical error");
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main() {
    auto T0 = chrono::steady_clock::now();
    ofstream fout("test_1e10_final_v3.txt");
    auto P = [&](const string& s){ cout << s; fout << s; };

    // ========== PART 0: Small-scale validation ==========
    cerr << "\n=== PART 0: Validation at 10^6 and 10^7 ===\n";

    struct ValCase { uint64_t lo, hi; uint64_t expected_primes; };
    vector<ValCase> val = {
        { 1000000ULL, 10000000ULL, 664579 - 78498 },     // = 586081
        { 10000000ULL, 100000000ULL, 5761455 - 664579 }  // = 5096876
    };

    for (auto& vc : val) {
        uint64_t n = 0;
        uint64_t first_p = 0, last_open = 0;
        stream_gaps(vc.lo, vc.hi, [&](uint64_t p, uint64_t g){
            if (n == 0) first_p = p;
            last_open = p;
            n++;
        });
        ostringstream ss;
        ss << "  [" << vc.lo << ", " << vc.hi << "): "
           << "n = " << n
           << "  expected = " << vc.expected_primes
           << "  " << (n == vc.expected_primes ? "PASS" : "FAIL") << "\n";
        ss << "    first = " << first_p << ", last_open = " << last_open << "\n";
        P(ss.str());
    }

    // ========== PART 1: Large window ==========
    P("\n============================================================\n");
    P("  WINDOW DEFINITION\n");
    P("  Gaps g_n = p_{n+1} - p_n with OPENING prime p_n in [LO, HI).\n");
    P("  Closing prime of last gap may be > HI.\n");
    P("============================================================\n\n");

    vector<uint64_t> Ms = { 6, 30, 210, 2310, 30030, 510510, 9699690 };
    size_t nM = Ms.size();

    // ---- PASS 1 ----
    cerr << "PASS 1: OLS + class sums\n";
    uint64_t n1 = 0;
    long double Sx=0, Sy=0, Sxx=0, Sxy=0, Syy=0;
    uint64_t first_p = 0, last_open = 0;
    uint32_t min_gap = UINT32_MAX, max_gap = 0;

    vector<vector<long double>> sum_g(nM), sum_lp(nM);
    vector<vector<uint32_t>>    cnt(nM);
    vector<uint64_t> phiM(nM);
    for (size_t k = 0; k < nM; ++k) {
        phiM[k] = phi_of(Ms[k]);
        sum_g[k].assign(Ms[k], 0.0L);
        sum_lp[k].assign(Ms[k], 0.0L);
        cnt[k].assign(Ms[k], 0);
    }

    stream_gaps(LO, HI, [&](uint64_t p, uint64_t g) {
        if (n1 == 0) first_p = p;
        last_open = p;
        if (g < min_gap) min_gap = (uint32_t)g;
        if (g > max_gap) max_gap = (uint32_t)g;

        double lp = (double)log((double)p);
        double gd = (double)g;
        n1++;
        Sx += lp; Sy += gd; Sxx += lp*lp; Sxy += lp*gd; Syy += gd*gd;

        for (size_t k = 0; k < nM; ++k) {
            uint32_t r = (uint32_t)(p % Ms[k]);
            sum_g[k][r] += gd;
            sum_lp[k][r] += lp;
            cnt[k][r]++;
        }
    });

    long double dn = (long double)n1;
    long double den = dn*Sxx - Sx*Sx;
    double beta  = (double)((dn*Sxy - Sx*Sy) / den);
    double alpha = (double)((Sy - beta*Sx) / dn);

    vector<vector<double>> mu(nM);
    for (size_t k = 0; k < nM; ++k) {
        mu[k].assign(Ms[k], 0.0);
        uint32_t min_class = UINT32_MAX;
        uint32_t non_empty = 0;
        for (uint64_t r = 0; r < Ms[k]; ++r) {
            if (cnt[k][r] == 0) continue;
            non_empty++;
            if (cnt[k][r] < min_class) min_class = cnt[k][r];
            if (cnt[k][r] < 2)
                throw runtime_error("Sparse class at M=" + to_string(Ms[k]));
            mu[k][r] = (double)((sum_g[k][r] - alpha*cnt[k][r]
                                - beta*sum_lp[k][r]) / (long double)cnt[k][r]);
        }
        cerr << "  M = " << Ms[k] << ": non-empty = " << non_empty
             << ", min class = " << min_class << "\n";
        if (min_class < 50)
            cerr << "  WARNING: sparse classes at M=" << Ms[k] << "\n";
    }
    for (auto& v : sum_g) v.clear();
    for (auto& v : sum_lp) v.clear();

    // ---- OLS decomposition consistency check ----
    long double mean_g   = Sy / dn;
    long double mean_lp  = Sx / dn;
    long double var_g    = (double)(Syy/dn - mean_g*mean_g);
    long double var_lp   = (double)(Sxx/dn - mean_lp*mean_lp);

    double beta_check    = (double)((dn*Sxy - Sx*Sy) / (dn*Sxx - Sx*Sx));
    double var_e_pred    = (double)(var_g - beta_check*beta_check*var_lp);

    cerr << "  Var(g)        = " << (double)var_g << "\n";
    cerr << "  Var(log p)    = " << (double)var_lp << "\n";
    cerr << "  beta          = " << beta_check << "\n";
    cerr << "  Var(g) - b^2*Var(log p) = " << var_e_pred << "\n";

    // ---- PASS 2 ----
    cerr << "PASS 2: Var, Cov, ABCD, rho_res\n";

    long double sum_e = 0, sum_e2 = 0, sum_ee = 0;
    uint64_t n2 = 0;
    double prev_e = 0;

    struct MA {
        long double A=0, B=0, C=0, D=0;
        long double sum_res2 = 0, sum_resres = 0;
        uint64_t cnt_ABCD = 0, cnt_res = 0;
        double prev_mu = 0, prev_res = 0;
        bool has_prev = false;
    };
    vector<MA> acc(nM);

    // ---- Block-level accumulators for the bootstrap -----------------------
    // 40 blocks x 7 M x 5 doubles ~= 11 KB. No residuals are stored.
    const int    N_BLOCKS  = 40;
    const int    N_BOOT    = 1000;
    const uint64_t SEED_BASE = 20261006ULL;   // independent seed per M
    const size_t block_size = (size_t)n1 / N_BLOCKS;   // n1 == n2 (asserted below)
    if (block_size == 0) throw runtime_error("block_size == 0");
    vector<array<double,40>> block_sum_u (nM), block_sum_u2(nM),
                             block_sum_uu(nM), block_first (nM), block_last (nM);
    for (size_t k = 0; k < nM; ++k) {
        block_sum_u[k].fill(0.0);
        block_sum_u2[k].fill(0.0);
        block_sum_uu[k].fill(0.0);
        block_first[k].fill(0.0);
        block_last[k].fill(0.0);
    }

    stream_gaps(LO, HI, [&](uint64_t p, uint64_t g) {
        double lp = (double)log((double)p);
        double e = (double)g - (alpha + beta*lp);
        n2++;
        sum_e += e; sum_e2 += e*e;
        if (n2 > 1) sum_ee += prev_e * e;
        prev_e = e;

        for (size_t k = 0; k < nM; ++k) {
            uint32_t r = (uint32_t)(p % Ms[k]);
            double mu_cur = mu[k][r];
            double res_cur = e - mu_cur;

            acc[k].sum_res2 += res_cur * res_cur;
            acc[k].cnt_res++;

            if (acc[k].has_prev) {
                acc[k].A += acc[k].prev_mu * mu_cur;
                acc[k].B += acc[k].prev_mu * res_cur;
                acc[k].C += acc[k].prev_res * mu_cur;
                acc[k].D += acc[k].prev_res * res_cur;
                acc[k].sum_resres += acc[k].prev_res * res_cur;
                acc[k].cnt_ABCD++;
            }

            // ---- block-level accumulators (5 per block, no residuals kept) ----
            size_t global_idx = n2 - 1;                 // current gap index
            size_t bidx = global_idx / block_size;
            if (bidx >= (size_t)N_BLOCKS) bidx = N_BLOCKS - 1;   // last block absorbs the remainder

            block_sum_u[k][bidx]  += res_cur;
            block_sum_u2[k][bidx] += res_cur * res_cur;
            if (acc[k].has_prev && ((global_idx - 1) / block_size) == bidx)
                block_sum_uu[k][bidx] += acc[k].prev_res * res_cur;  // within-block pairs only
            if (global_idx == bidx * block_size)
                block_first[k][bidx] = res_cur;
            block_last[k][bidx] = res_cur;   // blocks are contiguous: last write = block's last gap

            acc[k].prev_mu = mu_cur;
            acc[k].prev_res = res_cur;
            acc[k].has_prev = true;
        }
    });

    assert(n1 == n2);

    double mean_e = (double)(sum_e / dn);
    double var_e  = (double)(sum_e2 / dn) - mean_e*mean_e;
    cerr << "  Var(e) measured         = " << var_e << "\n";
    cerr << "  closure diff (should be < 1e-4) = " << var_e - var_e_pred << "\n";
    double cov_1  = (double)(sum_ee / (long double)(n2 - 1)) - mean_e*mean_e;
    double sd_g   = sqrt(var_e);
    double logp_center = 0.5 * (log((double)LO) + log((double)HI));
    double d_meas = logp_center - sd_g;

    // ---- variance split: var_e = var_between + var_within (law of total variance)
    vector<double> var_btw(nM), var_wth(nM);
    for (size_t k = 0; k < nM; ++k) {
        double vb = 0;
        for (uint64_t r = 0; r < Ms[k]; ++r) {
            if (cnt[k][r] == 0) continue;
            double w   = (double)cnt[k][r] / (double)n2;
            double dmu = mu[k][r] - mean_e;
            vb += w * dmu * dmu;
        }
        var_btw[k] = vb;
        var_wth[k] = var_e - vb;
    }

    // ---- Block bootstrap for rho_resid -----------------------------------
    // 40 blocks, 300 replicates, independent seed per M.
    // rho_rep = sum_uu_new / sum_u2_new (residuals are mean-zero, so no
    // centering is needed; cross-block terms are the 39 last*first pairs).
    vector<double> boot_lo(nM, 0.0), boot_med(nM, 0.0), boot_hi(nM, 0.0);
    for (size_t k = 0; k < nM; ++k) {
        mt19937_64 rng(SEED_BASE + k * 1000);        // independent seed per M
        uniform_int_distribution<int> pick(0, N_BLOCKS - 1);

        vector<double> rho_b;
        rho_b.reserve(N_BOOT);
        for (int b = 0; b < N_BOOT; ++b) {
            double su2 = 0, suu = 0, prev_last = 0;
            for (int j = 0; j < N_BLOCKS; ++j) {
                int bi = pick(rng);
                su2 += block_sum_u2[k][bi];
                suu += block_sum_uu[k][bi];
                if (j > 0) suu += prev_last * block_first[k][bi];
                prev_last = block_last[k][bi];
            }
            if (su2 > 0) rho_b.push_back(suu / su2);
        }
        if (rho_b.empty()) throw runtime_error("bootstrap produced no samples");
        sort(rho_b.begin(), rho_b.end());
        boot_lo[k]  = rho_b[(size_t)(0.025 * rho_b.size())];
        boot_med[k] = rho_b[rho_b.size() / 2];
        boot_hi[k]  = rho_b[(size_t)(0.975 * rho_b.size())];
        cerr << "  M = " << Ms[k] << ": bootstrap CI95 = ["
             << boot_lo[k] << ", " << boot_med[k] << ", " << boot_hi[k] << "]\n";
    }

    {
        ostringstream ss;
        ss << fixed << setprecision(6);
        ss << "\n============================================================\n";
        ss << "  [" << LO << ", " << HI << ")  RESULTS\n";
        ss << "============================================================\n";
        ss << "  n gaps (n1=n2)    = " << n2 << "\n";
        ss << "  first prime       = " << first_p << "\n";
        ss << "  last opening      = " << last_open << "\n";
        ss << "  min, max gap      = " << min_gap << ", " << max_gap << "\n";
        ss << "  log p center      = " << logp_center << "\n";
        ss << "  alpha, beta       = " << alpha << ", " << beta << "\n";
        ss << "  Var(g)            = " << var_e << "   (pred 438-445)\n";
        ss << "  Var(g) raw        = " << var_g << "\n";
        ss << "  Var(log p)        = " << var_lp << "\n";
        ss << "  beta^2 Var(log p) = " << beta_check*beta_check*var_lp << "\n";
        ss << "  sd(g)             = " << sd_g << "   (pred 20.9-21.1)\n";
        ss << "  d = L - sd(g)     = " << d_meas << "   (pred 2.5-2.7)\n";
        ss << "  Cov(g_i,g_i+1)    = " << cov_1 << "   (pred -11.0 to -11.8)\n";
        ss << "  |Cov|/(L/2)       = " << fabs(cov_1)/(logp_center/2)
           << "   (pred 0.93-1.00)\n";
        ss << "  rho_1             = " << (cov_1/var_e)
           << "   (pred -0.0248 to -0.0269)\n\n";
        P(ss.str());
    }

    for (size_t k = 0; k < nM; ++k) {
        double A = (double)(acc[k].A / acc[k].cnt_ABCD);
        double B = (double)(acc[k].B / acc[k].cnt_ABCD);
        double C = (double)(acc[k].C / acc[k].cnt_ABCD);
        double D = (double)(acc[k].D / acc[k].cnt_ABCD);
        double sum_ABCD = A + B + C + D;
        double var_res = (double)(acc[k].sum_res2 / acc[k].cnt_res);
        double cov_res = (double)(acc[k].sum_resres / (acc[k].cnt_res - 1));
        double rho_res_raw = cov_res / var_res;
        double n_c = (double)n2 / phiM[k];
        double corr = (n_c > 1.5) ? (1.0 - 1.0/n_c) : 1.0;
        double rho_res_corr = rho_res_raw / corr;

        ostringstream ss;
        ss << fixed << setprecision(6);
        ss << "  M = " << Ms[k] << "  phi = " << phiM[k]
           << "  n_c = " << n_c << "\n";
        ss << "    A, B, C, D       = " << A << ", " << B << ", " << C
           << ", " << D << "\n";
        ss << "    A+B+C+D          = " << sum_ABCD << "\n";
        ss << "    Cov              = " << cov_1 << "\n";
        ss << "    (ABCD - Cov)     = " << sum_ABCD - cov_1 << "\n";
        ss << "    (ABCD / Cov)     = " << sum_ABCD / cov_1 << "\n";
        ss << "    pairs_ABCD       = " << acc[k].cnt_ABCD << "\n";
        ss << "    rho_res (raw)    = " << rho_res_raw << "\n";
        ss << "    rho_res (corr)   = " << rho_res_corr << "\n";
        ss << "    var_between      = " << var_btw[k] << "\n";
        ss << "    var_within       = " << var_wth[k] << "\n";
        ss << "    boot95 lo/med/hi = " << boot_lo[k] << ", " << boot_med[k]
           << ", " << boot_hi[k] << "\n\n";
        P(ss.str());
    }

    double dt = chrono::duration<double>(chrono::steady_clock::now() - T0).count();
    {
        ostringstream ss;
        ss << "  total time = " << (int)dt << " s\n";
        P(ss.str());
    }

    // CSV output
    // NOTE: cnt[][] (Pass 1 class counts) is retained here (only sum_g/sum_lp
    //       were cleared), so min_class is computed directly from it.
    ofstream csv("phase_decomposition_v2.csv");
    csv << "window,M,n_gaps,logp_center,n_c,min_class,var_g,"
        << "cov_total,A,X,Y,D,ABCD_minus_cov,ABCD_over_cov,"
        << "rho_resid_raw,rho_resid_corrected,"
        << "var_between,var_within,boot_lo,boot_med,boot_hi,"
        << "var_g_raw,var_logp\n";

    for (size_t k = 0; k < nM; ++k) {
        double A = (double)(acc[k].A / acc[k].cnt_ABCD);
        double X = (double)(acc[k].B / acc[k].cnt_ABCD);  // rename B→X
        double Y = (double)(acc[k].C / acc[k].cnt_ABCD);  // rename C→Y
        double D = (double)(acc[k].D / acc[k].cnt_ABCD);
        double sum_ABCD = A + X + Y + D;

        double var_res = (double)(acc[k].sum_res2 / acc[k].cnt_res);
        double cov_res = (double)(acc[k].sum_resres / (acc[k].cnt_res - 1));
        double rho_res_raw = cov_res / var_res;
        double n_c = (double)n2 / phiM[k];
        double corr = (n_c > 1.5) ? (1.0 - 1.0/n_c) : 1.0;
        double rho_res_corr = rho_res_raw / corr;

        // min class calculation (skip empty classes)
        uint32_t min_class = UINT32_MAX;
        for (uint64_t r = 0; r < Ms[k]; ++r) {
            if (cnt[k][r] > 0 && cnt[k][r] < min_class)
                min_class = cnt[k][r];
        }

        // window field is quoted because it contains a comma (valid CSV)
        csv << "\"[" << LO << "," << HI << ")\""
            << "," << Ms[k] << ","
            << n2 << ","
            << fixed << setprecision(6) << logp_center << ","
            << setprecision(2) << n_c << ","
            << min_class << ","
            << setprecision(6) << var_e << ","
            << cov_1 << ","
            << A << "," << X << "," << Y << "," << D << ","
            << (sum_ABCD - cov_1) << ","
            << (sum_ABCD / cov_1) << ","
            << rho_res_raw << ","
            << rho_res_corr << ","
            << var_btw[k] << ","
            << var_wth[k] << ","
            << boot_lo[k] << ","
            << boot_med[k] << ","
            << boot_hi[k] << "," << var_g << "," << var_lp << "\n";
    }
    csv.close();

    fout.close();
    cout << "\n  [Saved to test_1e10_final_v3.txt]\n";
    return 0;
}
