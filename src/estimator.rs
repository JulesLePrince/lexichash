/// `1 / ln(4)`, independent of `k` and `n`.
const LN4_INV: f64 = 1.0 / (2.0 * std::f64::consts::LN_2);

/// `mu`/`sigma2` fit for `k` in `[19, 41]`, `n` in `[1e6, 1e9]`; refit if
/// that range changes. `mu`'s correction saturates with `k - log4(n)`
/// rather than being constant.
const MU_INF: f64 = 0.265_89;
const MU_A: f64 = 0.133_89;
const MU_TAU: f64 = 5.559_8;
const SIGMA2_SLOPE: f64 = 0.074_900;
const SIGMA2_INTERCEPT: f64 = 0.260_58;

/// `Sigma_Phi(n) ~= log4(n) + KAPPA` (Mellin/Frullani asymptotic, matches
/// the exact recurrence to ~1e-4). `KAPPA = gamma/ln(4) - 1 + (4/9)(log2(5) - 2)`.
const KAPPA: f64 = -0.4405477580782944;

/// Estimates the mutation rate from a LexicHash mean score.
#[derive(Debug, Clone)]
pub struct MutationRateEstimator {
    k: usize,
    m: f64,
    sigma_phi: f64,
    e0: f64,
    e0_inv: f64,
    mu: f64,
    sigma2: f64,
    b1: f64,
    b2: f64,
    b2_inv: f64,
}

impl MutationRateEstimator {
    /// `k` is the sketch size, `n` the number of k-mers.
    #[inline(always)]
    pub fn new(k: usize, num_kmers: usize) -> Self {
        let mut res = unsafe { Self::new_uninit() };
        res.rebuild(k, num_kmers);
        res
    }

    #[allow(clippy::missing_safety_doc)]
    #[inline(always)]
    pub const unsafe fn new_uninit() -> Self {
        Self {
            k: 0,
            m: 0.0,
            sigma_phi: 0.0,
            e0: 0.0,
            e0_inv: 0.0,
            mu: 0.0,
            sigma2: 0.0,
            b1: 0.0,
            b2: 0.0,
            b2_inv: 0.0,
        }
    }

    /// Rebuilds this estimator for a new `(k, n)`, in place, using
    /// closed-form asymptotic approximations.
    #[inline(always)]
    pub fn rebuild(&mut self, k: usize, num_kmers: usize) {
        assert!(k >= 1, "k must be at least 1");

        let kf = k as f64;
        let n = num_kmers as f64;
        let log4_n = n.ln() * LN4_INV;
        let m = log4_n - 0.8;

        let sigma_phi = log4_n + KAPPA;
        let e0 = kf - sigma_phi;

        let width = kf - log4_n;
        let mu = (log4_n + kf) / 2.0 + MU_INF - MU_A * (-width / MU_TAU).exp();
        let sigma2 = width * width / 12.0 + SIGMA2_SLOPE * width + SIGMA2_INTERCEPT;

        self.set_moments(k, m, sigma_phi, e0, mu, sigma2);
    }

    /// Exact `O(k)`-loop counterpart to `rebuild`.
    pub fn rebuild_precise(&mut self, k: usize, num_kmers: usize) {
        assert!(k >= 1, "k must be at least 1");

        let n = num_kmers as f64;
        let m = n.ln() * LN4_INV - 0.8;

        let mut p = 0.25_f64; // 4^-1, 4^-2, ... as we go
        let mut c_prev = 1.0 - (-n * p).exp();
        let mut s = 0.0_f64;
        let mut sigma_phi = 0.0_f64;
        let mut e0 = 0.0_f64;
        let mut sum_j_e = 0.0_f64;
        let mut sum_j2_e = 0.0_f64;

        for j in 1..=k {
            let jf = j as f64;

            let phi_j = c_prev * c_prev + s * (1.0 / 3.0);
            let e_j = 1.0 - phi_j;

            sigma_phi += phi_j;
            e0 += e_j;
            sum_j_e += jf * e_j;
            sum_j2_e += jf * jf * e_j;

            if j < k {
                p *= 0.25;
                let c_next = 1.0 - (-n * p).exp();
                let d = c_prev - c_next;
                s = d * d + s / 4.0;
                c_prev = c_next;
            }
        }

        let mu = sum_j_e / e0;
        let sigma2 = sum_j2_e / e0 - mu * mu;

        self.set_moments(k, m, sigma_phi, e0, mu, sigma2);
    }

    #[inline(always)]
    fn set_moments(&mut self, k: usize, m: f64, sigma_phi: f64, e0: f64, mu: f64, sigma2: f64) {
        let b1 = mu + m;
        let b2 = sigma2 + m * m;

        self.k = k;
        self.m = m;
        self.sigma_phi = sigma_phi;
        self.e0 = e0;
        self.e0_inv = 1.0 / e0;
        self.mu = mu;
        self.sigma2 = sigma2;
        self.b1 = b1;
        self.b2 = b2;
        self.b2_inv = 1.0 / b2;
    }

    /// Estimate the mutation rate from an observed mean score.
    ///
    /// `NEWTON_STEPS` corresponds to the number of refinement steps, 2 is usually enough.
    /// Returns `None` if `score` is out of the reachable range for this `(n, k)`.
    #[inline(always)]
    pub fn estimate_mut_rate<const NEWTON_STEPS: usize>(&self, score: f64) -> Option<f64> {
        if !(self.sigma_phi..=self.k as f64).contains(&score) {
            return None;
        }

        // Closed-form seed.
        let r = ((score - self.sigma_phi) * self.e0_inv).clamp(f64::MIN_POSITIVE, 1.0);
        let log_r = r.ln();
        let disc = (self.b1 * self.b1 + 2.0 * self.b2 * log_r).max(0.0);
        let mut rho = ((self.b1 - disc.sqrt()) * self.b2_inv).max(0.0);

        // Newton refinement.
        for _ in 0..NEWTON_STEPS {
            let w = 1.0 / (1.0 + self.m * rho);
            let t = self.e0 * (-self.mu * rho + 0.5 * self.sigma2 * rho * rho).exp();
            let f = self.sigma_phi + w * t - score;
            let fp = w * t * (self.sigma2 * rho - self.mu - self.m * w);
            rho = (rho - f / fp).max(0.0);
        }

        Some(rho)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values from the exact O(k^2) Python implementation,
    // n = 1_000_000, k = 31.
    const N: usize = 1_000_000;
    const K: usize = 31;
    const CASES: [(f64, f64); 7] = [
        (0.005, 28.0447),
        (0.010, 25.5424),
        (0.020, 21.6031),
        (0.030, 18.7229),
        (0.050, 15.0009),
        (0.070, 12.8864),
        (0.100, 11.2222),
    ];

    #[test]
    fn seed_only_is_in_the_right_ballpark() {
        let est = MutationRateEstimator::new(K, N);
        for (rho, score) in CASES {
            let recovered = est.estimate_mut_rate::<0>(score).unwrap();
            assert!(
                (recovered - rho).abs() / rho < 0.10,
                "rho={rho}: recovered={recovered}"
            );
        }
    }

    #[test]
    fn two_newton_steps_match_python_within_half_a_percent() {
        let est = MutationRateEstimator::new(K, N);
        for (rho, score) in CASES {
            let recovered = est.estimate_mut_rate::<2>(score).unwrap();
            let rel_err = (recovered - rho).abs() / rho;
            assert!(
                rel_err < 0.005,
                "rho={rho}: recovered={recovered}, rel_err={rel_err}"
            );
        }
    }

    #[test]
    fn newton_iteration_converges_and_does_not_diverge() {
        let est = MutationRateEstimator::new(K, N);
        for (rho, score) in CASES {
            let r2 = est.estimate_mut_rate::<2>(score).unwrap();
            let r5 = est.estimate_mut_rate::<5>(score).unwrap();
            let r8 = est.estimate_mut_rate::<8>(score).unwrap();
            assert!((r2 - r5).abs() / rho < 0.005, "rho={rho}: r2={r2}, r5={r5}");
            assert!(
                (r5 - r8).abs() < 1e-9,
                "rho={rho}: r5={r5}, r8={r8} (should have converged)"
            );
        }
    }

    #[test]
    fn out_of_range_score_returns_none() {
        let est = MutationRateEstimator::new(K, N);
        assert!(est.estimate_mut_rate::<2>(K as f64 + 1.0).is_none());
        assert!(est.estimate_mut_rate::<2>(0.0).is_none());
    }

    #[test]
    fn build_matches_fresh_construction() {
        let mut est = MutationRateEstimator::new(K, N / 2);
        est.rebuild(K, N);
        let fresh = MutationRateEstimator::new(K, N);
        for (_, score) in CASES {
            assert_eq!(
                est.estimate_mut_rate::<2>(score),
                fresh.estimate_mut_rate::<2>(score)
            );
        }
    }

    #[test]
    fn rebuild_precise_matches_python_within_half_a_percent() {
        let mut est = MutationRateEstimator::new(K, N);
        est.rebuild_precise(K, N);
        for (rho, score) in CASES {
            let recovered = est.estimate_mut_rate::<2>(score).unwrap();
            let rel_err = (recovered - rho).abs() / rho;
            assert!(
                rel_err < 0.005,
                "rho={rho}: recovered={recovered}, rel_err={rel_err}"
            );
        }
    }

    #[test]
    fn rebuild_matches_rebuild_precise() {
        // rebuild() should track rebuild_precise() closely.
        let mut precise_est = MutationRateEstimator::new(K, N);
        precise_est.rebuild_precise(K, N);
        let fast_est = MutationRateEstimator::new(K, N);
        for (rho, score) in CASES {
            let precise = precise_est.estimate_mut_rate::<2>(score).unwrap();
            let fast = fast_est.estimate_mut_rate::<2>(score).unwrap();
            let rel_err = (precise - fast).abs() / rho;
            assert!(
                rel_err < 0.005,
                "rho={rho}: precise={precise}, fast={fast}, rel_err={rel_err}"
            );
        }
    }

    #[test]
    #[ignore = "This is a benchmark, not a test"]
    fn bench_estimate_mut_rate() {
        use core::hint::black_box;

        let est = MutationRateEstimator::new(K, N);
        let rep = 20_000_000;

        let start = std::time::Instant::now();
        let mut acc = 0.0;
        for _ in 0..rep {
            acc += est.estimate_mut_rate::<2>(black_box(15.0)).unwrap();
        }
        black_box(acc);
        eprintln!(
            "estimate_mut_rate: {:.02} ns/call",
            start.elapsed().as_secs_f64() * 1e9 / rep as f64
        );
    }

    #[test]
    #[ignore = "This is a benchmark, not a test"]
    fn bench_rebuild() {
        use core::hint::black_box;

        let mut est = MutationRateEstimator::new(K, N);
        let rep = 20_000_000;

        let start = std::time::Instant::now();
        for _ in 0..rep {
            est.rebuild(K, black_box(N));
        }
        black_box(&est);
        eprintln!(
            "rebuild: {:.02} ns/call",
            start.elapsed().as_secs_f64() * 1e9 / rep as f64
        );
    }
}
