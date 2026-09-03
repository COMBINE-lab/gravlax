//! Statistics for 3'-site analysis: internal-priming detection against the genome, the
//! per-gene site×group G-test, permutation p-values, and BH-FDR.
//!
//! Internal priming: oligo(dT) also primes on genomic A-stretches inside transcripts, creating
//! 3' "sites" that are sequence artifacts, not cleavage sites. Field-standard filters look at the
//! genomic sequence just downstream of the putative cleavage position in transcript orientation
//! (an A-rich window or a long A-run means the "tail" may have been templated). Crucially this
//! artifact replicates across datasets — it is sequence-driven — so replication cannot rebut it;
//! only the genome can. Defaults follow the literature: ≥12 A in the 20 nt immediately
//! downstream, or a run of ≥8 A within 140 nt downstream.

/// A-content just downstream of a cleavage position, in transcript orientation.
/// `tp` is the molecule 3' coordinate as the archive defines it: on + the 0-based *exclusive*
/// end (so downstream starts at index `tp`), on − the 0-based *inclusive* start (downstream is
/// leftward; transcript A = genomic T there). Returns (a_in_20, longest_a_run_in_140).
pub fn ip_stats(seq: &[u8], tp: u32, rev: bool) -> (u32, u32) {
    let n = seq.len();
    let (window, base): (&[u8], u8) = if !rev {
        let s = (tp as usize).min(n);
        (&seq[s..(s + 140).min(n)], b'A')
    } else {
        let e = (tp as usize).min(n);
        (&seq[e.saturating_sub(140)..e], b'T')
    };
    // Orient so index 0 is the base immediately downstream in transcript direction.
    let it: Box<dyn Iterator<Item = u8> + '_> = if !rev {
        Box::new(window.iter().copied())
    } else {
        Box::new(window.iter().rev().copied())
    };
    let (mut a20, mut run, mut best) = (0u32, 0u32, 0u32);
    for (i, b) in it.enumerate() {
        let is_a = b == base;
        if is_a && i < 20 {
            a20 += 1;
        }
        run = if is_a { run + 1 } else { 0 };
        best = best.max(run);
    }
    (a20, best)
}

pub const IP_A20: u32 = 12;
pub const IP_RUN: u32 = 8;

pub fn is_internal_priming(a20: u32, run: u32) -> bool {
    a20 >= IP_A20 || run >= IP_RUN
}

/// Likelihood-ratio (G) statistic for an S×G contingency table of UMI counts.
/// Returns (G, df) over rows/columns with nonzero margins.
pub fn g_statistic(table: &[Vec<u64>]) -> (f64, u64) {
    let ngroups = table.first().map(|r| r.len()).unwrap_or(0);
    let row_sums: Vec<u64> = table.iter().map(|r| r.iter().sum()).collect();
    let col_sums: Vec<u64> = (0..ngroups).map(|j| table.iter().map(|r| r[j]).sum()).collect();
    let total: u64 = row_sums.iter().sum();
    if total == 0 {
        return (0.0, 0);
    }
    let mut g = 0.0f64;
    for (i, row) in table.iter().enumerate() {
        for (j, &o) in row.iter().enumerate() {
            if o > 0 && row_sums[i] > 0 && col_sums[j] > 0 {
                let e = row_sums[i] as f64 * col_sums[j] as f64 / total as f64;
                g += 2.0 * o as f64 * (o as f64 / e).ln();
            }
        }
    }
    let live_rows = row_sums.iter().filter(|&&s| s > 0).count() as u64;
    let live_cols = col_sums.iter().filter(|&&s| s > 0).count() as u64;
    let df = live_rows.saturating_sub(1) * live_cols.saturating_sub(1);
    (g.max(0.0), df)
}

/// Upper tail of the chi-square distribution: P(X ≥ x) with k degrees of freedom, via the
/// regularized incomplete gamma function (series for x < k+1, continued fraction otherwise).
pub fn chi2_sf(x: f64, k: u64) -> f64 {
    if k == 0 || x <= 0.0 {
        return 1.0;
    }
    let (a, x) = (k as f64 / 2.0, x / 2.0);
    // ln Γ(a), Lanczos.
    let gammaln = |z: f64| -> f64 {
        const C: [f64; 6] = [
            76.18009172947146,
            -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
            0.1208650973866179e-2,
            -0.5395239384953e-5,
        ];
        let mut ser = 1.000000000190015;
        for (j, c) in C.iter().enumerate() {
            ser += c / (z + 1.0 + j as f64);
        }
        let tmp = z + 5.5 - (z + 0.5) * (z + 5.5).ln();
        -tmp + (2.5066282746310005 * ser / z).ln()
    };
    if x < a + 1.0 {
        // P(a,x) by series; return 1 - P.
        let mut ap = a;
        let mut sum = 1.0 / a;
        let mut del = sum;
        for _ in 0..500 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-15 {
                break;
            }
        }
        1.0 - sum * (-x + a * x.ln() - gammaln(a)).exp()
    } else {
        // Q(a,x) by continued fraction (Lentz).
        let (mut b, mut c, mut d) = (x + 1.0 - a, 1e300f64, 1.0 / (x + 1.0 - a));
        let mut h = d;
        for i in 1..500 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = (an * d + b).recip_or(1e-300);
            c = b + an / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-15 {
                break;
            }
        }
        (h * (-x + a * x.ln() - gammaln(a)).exp()).clamp(0.0, 1.0)
    }
}

trait RecipOr {
    fn recip_or(self, floor: f64) -> f64;
}
impl RecipOr for f64 {
    fn recip_or(self, floor: f64) -> f64 {
        if self.abs() < floor {
            floor.recip()
        } else {
            self.recip()
        }
    }
}

/// Benjamini–Hochberg adjusted q-values (monotone), same order as input.
pub fn bh_fdr(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| pvals[a].partial_cmp(&pvals[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut q = vec![0.0f64; n];
    let mut running = 1.0f64;
    for (rank, &i) in idx.iter().enumerate().rev() {
        let val = (pvals[i] * n as f64 / (rank + 1) as f64).min(1.0);
        running = running.min(val);
        q[i] = running;
    }
    q
}

/// Two-sided paired Student t test and exact sign-flip calibration for paired differences.
/// The exact p-value enumerates every sign assignment, so callers cannot accidentally treat
/// cells or molecules as replicates when the supplied vector contains biological samples.
pub fn paired_test(differences: &[f64]) -> Option<(f64, f64, f64)> {
    if differences.len() < 2 || differences.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let n = differences.len() as f64;
    let mean = differences.iter().sum::<f64>() / n;
    let variance = differences
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    let t = if variance == 0.0 {
        if mean == 0.0 {
            0.0
        } else {
            mean.signum() * f64::INFINITY
        }
    } else {
        mean / (variance / n).sqrt()
    };
    let p_t = if t.is_infinite() {
        0.0
    } else if t == 0.0 {
        1.0
    } else {
        // For t with nu degrees of freedom, the two-sided tail is
        // I_{nu/(nu+t^2)}(nu/2, 1/2).
        let nu = n - 1.0;
        regularized_beta(nu / (nu + t * t), nu / 2.0, 0.5)
    };
    let p_flip = if differences.len() <= 20 {
        let observed = mean.abs();
        let assignments = 1u64 << differences.len();
        let mut extreme = 0u64;
        for mask in 0..assignments {
            let signed = differences
                .iter()
                .enumerate()
                .map(|(index, value)| {
                    if mask & (1u64 << index) == 0 {
                        *value
                    } else {
                        -*value
                    }
                })
                .sum::<f64>()
                / n;
            if signed.abs() + 1e-15 >= observed {
                extreme += 1;
            }
        }
        extreme as f64 / assignments as f64
    } else {
        f64::NAN
    };
    Some((t, p_t.clamp(0.0, 1.0), p_flip.clamp(0.0, 1.0)))
}

fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_fraction(x, a, b) / a
    } else {
        1.0 - front * beta_fraction(1.0 - x, b, a) / b
    }
}

fn beta_fraction(x: f64, a: f64, b: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=500 {
        let m2 = 2 * m;
        let mut aa = m as f64 * (b - m as f64) * x / ((qam + m2 as f64) * (a + m2 as f64));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m as f64) * (qab + m as f64) * x / ((a + m2 as f64) * (qap + m2 as f64));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < 3e-14 {
            break;
        }
    }
    h
}

/// Deterministic maximum-likelihood beta-binomial contrast for one path versus the remaining
/// same-strand path fragments. Each tuple is `(path_umis, strand_path_umis)` for one biological
/// sample. The two models each fit one shared concentration; the alternative has separate group
/// means and the null has one common mean.
#[derive(Clone, Copy, Debug)]
pub struct BetaBinomialContrast {
    pub null_mean: f64,
    pub null_concentration: f64,
    pub condition_a_mean: f64,
    pub condition_b_mean: f64,
    pub alternative_concentration: f64,
    pub null_log_likelihood: f64,
    pub alternative_log_likelihood: f64,
    pub likelihood_ratio: f64,
    pub p_value: f64,
}

const BB_MEAN_EPSILON: f64 = 1e-8;
const BB_MIN_CONCENTRATION: f64 = 1e-3;
const BB_MAX_CONCENTRATION: f64 = 1e6;

fn log_gamma(z: f64) -> f64 {
    // Numerical Recipes' six-coefficient Lanczos approximation. Every caller supplies z > 0.
    const COEFFICIENTS: [f64; 6] = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ];
    let mut series = 1.000000000190015;
    for (index, coefficient) in COEFFICIENTS.iter().enumerate() {
        series += coefficient / (z + 1.0 + index as f64);
    }
    let tmp = z + 5.5 - (z + 0.5) * (z + 5.5).ln();
    -tmp + (2.5066282746310005 * series / z).ln()
}

fn log_beta(a: f64, b: f64) -> f64 {
    log_gamma(a) + log_gamma(b) - log_gamma(a + b)
}

fn beta_binomial_log_likelihood(
    observations: &[(u64, u64)],
    mean: f64,
    concentration: f64,
) -> f64 {
    let mean = mean.clamp(BB_MEAN_EPSILON, 1.0 - BB_MEAN_EPSILON);
    let concentration = concentration.clamp(BB_MIN_CONCENTRATION, BB_MAX_CONCENTRATION);
    let alpha = mean * concentration;
    let beta = (1.0 - mean) * concentration;
    observations
        .iter()
        .map(|&(successes, total)| {
            log_gamma(total as f64 + 1.0)
                - log_gamma(successes as f64 + 1.0)
                - log_gamma((total - successes) as f64 + 1.0)
                + log_beta(
                    successes as f64 + alpha,
                    (total - successes) as f64 + beta,
                )
                - log_beta(alpha, beta)
        })
        .sum()
}

fn golden_max<F>(mut lower: f64, mut upper: f64, iterations: usize, objective: F) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    const RATIO: f64 = 0.6180339887498949;
    let mut left = upper - RATIO * (upper - lower);
    let mut right = lower + RATIO * (upper - lower);
    let mut left_value = objective(left);
    let mut right_value = objective(right);
    for _ in 0..iterations {
        if left_value < right_value {
            lower = left;
            left = right;
            left_value = right_value;
            right = lower + RATIO * (upper - lower);
            right_value = objective(right);
        } else {
            upper = right;
            right = left;
            right_value = left_value;
            left = upper - RATIO * (upper - lower);
            left_value = objective(left);
        }
    }
    let candidates = [
        (lower, objective(lower)),
        (left, left_value),
        (right, right_value),
        (upper, objective(upper)),
    ];
    candidates
        .into_iter()
        .max_by(|left, right| {
            left.1
                .partial_cmp(&right.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap()
}

fn fit_mean(observations: &[(u64, u64)], concentration: f64) -> (f64, f64) {
    let successes: u64 = observations.iter().map(|row| row.0).sum();
    let total: u64 = observations.iter().map(|row| row.1).sum();
    if successes == 0 {
        let mean = BB_MEAN_EPSILON;
        return (
            mean,
            beta_binomial_log_likelihood(observations, mean, concentration),
        );
    }
    if successes == total {
        let mean = 1.0 - BB_MEAN_EPSILON;
        return (
            mean,
            beta_binomial_log_likelihood(observations, mean, concentration),
        );
    }
    golden_max(
        BB_MEAN_EPSILON,
        1.0 - BB_MEAN_EPSILON,
        64,
        |mean| beta_binomial_log_likelihood(observations, mean, concentration),
    )
}

fn fit_profile_concentration<F>(profile: F) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    const GRID_POINTS: usize = 49;
    let lower = BB_MIN_CONCENTRATION.ln();
    let upper = BB_MAX_CONCENTRATION.ln();
    let step = (upper - lower) / (GRID_POINTS - 1) as f64;
    let grid: Vec<(f64, f64)> = (0..GRID_POINTS)
        .map(|index| {
            let log_concentration = lower + index as f64 * step;
            (log_concentration, profile(log_concentration.exp()))
        })
        .collect();
    let best = grid
        .iter()
        .enumerate()
        .max_by(|(_, left), (_, right)| {
            left.1
                .partial_cmp(&right.1)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(index, _)| index)
        .unwrap();
    let bracket_lower = grid[best.saturating_sub(1)].0;
    let bracket_upper = grid[(best + 1).min(GRID_POINTS - 1)].0;
    let (log_concentration, likelihood) = if bracket_lower == bracket_upper {
        grid[best]
    } else {
        golden_max(bracket_lower, bracket_upper, 48, |log_concentration| {
            profile(log_concentration.exp())
        })
    };
    (log_concentration.exp(), likelihood)
}

pub fn beta_binomial_contrast(
    condition_a: &[(u64, u64)],
    condition_b: &[(u64, u64)],
) -> Option<BetaBinomialContrast> {
    if condition_a.is_empty()
        || condition_b.is_empty()
        || condition_a
            .iter()
            .chain(condition_b)
            .any(|&(successes, total)| total == 0 || successes > total)
    {
        return None;
    }
    let combined: Vec<(u64, u64)> = condition_a
        .iter()
        .chain(condition_b)
        .copied()
        .collect();
    let (null_concentration, null_log_likelihood) =
        fit_profile_concentration(|concentration| fit_mean(&combined, concentration).1);
    let (null_mean, _) = fit_mean(&combined, null_concentration);
    let (alternative_concentration, alternative_log_likelihood) =
        fit_profile_concentration(|concentration| {
            fit_mean(condition_a, concentration).1 + fit_mean(condition_b, concentration).1
        });
    let (condition_a_mean, _) = fit_mean(condition_a, alternative_concentration);
    let (condition_b_mean, _) = fit_mean(condition_b, alternative_concentration);
    let likelihood_ratio =
        (2.0 * (alternative_log_likelihood - null_log_likelihood)).max(0.0);
    Some(BetaBinomialContrast {
        null_mean,
        null_concentration,
        condition_a_mean,
        condition_b_mean,
        alternative_concentration,
        null_log_likelihood,
        alternative_log_likelihood,
        likelihood_ratio,
        p_value: chi2_sf(likelihood_ratio, 1),
    })
}

/// Deterministic xorshift64* for permutation tests — reproducible runs, no rand dependency.
pub struct Rng(u64);
impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed.max(1))
    }
    pub fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
    }
}

/// Permutation p-value for the observed G: shuffle cell→group labels, rebuild the site×group
/// table from per-cell per-site counts, and count permuted G ≥ observed.
/// `cell_site_counts[c]` lists (site, umis) for grouped cell c; `groups[c]` is its group.
pub fn permutation_p(
    cell_site_counts: &[Vec<(usize, u64)>],
    groups: &[usize],
    n_sites: usize,
    n_groups: usize,
    g_obs: f64,
    n_perm: usize,
    seed: u64,
) -> f64 {
    let mut rng = Rng::new(seed);
    let mut labels: Vec<usize> = groups.to_vec();
    let mut ge = 0usize;
    for _ in 0..n_perm {
        rng.shuffle(&mut labels);
        let mut table = vec![vec![0u64; n_groups]; n_sites];
        for (c, sites) in cell_site_counts.iter().enumerate() {
            for &(s, n) in sites {
                table[s][labels[c]] += n;
            }
        }
        let (g, _) = g_statistic(&table);
        if g >= g_obs {
            ge += 1;
        }
    }
    (1 + ge) as f64 / (1 + n_perm) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chi2_sf_matches_known_values() {
        // 1 df: P(X ≥ 3.841) ≈ 0.05; 2 df: P(X ≥ 5.991) ≈ 0.05; 5 df: P(X ≥ 11.070) ≈ 0.05.
        assert!((chi2_sf(3.841, 1) - 0.05).abs() < 1e-3);
        assert!((chi2_sf(5.991, 2) - 0.05).abs() < 1e-3);
        assert!((chi2_sf(11.070, 5) - 0.05).abs() < 1e-3);
        assert!((chi2_sf(0.0, 3) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn g_test_null_and_signal() {
        // Perfectly proportional table: G = 0.
        let (g, df) = g_statistic(&[vec![10, 20], vec![30, 60]]);
        assert!(g < 1e-9);
        assert_eq!(df, 1);
        // Strong redistribution: large G, small p.
        let (g, df) = g_statistic(&[vec![90, 10], vec![10, 90]]);
        assert_eq!(df, 1);
        assert!(chi2_sf(g, df) < 1e-10);
    }

    #[test]
    fn ip_stats_orientation() {
        //                  0123456789
        let seq = b"CCCCCTTTTTAAAAAAAAAAAACGT".to_vec(); // A-run of 12 starting at index 10
        let (a20, run) = ip_stats(&seq, 10, false);
        assert_eq!(run, 12);
        assert!(a20 >= 12);
        assert!(is_internal_priming(a20, run));
        // Reverse strand at tp=10: upstream (leftward) is CCCCCTTTTT; transcript-A = genomic T.
        let (a20r, runr) = ip_stats(&seq, 10, true);
        assert_eq!(runr, 5); // TTTTT adjacent to the site
        assert_eq!(a20r, 5);
        assert!(!is_internal_priming(a20r, runr));
    }

    #[test]
    fn bh_is_monotone() {
        let q = bh_fdr(&[0.01, 0.04, 0.03, 0.5]);
        assert!(q[0] <= q[2] && q[2] <= q[1] && q[1] <= q[3]);
        assert!((q[3] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn paired_test_uses_all_exact_sign_assignments() {
        let (_, p_t, p_flip) = paired_test(&[0.3, 0.2, 0.4, 0.1, 0.5, 0.2, 0.3, 0.4]).unwrap();
        assert!(p_t < 1e-3);
        assert!((p_flip - 2.0 / 256.0).abs() < 1e-12);
        assert!(paired_test(&[0.1]).is_none());
        assert!(paired_test(&[0.1, f64::NAN]).is_none());
    }

    #[test]
    fn beta_binomial_contrast_orders_null_and_balanced_signal() {
        let null_a = [(10, 100), (22, 200), (14, 140), (18, 180)];
        let null_b = [(11, 100), (20, 200), (15, 140), (17, 180)];
        let null = beta_binomial_contrast(&null_a, &null_b).unwrap();
        assert!(null.p_value > 0.5);
        assert!((null.condition_b_mean - null.condition_a_mean).abs() < 0.03);

        let signal_a = [(2, 100), (4, 120), (3, 110), (2, 90)];
        let signal_b = [(88, 100), (106, 120), (97, 110), (80, 90)];
        let signal = beta_binomial_contrast(&signal_a, &signal_b).unwrap();
        assert!(signal.p_value < 1e-4);
        assert!(signal.condition_b_mean - signal.condition_a_mean > 0.75);
        assert!(signal.likelihood_ratio > null.likelihood_ratio);
    }

    #[test]
    fn beta_binomial_rejects_invalid_or_empty_rows() {
        assert!(beta_binomial_contrast(&[], &[(1, 2)]).is_none());
        assert!(beta_binomial_contrast(&[(1, 0)], &[(1, 2)]).is_none());
        assert!(beta_binomial_contrast(&[(3, 2)], &[(1, 2)]).is_none());
    }
}
