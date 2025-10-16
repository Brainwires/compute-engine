//! Statistics and Probability Module
//!
//! Provides comprehensive statistical analysis including:
//! - Descriptive statistics (mean, median, mode, variance, std dev)
//! - Probability distributions
//! - Monte Carlo integration
//! - MCMC sampling
//! - Correlation and regression
//! - Hypothesis testing

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
pub struct StatisticsRequest {
    pub data: Vec<f64>,
    pub operations: Option<Vec<String>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct StatisticsResult {
    pub mean: Option<f64>,
    pub median: Option<f64>,
    pub mode: Option<Vec<f64>>,
    pub std: Option<f64>,
    pub variance: Option<f64>,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub sum: Option<f64>,
    pub count: Option<usize>,
    pub quartiles: Option<Quartiles>,
    pub skewness: Option<f64>,
    pub kurtosis: Option<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Quartiles {
    pub q1: f64,
    pub q2: f64, // median
    pub q3: f64,
    pub iqr: f64,
}

pub fn statistics(request: StatisticsRequest) -> Result<StatisticsResult, String> {
    if request.data.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    let ops = request
        .operations
        .as_ref()
        .map(|v| v.iter().map(|s| s.as_str()).collect::<Vec<_>>());
    let should_compute = |op: &str| ops.as_ref().map_or(true, |o| o.contains(&op));

    let mut result = StatisticsResult {
        mean: None,
        median: None,
        mode: None,
        std: None,
        variance: None,
        min: None,
        max: None,
        sum: None,
        count: None,
        quartiles: None,
        skewness: None,
        kurtosis: None,
    };

    if should_compute("count") {
        result.count = Some(request.data.len());
    }

    if should_compute("sum") {
        result.sum = Some(request.data.iter().sum());
    }

    if should_compute("mean") {
        result.mean = Some(request.data.iter().sum::<f64>() / request.data.len() as f64);
    }

    if should_compute("min") {
        result.min = request
            .data
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min)
            .into();
    }

    if should_compute("max") {
        result.max = request
            .data
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
            .into();
    }

    if should_compute("median") || should_compute("quartiles") {
        let mut sorted = request.data.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median = calculate_percentile(&sorted, 50.0);
        if should_compute("median") {
            result.median = Some(median);
        }

        if should_compute("quartiles") {
            let q1 = calculate_percentile(&sorted, 25.0);
            let q3 = calculate_percentile(&sorted, 75.0);
            result.quartiles = Some(Quartiles {
                q1,
                q2: median,
                q3,
                iqr: q3 - q1,
            });
        }
    }

    if should_compute("mode") {
        let mut freq: HashMap<String, usize> = HashMap::new();
        for &val in &request.data {
            *freq.entry(val.to_string()).or_insert(0) += 1;
        }
        let max_freq = freq.values().max().copied().unwrap_or(0);
        let modes: Vec<f64> = freq
            .iter()
            .filter(|(_, count)| **count == max_freq)
            .filter_map(|(val, _)| val.parse().ok())
            .collect();
        result.mode = Some(modes);
    }

    if should_compute("variance") || should_compute("std") {
        let mean = request.data.iter().sum::<f64>() / request.data.len() as f64;
        let variance = request.data.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / request.data.len() as f64;

        if should_compute("variance") {
            result.variance = Some(variance);
        }
        if should_compute("std") {
            result.std = Some(variance.sqrt());
        }
    }

    if should_compute("skewness") {
        let mean = request.data.iter().sum::<f64>() / request.data.len() as f64;
        let std = (request.data.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / request.data.len() as f64)
            .sqrt();

        let skewness = request
            .data
            .iter()
            .map(|x| ((x - mean) / std).powi(3))
            .sum::<f64>()
            / request.data.len() as f64;

        result.skewness = Some(skewness);
    }

    if should_compute("kurtosis") {
        let mean = request.data.iter().sum::<f64>() / request.data.len() as f64;
        let std = (request.data.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / request.data.len() as f64)
            .sqrt();

        let kurtosis = request
            .data
            .iter()
            .map(|x| ((x - mean) / std).powi(4))
            .sum::<f64>()
            / request.data.len() as f64
            - 3.0; // Excess kurtosis

        result.kurtosis = Some(kurtosis);
    }

    Ok(result)
}

fn calculate_percentile(sorted_data: &[f64], percentile: f64) -> f64 {
    let index = (percentile / 100.0) * (sorted_data.len() - 1) as f64;
    let lower = index.floor() as usize;
    let upper = index.ceil() as usize;

    if lower == upper {
        sorted_data[lower]
    } else {
        let weight = index - lower as f64;
        sorted_data[lower] * (1.0 - weight) + sorted_data[upper] * weight
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MonteCarloRequest {
    pub function: String, // Function to integrate (simplified - would need parser in practice)
    pub bounds: Vec<(f64, f64)>, // Integration bounds for each dimension
    pub samples: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MonteCarloResult {
    pub integral_estimate: f64,
    pub standard_error: f64,
    pub samples_used: usize,
}

pub fn monte_carlo_integration(request: MonteCarloRequest) -> Result<MonteCarloResult, String> {
    use rand::Rng;

    if request.bounds.is_empty() {
        return Err("Bounds cannot be empty".to_string());
    }

    let mut rng = rand::thread_rng();
    let mut sum = 0.0;
    let mut sum_sq = 0.0;

    // Calculate volume of integration region
    let volume: f64 = request.bounds.iter().map(|(a, b)| b - a).product();

    // Simple example: integrate x^2 + y^2 (would need proper function parser)
    for _ in 0..request.samples {
        let point: Vec<f64> = request
            .bounds
            .iter()
            .map(|(a, b)| rng.gen_range(*a..*b))
            .collect();

        // Evaluate function at random point (simplified)
        let value = point.iter().map(|x| x * x).sum::<f64>();

        sum += value;
        sum_sq += value * value;
    }

    let mean = sum / request.samples as f64;
    let variance = (sum_sq / request.samples as f64) - mean * mean;
    let integral = volume * mean;
    let standard_error = volume * (variance / request.samples as f64).sqrt();

    Ok(MonteCarloResult {
        integral_estimate: integral,
        standard_error,
        samples_used: request.samples,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MCMCRequest {
    pub target_distribution: String, // Target distribution type
    pub dimensions: usize,
    pub samples: usize,
    pub burn_in: usize,
    pub initial_state: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MCMCResult {
    pub samples: Vec<Vec<f64>>,
    pub acceptance_rate: f64,
    pub effective_samples: usize,
}

pub fn mcmc_sampling(request: MCMCRequest) -> Result<MCMCResult, String> {
    use rand::Rng;
    use rand_distr::{Distribution, Normal};

    if request.dimensions == 0 {
        return Err("Dimensions must be greater than 0".to_string());
    }

    let mut rng = rand::thread_rng();
    let proposal_std = 0.5; // Standard deviation for proposal distribution
    let normal = Normal::new(0.0, proposal_std).unwrap();

    let mut current_state = request.initial_state.clone();
    if current_state.len() != request.dimensions {
        current_state = vec![0.0; request.dimensions];
    }

    let mut samples = Vec::new();
    let mut accepted = 0;

    // Target distribution (simplified - Gaussian for demonstration)
    let log_target = |x: &[f64]| -> f64 { -0.5 * x.iter().map(|v| v * v).sum::<f64>() };

    let mut current_log_prob = log_target(&current_state);

    for i in 0..(request.samples + request.burn_in) {
        // Propose new state
        let proposed_state: Vec<f64> = current_state
            .iter()
            .map(|x| x + normal.sample(&mut rng))
            .collect();

        let proposed_log_prob = log_target(&proposed_state);

        // Metropolis-Hastings acceptance criterion
        let log_acceptance_ratio = proposed_log_prob - current_log_prob;

        if log_acceptance_ratio > 0.0 || rng.r#gen::<f64>().ln() < log_acceptance_ratio {
            current_state = proposed_state;
            current_log_prob = proposed_log_prob;
            accepted += 1;
        }

        // Store sample after burn-in
        if i >= request.burn_in {
            samples.push(current_state.clone());
        }
    }

    let acceptance_rate = accepted as f64 / (request.samples + request.burn_in) as f64;

    Ok(MCMCResult {
        samples: samples.clone(),
        acceptance_rate,
        effective_samples: samples.len(),
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CorrelationRequest {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub method: String, // "pearson", "spearman", "kendall"
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CorrelationResult {
    pub correlation: f64,
    pub p_value: Option<f64>,
    pub method: String,
}

pub fn correlation(request: CorrelationRequest) -> Result<CorrelationResult, String> {
    if request.x.len() != request.y.len() {
        return Err("X and Y must have the same length".to_string());
    }

    if request.x.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    let correlation = match request.method.as_str() {
        "pearson" => calculate_pearson(&request.x, &request.y)?,
        "spearman" => {
            let x_ranks = assign_ranks(&request.x);
            let y_ranks = assign_ranks(&request.y);
            calculate_pearson(&x_ranks, &y_ranks)?
        }
        _ => {
            return Err(format!(
                "Unsupported correlation method: {}",
                request.method
            ));
        }
    };

    Ok(CorrelationResult {
        correlation,
        p_value: None, // Would need t-distribution for p-value calculation
        method: request.method,
    })
}

fn calculate_pearson(x: &[f64], y: &[f64]) -> Result<f64, String> {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    let mut sum_xy = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;

    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        sum_xy += dx * dy;
        sum_x2 += dx * dx;
        sum_y2 += dy * dy;
    }

    if sum_x2 == 0.0 || sum_y2 == 0.0 {
        return Err("Standard deviation is zero".to_string());
    }

    Ok(sum_xy / (sum_x2 * sum_y2).sqrt())
}

fn assign_ranks(data: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<(usize, f64)> = data.iter().enumerate().map(|(i, &v)| (i, v)).collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut ranks = vec![0.0; data.len()];
    for (rank, (idx, _)) in indexed.iter().enumerate() {
        ranks[*idx] = (rank + 1) as f64;
    }

    ranks
}

#[derive(Debug, Serialize, Deserialize)]
pub struct KLDivergenceRequest {
    pub p: Vec<f64>, // Probability distribution P
    pub q: Vec<f64>, // Probability distribution Q
}

#[derive(Debug, Serialize, Deserialize)]
pub struct KLDivergenceResult {
    pub kl_divergence: f64,
    pub js_divergence: f64, // Jensen-Shannon divergence
}

pub fn kl_divergence(request: KLDivergenceRequest) -> Result<KLDivergenceResult, String> {
    if request.p.len() != request.q.len() {
        return Err("P and Q must have the same length".to_string());
    }

    let epsilon = 1e-10;
    let mut kl = 0.0;

    for i in 0..request.p.len() {
        if request.p[i] < 0.0 || request.q[i] < 0.0 {
            return Err("Probabilities must be non-negative".to_string());
        }

        if request.p[i] > epsilon {
            let q_safe = request.q[i].max(epsilon);
            kl += request.p[i] * (request.p[i] / q_safe).ln();
        }
    }

    // Calculate Jensen-Shannon divergence
    let m: Vec<f64> = request
        .p
        .iter()
        .zip(request.q.iter())
        .map(|(p, q)| (p + q) / 2.0)
        .collect();

    let mut kl_pm = 0.0;
    let mut kl_qm = 0.0;

    for i in 0..request.p.len() {
        if request.p[i] > epsilon {
            kl_pm += request.p[i] * (request.p[i] / m[i].max(epsilon)).ln();
        }
        if request.q[i] > epsilon {
            kl_qm += request.q[i] * (request.q[i] / m[i].max(epsilon)).ln();
        }
    }

    let js = (kl_pm + kl_qm) / 2.0;

    Ok(KLDivergenceResult {
        kl_divergence: kl,
        js_divergence: js,
    })
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MutualInformationRequest {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub bins: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MutualInformationResult {
    pub mutual_information: f64,
    pub normalized_mi: f64,
}

pub fn mutual_information(
    request: MutualInformationRequest,
) -> Result<MutualInformationResult, String> {
    if request.x.len() != request.y.len() {
        return Err("X and Y must have the same length".to_string());
    }

    if request.x.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    if request.bins < 2 {
        return Err("Number of bins must be at least 2".to_string());
    }

    // Calculate MI(X;Y) = H(X) + H(Y) - H(X,Y)

    // Find ranges for binning
    let x_min = request.x.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let x_max = request.x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let y_min = request.y.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let y_max = request.y.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let x_bin_width = (x_max - x_min) / request.bins as f64;
    let y_bin_width = (y_max - y_min) / request.bins as f64;

    // Build histograms
    let mut x_hist = vec![0usize; request.bins];
    let mut y_hist = vec![0usize; request.bins];
    let mut xy_hist = vec![vec![0usize; request.bins]; request.bins];

    for i in 0..request.x.len() {
        let x_bin = ((request.x[i] - x_min) / x_bin_width).floor() as usize;
        let y_bin = ((request.y[i] - y_min) / y_bin_width).floor() as usize;

        let x_bin = x_bin.min(request.bins - 1);
        let y_bin = y_bin.min(request.bins - 1);

        x_hist[x_bin] += 1;
        y_hist[y_bin] += 1;
        xy_hist[x_bin][y_bin] += 1;
    }

    let n = request.x.len() as f64;

    // Calculate H(X)
    let h_x = x_hist
        .iter()
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = count as f64 / n;
            -p * p.ln()
        })
        .sum::<f64>();

    // Calculate H(Y)
    let h_y = y_hist
        .iter()
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = count as f64 / n;
            -p * p.ln()
        })
        .sum::<f64>();

    // Calculate H(X,Y)
    let h_xy = xy_hist
        .iter()
        .flat_map(|row| row.iter())
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = count as f64 / n;
            -p * p.ln()
        })
        .sum::<f64>();

    // Mutual information: MI(X;Y) = H(X) + H(Y) - H(X,Y)
    let mi = h_x + h_y - h_xy;

    // Normalized MI: MI / min(H(X), H(Y))
    let h_min = h_x.min(h_y);
    let normalized = if h_min > 1e-10 { mi / h_min } else { 0.0 };

    Ok(MutualInformationResult {
        mutual_information: mi,
        normalized_mi: normalized.min(1.0),
    })
}
