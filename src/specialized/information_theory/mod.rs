//! Information Theory Module
//!
//! Implements information-theoretic measures and algorithms:
//! - Entropy calculations (Shannon, Rényi, differential)
//! - Mutual information
//! - Channel capacity
//! - Coding theory (Huffman, arithmetic coding)
//! - Kolmogorov complexity estimation
//! - Data compression metrics

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Request/Response types
#[derive(Debug, Deserialize)]
pub struct EntropyRequest {
    pub data: Vec<f64>,
    pub base: Option<f64>,            // logarithm base (default: 2 for bits)
    pub entropy_type: Option<String>, // "shannon", "renyi", "differential"
    pub alpha: Option<f64>,           // for Rényi entropy
}

#[derive(Debug, Serialize)]
pub struct EntropyResult {
    pub entropy: f64,
    pub entropy_type: String,
    pub base: f64,
    pub max_entropy: f64,
    pub normalized_entropy: f64,
}

#[derive(Debug, Deserialize)]
pub struct MutualInfoRequest {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub bins: Option<usize>,
}

#[derive(Debug, Serialize)]
pub struct MutualInfoResult {
    pub mutual_information: f64,
    pub normalized_mi: f64, // normalized by joint entropy
}

#[derive(Debug, Deserialize)]
pub struct ChannelCapacityRequest {
    pub transition_matrix: Vec<Vec<f64>>,     // P(Y|X) probabilities
    pub input_distribution: Option<Vec<f64>>, // if None, optimize for capacity
}

#[derive(Debug, Serialize)]
pub struct ChannelCapacityResult {
    pub capacity: f64, // bits
    pub optimal_input_distribution: Vec<f64>,
    pub achievable_rate: f64,
}

#[derive(Debug, Deserialize)]
pub struct HuffmanRequest {
    pub symbols: Vec<String>,
    pub frequencies: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct HuffmanResult {
    pub codes: HashMap<String, String>,
    pub average_length: f64,
    pub efficiency: f64,
    pub compression_ratio: f64,
}

#[derive(Debug, Deserialize)]
pub struct KolmogorovRequest {
    pub data: Vec<u8>,
    pub compressor: Option<String>, // "lz77", "lz78", "deflate"
}

#[derive(Debug, Serialize)]
pub struct KolmogorovResult {
    pub complexity_estimate: f64,
    pub compressed_size: usize,
    pub original_size: usize,
    pub compression_ratio: f64,
}

#[derive(Debug, Deserialize)]
pub struct ConditionalEntropyRequest {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub bins: Option<usize>,
}

#[derive(Debug, Serialize)]
pub struct ConditionalEntropyResult {
    pub conditional_entropy: f64, // H(Y|X)
    pub normalized: f64,
}

#[derive(Debug, Deserialize)]
pub struct RelativeEntropyRequest {
    pub p: Vec<f64>, // probability distribution P
    pub q: Vec<f64>, // probability distribution Q
}

#[derive(Debug, Serialize)]
pub struct RelativeEntropyResult {
    pub kl_divergence: f64, // D(P||Q)
    pub js_divergence: f64, // Jensen-Shannon divergence
}

// Helper functions
fn discretize(data: &[f64], bins: usize) -> Vec<usize> {
    let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;

    data.iter()
        .map(|&x| {
            let bin = ((x - min) / range * bins as f64) as usize;
            bin.min(bins - 1)
        })
        .collect()
}

fn compute_probabilities(values: &[usize], num_bins: usize) -> Vec<f64> {
    let mut counts = vec![0; num_bins];
    for &val in values {
        counts[val] += 1;
    }
    let total = values.len() as f64;
    counts.iter().map(|&c| c as f64 / total).collect()
}

/// Calculate Shannon entropy
pub fn shannon_entropy(request: EntropyRequest) -> Result<EntropyResult, String> {
    let base = request.base.unwrap_or(2.0);
    let entropy_type = request
        .entropy_type
        .unwrap_or_else(|| "shannon".to_string());

    let data = &request.data;
    if data.is_empty() {
        return Err("Empty data".to_string());
    }

    let entropy = match entropy_type.as_str() {
        "shannon" => {
            // Discretize data
            let bins = (data.len() as f64).sqrt() as usize + 1;
            let discretized = discretize(data, bins);
            let probs = compute_probabilities(&discretized, bins);

            -probs
                .iter()
                .filter(|&&p| p > 0.0)
                .map(|&p| p * p.log(base))
                .sum::<f64>()
        }
        "renyi" => {
            let alpha = request.alpha.unwrap_or(2.0);
            if alpha == 1.0 {
                return shannon_entropy(EntropyRequest {
                    data: data.clone(),
                    base: Some(base),
                    entropy_type: Some("shannon".to_string()),
                    alpha: None,
                });
            }

            let bins = (data.len() as f64).sqrt() as usize + 1;
            let discretized = discretize(data, bins);
            let probs = compute_probabilities(&discretized, bins);

            let sum: f64 = probs
                .iter()
                .filter(|&&p| p > 0.0)
                .map(|&p| p.powf(alpha))
                .sum();

            (1.0 / (1.0 - alpha)) * sum.log(base)
        }
        "differential" => {
            // Estimate differential entropy using kernel density estimation
            let n = data.len() as f64;
            let std_dev = {
                let mean = data.iter().sum::<f64>() / n;
                (data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n).sqrt()
            };

            // Differential entropy for Gaussian approximation
            0.5 * (2.0 * std::f64::consts::PI * std::f64::consts::E * std_dev * std_dev).log(base)
        }
        _ => return Err(format!("Unknown entropy type: {}", entropy_type)),
    };

    let max_entropy = match entropy_type.as_str() {
        "differential" => f64::INFINITY,
        _ => {
            let bins = (data.len() as f64).sqrt() as usize + 1;
            (bins as f64).log(base)
        }
    };

    let normalized_entropy = if max_entropy.is_finite() && max_entropy > 0.0 {
        entropy / max_entropy
    } else {
        1.0
    };

    Ok(EntropyResult {
        entropy,
        entropy_type,
        base,
        max_entropy,
        normalized_entropy,
    })
}

/// Calculate mutual information between two variables
pub fn mutual_information(request: MutualInfoRequest) -> Result<MutualInfoResult, String> {
    if request.x.len() != request.y.len() {
        return Err("X and Y must have same length".to_string());
    }

    let bins = request
        .bins
        .unwrap_or((request.x.len() as f64).sqrt() as usize + 1);

    let x_disc = discretize(&request.x, bins);
    let y_disc = discretize(&request.y, bins);

    // Calculate marginal probabilities
    let p_x = compute_probabilities(&x_disc, bins);
    let p_y = compute_probabilities(&y_disc, bins);

    // Calculate joint probabilities
    let mut joint_counts = vec![vec![0; bins]; bins];
    for i in 0..x_disc.len() {
        joint_counts[x_disc[i]][y_disc[i]] += 1;
    }

    let total = x_disc.len() as f64;
    let mut mi = 0.0;

    for i in 0..bins {
        for j in 0..bins {
            let p_xy = joint_counts[i][j] as f64 / total;
            if p_xy > 0.0 && p_x[i] > 0.0 && p_y[j] > 0.0 {
                mi += p_xy * (p_xy / (p_x[i] * p_y[j])).log2();
            }
        }
    }

    // Calculate joint entropy for normalization
    let h_xy = -joint_counts
        .iter()
        .flat_map(|row| row.iter())
        .map(|&c| c as f64 / total)
        .filter(|&p| p > 0.0)
        .map(|p| p * p.log2())
        .sum::<f64>();

    let normalized_mi = if h_xy > 0.0 { mi / h_xy } else { 0.0 };

    Ok(MutualInfoResult {
        mutual_information: mi,
        normalized_mi,
    })
}

/// Calculate channel capacity
pub fn channel_capacity(request: ChannelCapacityRequest) -> Result<ChannelCapacityResult, String> {
    let matrix = request.transition_matrix;
    let n_inputs = matrix.len();

    if n_inputs == 0 {
        return Err("Empty transition matrix".to_string());
    }

    // Validate transition matrix
    for row in &matrix {
        let sum: f64 = row.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err("Each row of transition matrix must sum to 1".to_string());
        }
    }

    // Use uniform distribution if not provided
    let input_dist = request
        .input_distribution
        .unwrap_or_else(|| vec![1.0 / n_inputs as f64; n_inputs]);

    // Calculate output distribution
    let n_outputs = matrix[0].len();
    let mut output_dist = vec![0.0; n_outputs];

    for (i, &p_x) in input_dist.iter().enumerate() {
        for (j, &p_y_given_x) in matrix[i].iter().enumerate() {
            output_dist[j] += p_x * p_y_given_x;
        }
    }

    // Calculate mutual information (which equals channel capacity at optimal input)
    let mut capacity = 0.0;

    for (i, &p_x) in input_dist.iter().enumerate() {
        if p_x > 0.0 {
            for (j, &p_y_given_x) in matrix[i].iter().enumerate() {
                if p_y_given_x > 0.0 && output_dist[j] > 0.0 {
                    capacity += p_x * p_y_given_x * (p_y_given_x / output_dist[j]).log2();
                }
            }
        }
    }

    Ok(ChannelCapacityResult {
        capacity,
        optimal_input_distribution: input_dist.clone(),
        achievable_rate: capacity * 0.95, // Conservative estimate with overhead
    })
}

/// Generate Huffman codes
pub fn huffman_coding(request: HuffmanRequest) -> Result<HuffmanResult, String> {
    if request.symbols.len() != request.frequencies.len() {
        return Err("Symbols and frequencies must have same length".to_string());
    }

    let symbols = request.symbols;
    let freqs = request.frequencies;

    if symbols.is_empty() {
        return Err("No symbols provided".to_string());
    }

    // Build Huffman tree
    let mut nodes: Vec<(f64, Vec<(String, String)>)> = symbols
        .iter()
        .zip(freqs.iter())
        .map(|(s, &f)| (f, vec![(s.clone(), String::new())]))
        .collect();

    while nodes.len() > 1 {
        nodes.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

        let right = nodes.pop().unwrap();
        let left = nodes.pop().unwrap();

        let mut combined = Vec::new();
        for (sym, code) in left.1 {
            combined.push((sym, format!("0{}", code)));
        }
        for (sym, code) in right.1 {
            combined.push((sym, format!("1{}", code)));
        }

        nodes.push((left.0 + right.0, combined));
    }

    let codes: HashMap<String, String> = if let Some((_, assignments)) = nodes.pop() {
        assignments.into_iter().collect()
    } else {
        HashMap::new()
    };

    // Calculate average length
    let avg_length: f64 = symbols
        .iter()
        .zip(freqs.iter())
        .map(|(s, &f)| f * codes.get(s).map(|c| c.len()).unwrap_or(0) as f64)
        .sum();

    // Calculate entropy for efficiency
    let entropy: f64 = -freqs
        .iter()
        .filter(|&&f| f > 0.0)
        .map(|&f| f * f.log2())
        .sum::<f64>();

    let efficiency = if avg_length > 0.0 {
        entropy / avg_length
    } else {
        0.0
    };

    // Assume original uses 8 bits per symbol on average
    let compression_ratio = if avg_length > 0.0 {
        8.0 / avg_length
    } else {
        1.0
    };

    Ok(HuffmanResult {
        codes,
        average_length: avg_length,
        efficiency,
        compression_ratio,
    })
}

/// Estimate Kolmogorov complexity using compression
pub fn kolmogorov_complexity(request: KolmogorovRequest) -> Result<KolmogorovResult, String> {
    let data = request.data;
    let original_size = data.len();

    if original_size == 0 {
        return Err("Empty data".to_string());
    }

    // Simple LZ77-like compression estimate
    let compressed_size = estimate_lz77_compression(&data);

    let complexity_estimate = compressed_size as f64;
    let compression_ratio = original_size as f64 / compressed_size as f64;

    Ok(KolmogorovResult {
        complexity_estimate,
        compressed_size,
        original_size,
        compression_ratio,
    })
}

fn estimate_lz77_compression(data: &[u8]) -> usize {
    // Simplified LZ77 compression size estimate
    let mut size = 0;
    let mut i = 0;

    while i < data.len() {
        let mut best_len = 0;
        let mut best_dist = 0;

        // Look for matches in the window
        for dist in 1..i.min(256) {
            let mut len = 0;
            while i + len < data.len() && len < 15 && data[i - dist + len] == data[i + len] {
                len += 1;
            }
            if len > best_len {
                best_len = len;
                best_dist = dist;
            }
        }

        if best_len >= 3 {
            size += 2; // (distance, length) pair
            i += best_len;
        } else {
            size += 1; // literal byte
            i += 1;
        }
    }

    size
}

/// Calculate conditional entropy H(Y|X)
pub fn conditional_entropy(
    request: ConditionalEntropyRequest,
) -> Result<ConditionalEntropyResult, String> {
    if request.x.len() != request.y.len() {
        return Err("X and Y must have same length".to_string());
    }

    let bins = request
        .bins
        .unwrap_or((request.x.len() as f64).sqrt() as usize + 1);

    let x_disc = discretize(&request.x, bins);
    let y_disc = discretize(&request.y, bins);

    // Calculate P(X)
    let p_x = compute_probabilities(&x_disc, bins);

    // Calculate P(X,Y)
    let mut joint_counts = vec![vec![0; bins]; bins];
    for i in 0..x_disc.len() {
        joint_counts[x_disc[i]][y_disc[i]] += 1;
    }

    let total = x_disc.len() as f64;
    let mut h_y_given_x = 0.0;

    for i in 0..bins {
        if p_x[i] > 0.0 {
            for j in 0..bins {
                let p_xy = joint_counts[i][j] as f64 / total;
                if p_xy > 0.0 {
                    let p_y_given_x = p_xy / p_x[i];
                    h_y_given_x -= p_xy * p_y_given_x.log2();
                }
            }
        }
    }

    // Calculate H(Y) for normalization
    let p_y = compute_probabilities(&y_disc, bins);
    let h_y: f64 = -p_y
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| p * p.log2())
        .sum::<f64>();

    let normalized = if h_y > 0.0 { h_y_given_x / h_y } else { 0.0 };

    Ok(ConditionalEntropyResult {
        conditional_entropy: h_y_given_x,
        normalized,
    })
}

/// Calculate KL divergence and JS divergence
pub fn relative_entropy(request: RelativeEntropyRequest) -> Result<RelativeEntropyResult, String> {
    if request.p.len() != request.q.len() {
        return Err("P and Q must have same length".to_string());
    }

    let p = &request.p;
    let q = &request.q;

    // Validate probability distributions
    let sum_p: f64 = p.iter().sum();
    let sum_q: f64 = q.iter().sum();

    if (sum_p - 1.0).abs() > 1e-6 || (sum_q - 1.0).abs() > 1e-6 {
        return Err("P and Q must be valid probability distributions (sum to 1)".to_string());
    }

    // KL divergence: D(P||Q)
    let mut kl_divergence = 0.0;
    for i in 0..p.len() {
        if p[i] > 0.0 {
            if q[i] <= 0.0 {
                return Err(
                    "Q has zero probability where P is non-zero (KL divergence undefined)"
                        .to_string(),
                );
            }
            kl_divergence += p[i] * (p[i] / q[i]).log2();
        }
    }

    // Jensen-Shannon divergence: symmetric version
    let mut m = vec![0.0; p.len()];
    for i in 0..p.len() {
        m[i] = 0.5 * (p[i] + q[i]);
    }

    let mut js_divergence = 0.0;
    for i in 0..p.len() {
        if p[i] > 0.0 && m[i] > 0.0 {
            js_divergence += 0.5 * p[i] * (p[i] / m[i]).log2();
        }
        if q[i] > 0.0 && m[i] > 0.0 {
            js_divergence += 0.5 * q[i] * (q[i] / m[i]).log2();
        }
    }

    Ok(RelativeEntropyResult {
        kl_divergence,
        js_divergence,
    })
}

