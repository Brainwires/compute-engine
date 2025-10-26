//! Clustering Algorithms
//!
//! Unsupervised learning methods for grouping data:
//! - K-means clustering
//! - K-means++ initialization
//! - Hierarchical clustering (future)

use rand::Rng;

/// Result of K-means clustering
#[derive(Debug, Clone)]
pub struct KMeansResult {
    /// Cluster centers (k centroids, each with n_features dimensions)
    pub centroids: Vec<Vec<f64>>,
    /// Cluster assignment for each sample
    pub labels: Vec<usize>,
    /// Sum of squared distances to nearest centroid (inertia)
    pub inertia: f64,
    /// Number of iterations performed
    pub n_iterations: usize,
    /// Whether the algorithm converged
    pub converged: bool,
}

/// K-means clustering
///
/// Partitions data into k clusters by minimizing within-cluster variance.
///
/// # Arguments
/// * `data` - Samples to cluster (each row is a sample)
/// * `k` - Number of clusters
/// * `max_iterations` - Maximum iterations (default: 300)
/// * `tolerance` - Convergence tolerance (default: 1e-4)
/// * `init_method` - Initialization: "random" or "kmeans++" (default: "kmeans++")
///
/// # Algorithm
/// 1. Initialize k centroids
/// 2. Assign each point to nearest centroid
/// 3. Update centroids as mean of assigned points
/// 4. Repeat until convergence or max iterations
///
/// # Example
/// ```
/// use computational_engine::specialized::machine_learning::clustering::kmeans;
///
/// let data = vec![
///     vec![1.0, 2.0],
///     vec![1.5, 1.8],
///     vec![5.0, 8.0],
///     vec![8.0, 8.0],
///     vec![1.0, 0.6],
///     vec![9.0, 11.0],
/// ];
///
/// let result = kmeans(&data, 2, Some(100), Some(1e-4), Some("kmeans++")).unwrap();
/// assert_eq!(result.centroids.len(), 2);
/// assert_eq!(result.labels.len(), 6);
/// ```
pub fn kmeans(
    data: &[Vec<f64>],
    k: usize,
    max_iterations: Option<usize>,
    tolerance: Option<f64>,
    init_method: Option<&str>,
) -> Result<KMeansResult, String> {
    if data.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    if k == 0 {
        return Err("Number of clusters must be positive".to_string());
    }

    if k > data.len() {
        return Err(format!(
            "Number of clusters ({}) cannot exceed number of samples ({})",
            k,
            data.len()
        ));
    }

    let n_samples = data.len();
    let n_features = data[0].len();
    let max_iter = max_iterations.unwrap_or(300);
    let tol = tolerance.unwrap_or(1e-4);
    let method = init_method.unwrap_or("kmeans++");

    // Initialize centroids
    let mut centroids = match method {
        "kmeans++" => kmeans_plus_plus_init(data, k),
        "random" => random_init(data, k),
        _ => return Err(format!("Unknown initialization method: {}", method)),
    };

    let mut labels = vec![0; n_samples];
    let mut converged = false;

    for iteration in 0..max_iter {
        // Assignment step: assign each point to nearest centroid
        let mut new_labels = vec![0; n_samples];
        for (i, sample) in data.iter().enumerate() {
            let mut min_dist = f64::INFINITY;
            let mut best_cluster = 0;

            for (j, centroid) in centroids.iter().enumerate() {
                let dist = euclidean_distance(sample, centroid);
                if dist < min_dist {
                    min_dist = dist;
                    best_cluster = j;
                }
            }

            new_labels[i] = best_cluster;
        }

        // Update step: compute new centroids
        let mut new_centroids = vec![vec![0.0; n_features]; k];
        let mut counts = vec![0; k];

        for (i, sample) in data.iter().enumerate() {
            let cluster = new_labels[i];
            counts[cluster] += 1;
            for (j, &val) in sample.iter().enumerate() {
                new_centroids[cluster][j] += val;
            }
        }

        // Average to get centroids
        for (cluster, count) in counts.iter().enumerate() {
            if *count > 0 {
                for j in 0..n_features {
                    new_centroids[cluster][j] /= *count as f64;
                }
            } else {
                // Empty cluster: reinitialize randomly
                let random_idx = rand::thread_rng().gen_range(0..n_samples);
                new_centroids[cluster] = data[random_idx].clone();
            }
        }

        // Check convergence
        let max_shift = centroids
            .iter()
            .zip(&new_centroids)
            .map(|(old, new)| euclidean_distance(old, new))
            .fold(0.0, f64::max);

        centroids = new_centroids;
        labels = new_labels;

        if max_shift < tol {
            converged = true;
            // Compute final inertia
            let inertia = compute_inertia(data, &centroids, &labels);
            return Ok(KMeansResult {
                centroids,
                labels,
                inertia,
                n_iterations: iteration + 1,
                converged,
            });
        }
    }

    // Max iterations reached
    let inertia = compute_inertia(data, &centroids, &labels);
    Ok(KMeansResult {
        centroids,
        labels,
        inertia,
        n_iterations: max_iter,
        converged,
    })
}

/// K-means++ initialization
///
/// Chooses initial centroids that are spread out, improving convergence.
fn kmeans_plus_plus_init(data: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let mut rng = rand::thread_rng();
    let mut centroids = Vec::new();

    // Choose first centroid randomly
    let first_idx = rng.gen_range(0..data.len());
    centroids.push(data[first_idx].clone());

    // Choose remaining centroids
    for target_count in 2..=k {
        // Compute distance to nearest centroid for each point
        let mut distances = Vec::new();
        let mut total_distance = 0.0;

        for sample in data {
            let min_dist = centroids
                .iter()
                .map(|c| euclidean_distance_squared(sample, c))
                .fold(f64::INFINITY, f64::min);

            distances.push(min_dist);
            total_distance += min_dist;
        }

        // Choose next centroid with probability proportional to distance^2
        let threshold = rng.gen_range(0.0..total_distance);
        let mut cumsum = 0.0;
        let mut added = false;

        for (i, &dist) in distances.iter().enumerate() {
            cumsum += dist;
            if cumsum >= threshold {
                centroids.push(data[i].clone());
                added = true;
                break;
            }
        }

        // Fallback if loop didn't break (shouldn't happen)
        if !added {
            centroids.push(data[rng.gen_range(0..data.len())].clone());
        }
    }

    centroids
}

/// Random initialization
fn random_init(data: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let mut rng = rand::thread_rng();
    let mut centroids = Vec::new();
    let mut used_indices = std::collections::HashSet::new();

    while centroids.len() < k {
        let idx = rng.gen_range(0..data.len());
        if used_indices.insert(idx) {
            centroids.push(data[idx].clone());
        }
    }

    centroids
}

/// Euclidean distance between two points
pub fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    euclidean_distance_squared(a, b).sqrt()
}

/// Squared Euclidean distance
pub fn euclidean_distance_squared(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(&x, &y)| (x - y) * (x - y))
        .sum()
}

/// Compute inertia (sum of squared distances to nearest centroid)
fn compute_inertia(data: &[Vec<f64>], centroids: &[Vec<f64>], labels: &[usize]) -> f64 {
    data.iter()
        .zip(labels)
        .map(|(sample, &label)| euclidean_distance_squared(sample, &centroids[label]))
        .sum()
}

/// Silhouette score for clustering quality
///
/// Measures how similar a point is to its own cluster compared to other clusters.
/// Range: [-1, 1] where 1 is best, 0 is overlapping clusters, -1 is misclassified.
///
/// # Arguments
/// * `data` - Original data
/// * `labels` - Cluster assignments
///
/// # Returns
/// Average silhouette score across all samples
pub fn silhouette_score(data: &[Vec<f64>], labels: &[usize]) -> Result<f64, String> {
    if data.len() != labels.len() {
        return Err("Data and labels must have same length".to_string());
    }

    if data.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    let n_samples = data.len();
    let n_clusters = labels.iter().max().unwrap_or(&0) + 1;

    if n_clusters < 2 {
        return Err("Need at least 2 clusters for silhouette score".to_string());
    }

    let mut total_score = 0.0;

    for i in 0..n_samples {
        let cluster_i = labels[i];

        // Compute a(i): mean distance to points in same cluster
        let mut same_cluster_dist = 0.0;
        let mut same_cluster_count = 0;

        for j in 0..n_samples {
            if i != j && labels[j] == cluster_i {
                same_cluster_dist += euclidean_distance(&data[i], &data[j]);
                same_cluster_count += 1;
            }
        }

        let a_i = if same_cluster_count > 0 {
            same_cluster_dist / same_cluster_count as f64
        } else {
            0.0
        };

        // Compute b(i): min mean distance to points in other clusters
        let mut b_i = f64::INFINITY;

        for cluster in 0..n_clusters {
            if cluster != cluster_i {
                let mut other_cluster_dist = 0.0;
                let mut other_cluster_count = 0;

                for j in 0..n_samples {
                    if labels[j] == cluster {
                        other_cluster_dist += euclidean_distance(&data[i], &data[j]);
                        other_cluster_count += 1;
                    }
                }

                if other_cluster_count > 0 {
                    let mean_dist = other_cluster_dist / other_cluster_count as f64;
                    b_i = b_i.min(mean_dist);
                }
            }
        }

        // Silhouette score for point i
        let s_i = if a_i < b_i {
            1.0 - a_i / b_i
        } else if a_i > b_i {
            b_i / a_i - 1.0
        } else {
            0.0
        };

        total_score += s_i;
    }

    Ok(total_score / n_samples as f64)
}

