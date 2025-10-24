//! Dimensionality Reduction Techniques
//!
//! Methods for reducing the number of features while preserving variance:
//! - PCA (Principal Component Analysis)
//! - Explained variance ratio
//! - Data projection and reconstruction

use nalgebra::{DMatrix, DVector};

/// Result of PCA computation
#[derive(Debug, Clone)]
pub struct PCAResult {
    /// Principal components (eigenvectors as columns)
    pub components: Vec<Vec<f64>>,
    /// Explained variance for each component
    pub explained_variance: Vec<f64>,
    /// Explained variance ratio for each component (sums to 1.0)
    pub explained_variance_ratio: Vec<f64>,
    /// Mean of the original data (used for centering)
    pub mean: Vec<f64>,
    /// Singular values
    pub singular_values: Vec<f64>,
}

/// Perform Principal Component Analysis (PCA)
///
/// Reduces dimensionality by finding orthogonal directions of maximum variance.
///
/// # Arguments
/// * `data` - Matrix where each row is a sample, each column is a feature
/// * `n_components` - Number of principal components to retain (None = keep all)
///
/// # Algorithm
/// 1. Center the data by subtracting the mean
/// 2. Compute SVD: X = U Σ V^T
/// 3. Principal components are columns of V
/// 4. Explained variance = (σ_i^2) / (n - 1)
///
/// # Example
/// ```
/// use computational_engine::specialized::machine_learning::dimensionality_reduction::pca;
///
/// let data = vec![
///     vec![2.5, 2.4],
///     vec![0.5, 0.7],
///     vec![2.2, 2.9],
///     vec![1.9, 2.2],
///     vec![3.1, 3.0],
/// ];
///
/// let result = pca(&data, Some(1)).unwrap();
/// assert_eq!(result.components[0].len(), 2); // First PC has 2 coefficients
/// ```
pub fn pca(data: &[Vec<f64>], n_components: Option<usize>) -> Result<PCAResult, String> {
    if data.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    let n_samples = data.len();
    let n_features = data[0].len();

    if n_samples < 2 {
        return Err("Need at least 2 samples for PCA".to_string());
    }

    // Determine number of components
    let k = n_components.unwrap_or(n_features).min(n_features).min(n_samples);

    // Step 1: Compute mean for each feature
    let mut mean = vec![0.0; n_features];
    for sample in data {
        for (i, &val) in sample.iter().enumerate() {
            mean[i] += val;
        }
    }
    for m in &mut mean {
        *m /= n_samples as f64;
    }

    // Step 2: Center the data
    let mut centered_data = Vec::new();
    for sample in data {
        let mut centered_sample = Vec::new();
        for (i, &val) in sample.iter().enumerate() {
            centered_sample.push(val - mean[i]);
        }
        centered_data.push(centered_sample);
    }

    // Step 3: Convert to nalgebra matrix (samples as rows)
    let mut mat_data = Vec::new();
    for row in &centered_data {
        mat_data.extend(row);
    }
    let x_matrix = DMatrix::from_row_slice(n_samples, n_features, &mat_data);

    // Step 4: Compute SVD
    let svd = x_matrix.svd(true, true);

    // Extract V (right singular vectors - these are the principal components)
    let v_matrix = svd
        .v_t
        .ok_or("SVD failed to compute V matrix")?
        .transpose();

    // Extract singular values
    let singular_values: Vec<f64> = svd.singular_values.iter().copied().collect();

    // Step 5: Compute explained variance
    let total_variance: f64 = singular_values
        .iter()
        .take(k)
        .map(|&s| s * s / (n_samples - 1) as f64)
        .sum();

    let mut explained_variance = Vec::new();
    let mut explained_variance_ratio = Vec::new();

    for i in 0..k {
        let variance = singular_values[i] * singular_values[i] / (n_samples - 1) as f64;
        explained_variance.push(variance);
        explained_variance_ratio.push(variance / total_variance);
    }

    // Step 6: Extract principal components (columns of V)
    let mut components = Vec::new();
    for i in 0..k {
        let mut component = Vec::new();
        for j in 0..n_features {
            component.push(v_matrix[(j, i)]);
        }
        components.push(component);
    }

    Ok(PCAResult {
        components,
        explained_variance,
        explained_variance_ratio,
        mean,
        singular_values: singular_values.into_iter().take(k).collect(),
    })
}

/// Transform data to PCA space
///
/// Projects data onto the principal components.
///
/// # Arguments
/// * `data` - Original data (same features as training data)
/// * `pca_result` - Result from PCA computation
///
/// # Returns
/// Transformed data in reduced dimensionality space
pub fn transform(data: &[Vec<f64>], pca_result: &PCAResult) -> Result<Vec<Vec<f64>>, String> {
    if data.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    let n_features = pca_result.mean.len();
    let n_components = pca_result.components.len();

    let mut transformed = Vec::new();

    for sample in data {
        if sample.len() != n_features {
            return Err(format!(
                "Sample has {} features, expected {}",
                sample.len(),
                n_features
            ));
        }

        // Center the sample
        let centered: Vec<f64> = sample
            .iter()
            .zip(&pca_result.mean)
            .map(|(&x, &m)| x - m)
            .collect();

        // Project onto each principal component
        let mut transformed_sample = Vec::new();
        for component in &pca_result.components {
            let projection: f64 = centered
                .iter()
                .zip(component)
                .map(|(&x, &c)| x * c)
                .sum();
            transformed_sample.push(projection);
        }

        transformed.push(transformed_sample);
    }

    Ok(transformed)
}

/// Inverse transform from PCA space back to original space
///
/// Reconstructs original data from principal component representation.
///
/// # Arguments
/// * `transformed_data` - Data in PCA space
/// * `pca_result` - Result from PCA computation
///
/// # Returns
/// Reconstructed data in original feature space
pub fn inverse_transform(
    transformed_data: &[Vec<f64>],
    pca_result: &PCAResult,
) -> Result<Vec<Vec<f64>>, String> {
    if transformed_data.is_empty() {
        return Err("Data cannot be empty".to_string());
    }

    let n_features = pca_result.mean.len();
    let n_components = pca_result.components.len();

    let mut reconstructed = Vec::new();

    for sample in transformed_data {
        if sample.len() != n_components {
            return Err(format!(
                "Sample has {} components, expected {}",
                sample.len(),
                n_components
            ));
        }

        // Reconstruct by weighted sum of principal components
        let mut reconstructed_sample = vec![0.0; n_features];

        for (i, component) in pca_result.components.iter().enumerate() {
            let weight = sample[i];
            for (j, &c) in component.iter().enumerate() {
                reconstructed_sample[j] += weight * c;
            }
        }

        // Add back the mean
        for (j, &m) in pca_result.mean.iter().enumerate() {
            reconstructed_sample[j] += m;
        }

        reconstructed.push(reconstructed_sample);
    }

    Ok(reconstructed)
}

/// Compute cumulative explained variance
///
/// Useful for determining how many components to retain.
///
/// # Returns
/// Vector where cumulative_variance[i] = sum of explained_variance_ratio[0..=i]
pub fn cumulative_explained_variance(pca_result: &PCAResult) -> Vec<f64> {
    let mut cumulative = Vec::new();
    let mut sum = 0.0;

    for &ratio in &pca_result.explained_variance_ratio {
        sum += ratio;
        cumulative.push(sum);
    }

    cumulative
}

/// Determine number of components needed to explain target variance
///
/// # Arguments
/// * `pca_result` - Result from PCA computation
/// * `target_variance` - Target cumulative explained variance (e.g., 0.95 for 95%)
///
/// # Returns
/// Number of components needed to reach target variance
pub fn components_for_variance(pca_result: &PCAResult, target_variance: f64) -> usize {
    let cumulative = cumulative_explained_variance(pca_result);

    for (i, &var) in cumulative.iter().enumerate() {
        if var >= target_variance {
            return i + 1;
        }
    }

    pca_result.components.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pca_basic() {
        // Simple 2D dataset
        let data = vec![
            vec![2.5, 2.4],
            vec![0.5, 0.7],
            vec![2.2, 2.9],
            vec![1.9, 2.2],
            vec![3.1, 3.0],
        ];

        let result = pca(&data, None).unwrap();

        // Should have 2 components (min of samples and features)
        assert_eq!(result.components.len(), 2);

        // Each component should have 2 features
        assert_eq!(result.components[0].len(), 2);

        // Explained variance ratios should sum to approximately 1.0
        let total_ratio: f64 = result.explained_variance_ratio.iter().sum();
        assert!((total_ratio - 1.0).abs() < 1e-10);

        // First component should explain more variance than second
        assert!(result.explained_variance[0] > result.explained_variance[1]);
    }

    #[test]
    fn test_pca_transform() {
        let data = vec![
            vec![2.5, 2.4],
            vec![0.5, 0.7],
            vec![2.2, 2.9],
            vec![1.9, 2.2],
            vec![3.1, 3.0],
        ];

        let result = pca(&data, Some(1)).unwrap();

        // Transform to 1D
        let transformed = transform(&data, &result).unwrap();

        assert_eq!(transformed.len(), data.len());
        assert_eq!(transformed[0].len(), 1);
    }

    #[test]
    fn test_pca_reconstruction() {
        let data = vec![
            vec![2.5, 2.4],
            vec![0.5, 0.7],
            vec![2.2, 2.9],
        ];

        // Use all components for perfect reconstruction
        let result = pca(&data, None).unwrap();
        let transformed = transform(&data, &result).unwrap();
        let reconstructed = inverse_transform(&transformed, &result).unwrap();

        // Reconstruction should be very close to original
        for (original, recon) in data.iter().zip(reconstructed.iter()) {
            for (o, r) in original.iter().zip(recon.iter()) {
                assert!((o - r).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_cumulative_variance() {
        let data = vec![
            vec![2.5, 2.4, 1.0],
            vec![0.5, 0.7, 0.2],
            vec![2.2, 2.9, 1.1],
            vec![1.9, 2.2, 0.9],
            vec![3.1, 3.0, 1.3],
        ];

        let result = pca(&data, None).unwrap();
        let cumulative = cumulative_explained_variance(&result);

        // Cumulative should be monotonically increasing
        for i in 1..cumulative.len() {
            assert!(cumulative[i] >= cumulative[i - 1]);
        }

        // Last value should be approximately 1.0
        assert!((cumulative.last().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_components_for_variance() {
        let data = vec![
            vec![2.5, 2.4, 1.0],
            vec![0.5, 0.7, 0.2],
            vec![2.2, 2.9, 1.1],
            vec![1.9, 2.2, 0.9],
            vec![3.1, 3.0, 1.3],
        ];

        let result = pca(&data, None).unwrap();

        // Should need at least 1 component for any variance
        assert!(components_for_variance(&result, 0.5) >= 1);

        // Should need all components for 100% variance
        assert_eq!(
            components_for_variance(&result, 1.0),
            result.components.len()
        );
    }
}
