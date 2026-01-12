// Unit tests for specialized::machine_learning::dimensionality_reduction
use crate::ml::*;

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
