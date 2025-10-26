// Unit tests for specialized::machine_learning::clustering
use crate::specialized::machine_learning::*;

use super::*;

// Import private functions from clustering module
use super::clustering::euclidean_distance;

    #[test]
    fn test_kmeans_basic() {
        // Two clear clusters
        let data = vec![
            vec![1.0, 2.0],
            vec![1.5, 1.8],
            vec![5.0, 8.0],
            vec![8.0, 8.0],
            vec![1.0, 0.6],
            vec![9.0, 11.0],
        ];

        let result = kmeans(&data, 2, Some(100), Some(1e-4), Some("kmeans++")).unwrap();

        assert_eq!(result.centroids.len(), 2);
        assert_eq!(result.labels.len(), 6);
        assert!(result.converged);

        // First 2 points should be in same cluster
        assert_eq!(result.labels[0], result.labels[1]);

        // Last 2 points should be in same cluster
        assert_eq!(result.labels[3], result.labels[5]);
    }

    #[test]
    fn test_kmeans_convergence() {
        let data = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.1],
            vec![10.0, 10.0],
            vec![10.1, 10.1],
        ];

        let result = kmeans(&data, 2, Some(100), Some(1e-4), Some("random")).unwrap();

        assert!(result.converged);
        assert!(result.n_iterations < 100);
    }

    #[test]
    fn test_inertia_decreases() {
        let data = vec![
            vec![1.0, 2.0],
            vec![1.5, 1.8],
            vec![5.0, 8.0],
            vec![8.0, 8.0],
        ];

        let result1 = kmeans(&data, 1, Some(100), Some(1e-4), Some("random")).unwrap();
        let result2 = kmeans(&data, 2, Some(100), Some(1e-4), Some("random")).unwrap();

        // More clusters should have lower inertia
        assert!(result2.inertia < result1.inertia);
    }

    #[test]
    fn test_silhouette_score() {
        let data = vec![
            vec![0.0, 0.0],
            vec![0.1, 0.1],
            vec![10.0, 10.0],
            vec![10.1, 10.1],
        ];

        let labels = vec![0, 0, 1, 1]; // Perfect clustering

        let score = silhouette_score(&data, &labels).unwrap();

        // Should be close to 1.0 for well-separated clusters
        assert!(score > 0.8);
    }

    #[test]
    fn test_euclidean_distance() {
        let a = vec![0.0, 0.0];
        let b = vec![3.0, 4.0];

        let dist = euclidean_distance(&a, &b);
        assert!((dist - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_cluster_handling() {
        // Test that empty clusters are handled gracefully
        let data = vec![
            vec![1.0, 1.0],
            vec![1.1, 1.1],
            vec![1.2, 1.2],
        ];

        // Requesting many clusters might create empty ones
        let result = kmeans(&data, 3, Some(100), Some(1e-4), Some("random"));

        assert!(result.is_ok());
    }
