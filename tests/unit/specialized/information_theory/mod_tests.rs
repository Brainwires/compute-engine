// Unit tests for specialized::information_theory::mod
use computational_engine::specialized::information_theory::mod::*;

use super::*;

    #[test]
    fn test_shannon_entropy() {
        let data = vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0];
        let result = shannon_entropy(EntropyRequest {
            data,
            base: Some(2.0),
            entropy_type: Some("shannon".to_string()),
            alpha: None,
        })
        .unwrap();

        assert!(result.entropy > 0.0);
        assert!(result.normalized_entropy <= 1.0);
    }

    #[test]
    fn test_mutual_information() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0]; // y = 2x, should have high MI

        let result = mutual_information(MutualInfoRequest {
            x,
            y,
            bins: Some(3),
        })
        .unwrap();

        assert!(result.mutual_information > 0.0);
    }
