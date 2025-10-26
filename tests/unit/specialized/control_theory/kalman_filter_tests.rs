// Unit tests for specialized::control_theory::kalman_filter
use computational_engine::specialized::control_theory::kalman_filter::*;

use super::*;

    #[test]
    fn test_kalman_creation() {
        let kf = KalmanFilter::new(
            vec![0.0],
            vec![vec![1.0]],
            vec![vec![0.1]],
            vec![vec![1.0]],
        );
        assert_eq!(kf.x.len(), 1);
    }
