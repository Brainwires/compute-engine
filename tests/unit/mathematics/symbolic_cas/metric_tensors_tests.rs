// Unit tests for mathematics::symbolic_cas::metric_tensors
use computational_engine::analyze::symbolic::metric_tensors::*;

use super::*;

    #[test]
    fn test_minkowski_metric() {
        let metric = minkowski_metric(true).unwrap();
        assert_eq!(metric.rows(), 4);
        assert_eq!(metric.cols(), 4);

        // Check signature
        assert_eq!(metric.get(0, 0), Some(&Expr::num(1)));
        assert_eq!(metric.get(1, 1), Some(&Expr::num(-1)));
    }

    #[test]
    fn test_euclidean_metric() {
        let metric = euclidean_metric(3);
        assert_eq!(metric.rows(), 3);

        // Should be identity
        assert_eq!(metric.get(0, 0), Some(&Expr::num(1)));
        assert_eq!(metric.get(0, 1), Some(&Expr::num(0)));
    }

    #[test]
    fn test_schwarzschild_metric() {
        let metric = schwarzschild_metric().unwrap();
        assert_eq!(metric.rows(), 4);
        println!("Schwarzschild metric:");
        println!("{}", metric);
    }

    #[test]
    fn test_kerr_metric() {
        let metric = kerr_metric_2d().unwrap();
        println!("Kerr metric (2D):");
        println!("{}", metric);
    }

    #[test]
    fn test_flrw_metric() {
        let metric = flrw_metric_2d().unwrap();
        println!("FLRW metric (2D):");
        println!("{}", metric);
    }

    #[test]
    fn test_ads_metric() {
        let metric = ads_metric_2d().unwrap();
        println!("AdS metric (2D):");
        println!("{}", metric);
    }
