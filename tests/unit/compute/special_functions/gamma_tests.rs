// Unit tests for mathematics::special_functions::gamma
use computational_engine::compute::special_functions::gamma::*;

use super::*;

    // ==================== Gamma Function Tests ====================

    #[test]
    fn test_gamma_integer() {
        let result = gamma_function(GammaRequest {
            x: 5.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(5) = 4! = 24
        assert!((result.value - 24.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_half() {
        let result = gamma_function(GammaRequest {
            x: 0.5,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(0.5) = √π ≈ 1.77245
        assert!((result.value - 1.77245).abs() < 0.01);
    }

    #[test]
    fn test_gamma_one() {
        let result = gamma_function(GammaRequest {
            x: 1.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(1) = 0! = 1
        assert!((result.value - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_two() {
        let result = gamma_function(GammaRequest {
            x: 2.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(2) = 1! = 1
        assert!((result.value - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_three() {
        let result = gamma_function(GammaRequest {
            x: 3.0,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(3) = 2! = 2
        assert!((result.value - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_seven_halves() {
        let result = gamma_function(GammaRequest {
            x: 3.5,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        // Γ(3.5) = 2.5 * 1.5 * 0.5 * Γ(0.5) ≈ 3.32335
        assert!((result.value - 3.32335).abs() < 0.01);
    }

    #[test]
    fn test_log_gamma() {
        let result = gamma_function(GammaRequest {
            x: 5.0,
            function: "log_gamma".to_string(),
            y: None,
        })
        .unwrap();
        // ln(Γ(5)) = ln(24) ≈ 3.17805
        assert!((result.value - 3.17805).abs() < 0.01);
    }

    #[test]
    fn test_digamma() {
        let result = gamma_function(GammaRequest {
            x: 1.0,
            function: "digamma".to_string(),
            y: None,
        })
        .unwrap();
        // ψ(1) ≈ -0.5772 (negative Euler-Mascheroni constant)
        assert!((result.value - (-0.5772)).abs() < 0.1);
    }

    #[test]
    fn test_beta_function() {
        let result = gamma_function(GammaRequest {
            x: 2.0,
            function: "beta".to_string(),
            y: Some(3.0),
        })
        .unwrap();
        // B(2,3) = Γ(2)Γ(3)/Γ(5) = 1*2/24 = 1/12 ≈ 0.0833
        assert!((result.value - 0.0833).abs() < 0.01);
    }

    #[test]
    fn test_beta_symmetric() {
        let result1 = gamma_function(GammaRequest {
            x: 3.0,
            function: "beta".to_string(),
            y: Some(5.0),
        })
        .unwrap();
        let result2 = gamma_function(GammaRequest {
            x: 5.0,
            function: "beta".to_string(),
            y: Some(3.0),
        })
        .unwrap();
        // B(3,5) = B(5,3) (symmetry)
        assert!((result1.value - result2.value).abs() < 1e-6);
    }

    #[test]
    fn test_beta_halves() {
        let result = gamma_function(GammaRequest {
            x: 0.5,
            function: "beta".to_string(),
            y: Some(0.5),
        })
        .unwrap();
        // B(0.5,0.5) = π
        assert!((result.value - PI).abs() < 0.01);
    }

    // ==================== Additional Comprehensive Tests ====================

    #[test]
    fn test_gamma_reflection_formula() {
        // Γ(z)Γ(1-z) = π/sin(πz)
        let z = 0.3;
        let gamma_z = gamma_function(GammaRequest {
            x: z,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        let gamma_1_minus_z = gamma_function(GammaRequest {
            x: 1.0 - z,
            function: "gamma".to_string(),
            y: None,
        })
        .unwrap();
        let product = gamma_z.value * gamma_1_minus_z.value;
        let expected = PI / (PI * z).sin();
        assert!((product - expected).abs() < 0.1);
    }

    #[test]
    fn test_beta_integral_relation() {
        // B(a,b) = ∫[0 to 1] t^(a-1) * (1-t)^(b-1) dt
        // For a=2, b=2: B(2,2) = 1/6
        let result = gamma_function(GammaRequest {
            x: 2.0,
            function: "beta".to_string(),
            y: Some(2.0),
        })
        .unwrap();
        assert!((result.value - 1.0 / 6.0).abs() < 0.01);
    }

    #[test]
    fn test_gamma_factorial_sequence() {
        // Γ(n+1) = n! for positive integers
        fn factorial(n: usize) -> f64 {
            (1..=n).map(|i| i as f64).product()
        }

        for n in 1..8 {
            let gamma_n_plus_1 = gamma_function(GammaRequest {
                x: (n + 1) as f64,
                function: "gamma".to_string(),
                y: None,
            })
            .unwrap();
            let factorial_n = factorial(n);
            assert!((gamma_n_plus_1.value - factorial_n).abs() < 0.01);
        }
    }
