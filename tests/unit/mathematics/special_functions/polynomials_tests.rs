// Unit tests for mathematics::special_functions::polynomials
use computational_engine::compute::special_functions::polynomials::*;

use super::*;

    // ==================== Legendre Polynomial Tests ====================

    #[test]
    fn test_legendre_p0() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "legendre".to_string(),
            n: 0,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // P_0(x) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_legendre_p1() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "legendre".to_string(),
            n: 1,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // P_1(x) = x
        assert!((result.value - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_legendre_p2() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "legendre".to_string(),
            n: 2,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // P_2(0.5) = (3*0.5^2 - 1)/2 = -0.125
        assert!((result.value - (-0.125)).abs() < 1e-6);
    }

    #[test]
    fn test_legendre_p3() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "legendre".to_string(),
            n: 3,
            x: 1.0,
            alpha: None,
        }).unwrap();
        // P_3(1) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_legendre_p4_zero() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "legendre".to_string(),
            n: 4,
            x: 0.0,
            alpha: None,
        }).unwrap();
        // P_4(0) = 3/8
        assert!((result.value - 0.375).abs() < 1e-6);
    }

    #[test]
    fn test_legendre_orthogonality() {
        // P_n(1) = 1 for all n
        for n in 0..5 {
            let result = orthogonal_polynomial(PolynomialRequest {
                polynomial_type: "legendre".to_string(),
                n,
                x: 1.0,
                alpha: None,
            }).unwrap();
            assert!((result.value - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_legendre_normalization() {
        // P_n(1) = 1, P_n(-1) = (-1)^n
        for n in 0..6 {
            let p_at_1 = orthogonal_polynomial(PolynomialRequest {
                polynomial_type: "legendre".to_string(),
                n,
                x: 1.0,
                alpha: None,
            }).unwrap();
            assert!((p_at_1.value - 1.0).abs() < 1e-6);

            let p_at_minus_1 = orthogonal_polynomial(PolynomialRequest {
                polynomial_type: "legendre".to_string(),
                n,
                x: -1.0,
                alpha: None,
            }).unwrap();
            let expected = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!((p_at_minus_1.value - expected).abs() < 1e-6);
        }
    }

    // ==================== Hermite Polynomial Tests ====================

    #[test]
    fn test_hermite_h0() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "hermite".to_string(),
            n: 0,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // H_0(x) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_hermite_h1() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "hermite".to_string(),
            n: 1,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // H_1(x) = 2x
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_hermite_h2() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "hermite".to_string(),
            n: 2,
            x: 1.0,
            alpha: None,
        }).unwrap();
        // H_2(1) = 4*1^2 - 2 = 2
        assert!((result.value - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_hermite_h3() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "hermite".to_string(),
            n: 3,
            x: 0.0,
            alpha: None,
        }).unwrap();
        // H_3(0) = 0 (odd function)
        assert!(result.value.abs() < 1e-6);
    }

    #[test]
    fn test_hermite_h4() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "hermite".to_string(),
            n: 4,
            x: 1.0,
            alpha: None,
        }).unwrap();
        // H_4(1) = 16*1^4 - 48*1^2 + 12 = -20
        assert!((result.value - (-20.0)).abs() < 1e-6);
    }

    #[test]
    fn test_hermite_h5() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "hermite".to_string(),
            n: 5,
            x: 2.0,
            alpha: None,
        }).unwrap();
        // H_5(2) = 32*2^5 - 160*2^3 + 120*2 = 1024 - 1280 + 240 = -16
        assert!((result.value - (-16.0)).abs() < 1e-3);
    }

    #[test]
    fn test_hermite_physicist_normalization() {
        // For physicist's Hermite: H_n(0) = 0 for odd n
        for n in 1..6 {
            if n % 2 == 1 {
                let result = orthogonal_polynomial(PolynomialRequest {
                    polynomial_type: "hermite".to_string(),
                    n,
                    x: 0.0,
                    alpha: None,
                }).unwrap();
                assert!(result.value.abs() < 1e-6);
            }
        }
    }

    // ==================== Laguerre Polynomial Tests ====================

    #[test]
    fn test_laguerre_l0() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "laguerre".to_string(),
            n: 0,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // L_0(x) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_laguerre_l1() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "laguerre".to_string(),
            n: 1,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // L_1(x) = 1 - x
        assert!((result.value - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_laguerre_l2() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "laguerre".to_string(),
            n: 2,
            x: 1.0,
            alpha: None,
        }).unwrap();
        // L_2(x) = (x^2 - 4x + 2)/2, L_2(1) = (1 - 4 + 2)/2 = -0.5
        assert!((result.value - (-0.5)).abs() < 1e-6);
    }

    #[test]
    fn test_laguerre_l3() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "laguerre".to_string(),
            n: 3,
            x: 0.0,
            alpha: None,
        }).unwrap();
        // L_3(0) = 1 (all Laguerre polynomials have L_n(0) = 1)
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_laguerre_l4() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "laguerre".to_string(),
            n: 4,
            x: 2.0,
            alpha: None,
        }).unwrap();
        // L_4(2) = 1/3 ≈ 0.333 (using recurrence relation)
        assert!((result.value - 0.333).abs() < 0.01);
    }

    // ==================== Chebyshev Polynomial Tests ====================

    #[test]
    fn test_chebyshev_t0() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "chebyshev".to_string(),
            n: 0,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // T_0(x) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_chebyshev_t1() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "chebyshev".to_string(),
            n: 1,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // T_1(x) = x
        assert!((result.value - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_chebyshev_t2() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "chebyshev".to_string(),
            n: 2,
            x: 0.5,
            alpha: None,
        }).unwrap();
        // T_2(x) = 2x^2 - 1, T_2(0.5) = 2*0.25 - 1 = -0.5
        assert!((result.value - (-0.5)).abs() < 1e-6);
    }

    #[test]
    fn test_chebyshev_t3() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "chebyshev".to_string(),
            n: 3,
            x: 1.0,
            alpha: None,
        }).unwrap();
        // T_3(1) = 1
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_chebyshev_t4() {
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "chebyshev".to_string(),
            n: 4,
            x: 0.0,
            alpha: None,
        }).unwrap();
        // T_4(0) = 1 (T_n(0) alternates 1, 0, -1, 0, 1, ...)
        assert!((result.value - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_chebyshev_bound() {
        // |T_n(x)| <= 1 for |x| <= 1
        let result = orthogonal_polynomial(PolynomialRequest {
            polynomial_type: "chebyshev".to_string(),
            n: 10,
            x: 0.8,
            alpha: None,
        }).unwrap();
        assert!(result.value.abs() <= 1.0 + 1e-6);
    }

    #[test]
    fn test_chebyshev_cosine_relation() {
        // T_n(cos(θ)) = cos(nθ)
        let theta = PI / 3.0;
        let x = theta.cos();
        for n in 0..5 {
            let t_n = orthogonal_polynomial(PolynomialRequest {
                polynomial_type: "chebyshev".to_string(),
                n,
                x,
                alpha: None,
            }).unwrap();
            let expected = (n as f64 * theta).cos();
            assert!((t_n.value - expected).abs() < 1e-6);
        }
    }
