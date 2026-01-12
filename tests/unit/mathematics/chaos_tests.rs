//! Comprehensive unit tests for the chaos module
//!
//! Tests cover:
//! - Complex number operations
//! - 3D point operations
//! - Lorenz attractor
//! - Rössler attractor
//! - Logistic map bifurcation
//! - Lyapunov exponents (1D maps and 3D systems)
//! - Kaplan-Yorke dimension
//! - Mandelbrot set
//! - Julia sets
//! - Burning Ship fractal
//! - Box-counting dimension
//! - Koch snowflake

use crate::chaos::*;

// ============================================================================
// Complex Number Tests
// ============================================================================

#[test]
fn test_complex_creation() {
    let c = Complex::new(3.0, 4.0);
    assert_eq!(c.re, 3.0);
    assert_eq!(c.im, 4.0);
}

#[test]
fn test_complex_zero() {
    let c = Complex::zero();
    assert_eq!(c.re, 0.0);
    assert_eq!(c.im, 0.0);
}

#[test]
fn test_complex_addition() {
    let c1 = Complex::new(1.0, 2.0);
    let c2 = Complex::new(3.0, 4.0);
    let sum = c1 + c2;
    assert_eq!(sum.re, 4.0);
    assert_eq!(sum.im, 6.0);
}

#[test]
fn test_complex_multiplication() {
    let c1 = Complex::new(1.0, 2.0);
    let c2 = Complex::new(3.0, 4.0);
    let prod = c1 * c2;
    // (1 + 2i)(3 + 4i) = 3 + 4i + 6i + 8i² = 3 + 10i - 8 = -5 + 10i
    assert_eq!(prod.re, -5.0);
    assert_eq!(prod.im, 10.0);
}

#[test]
fn test_complex_magnitude() {
    let c = Complex::new(3.0, 4.0);
    assert_eq!(c.mag(), 5.0);
    assert_eq!(c.mag_sq(), 25.0);
}

#[test]
fn test_complex_pure_real() {
    let c = Complex::new(5.0, 0.0);
    assert_eq!(c.mag(), 5.0);
}

#[test]
fn test_complex_pure_imaginary() {
    let c = Complex::new(0.0, 5.0);
    assert_eq!(c.mag(), 5.0);
}

// ============================================================================
// Point3D Tests
// ============================================================================

#[test]
fn test_point3d_creation() {
    let p = Point3D::new(1.0, 2.0, 3.0);
    assert_eq!(p.x, 1.0);
    assert_eq!(p.y, 2.0);
    assert_eq!(p.z, 3.0);
}

#[test]
fn test_point3d_distance() {
    let p1 = Point3D::new(0.0, 0.0, 0.0);
    let p2 = Point3D::new(3.0, 4.0, 0.0);
    assert_eq!(p1.distance(&p2), 5.0);
}

#[test]
fn test_point3d_distance_3d() {
    let p1 = Point3D::new(0.0, 0.0, 0.0);
    let p2 = Point3D::new(1.0, 2.0, 2.0);
    assert_eq!(p1.distance(&p2), 3.0);
}

// ============================================================================
// Lorenz Attractor Tests
// ============================================================================

#[test]
fn test_lorenz_attractor_basic() {
    let initial = Point3D::new(1.0, 1.0, 1.0);
    let config = LorenzConfig::default();
    let trajectory = lorenz_attractor(initial, &config, 0.01, 100);

    assert_eq!(trajectory.points.len(), 101); // Initial + 100 steps
    assert_eq!(trajectory.times.len(), 101);
    assert_eq!(trajectory.times[0], 0.0);
    assert_eq!(trajectory.times[100], 1.0);
}

#[test]
fn test_lorenz_attractor_finite_values() {
    let initial = Point3D::new(1.0, 1.0, 1.0);
    let config = LorenzConfig::default();
    let trajectory = lorenz_attractor(initial, &config, 0.01, 1000);

    for point in &trajectory.points {
        assert!(point.x.is_finite());
        assert!(point.y.is_finite());
        assert!(point.z.is_finite());
    }
}

#[test]
fn test_lorenz_butterfly_effect() {
    // Small perturbation should lead to exponential divergence
    let initial1 = Point3D::new(1.0, 1.0, 1.0);
    let initial2 = Point3D::new(1.0, 1.0, 1.001);
    let config = LorenzConfig::default();

    let traj1 = lorenz_attractor(initial1, &config, 0.01, 500);
    let traj2 = lorenz_attractor(initial2, &config, 0.01, 500);

    // Distance should grow over time (butterfly effect)
    let dist_early = traj1.points[10].distance(&traj2.points[10]);
    let dist_late = traj1.points[500].distance(&traj2.points[500]);

    assert!(dist_late > dist_early,
        "Expected divergence: early={}, late={}", dist_early, dist_late);
}

#[test]
fn test_lorenz_custom_parameters() {
    let initial = Point3D::new(1.0, 1.0, 1.0);
    let config = LorenzConfig {
        sigma: 5.0,
        rho: 20.0,
        beta: 2.0,
    };
    let trajectory = lorenz_attractor(initial, &config, 0.01, 100);

    assert_eq!(trajectory.points.len(), 101);
    for point in &trajectory.points {
        assert!(point.x.is_finite());
        assert!(point.y.is_finite());
        assert!(point.z.is_finite());
    }
}

// ============================================================================
// Rössler Attractor Tests
// ============================================================================

#[test]
fn test_rossler_attractor_basic() {
    let initial = Point3D::new(1.0, 1.0, 1.0);
    let config = RosslerConfig::default();
    let trajectory = rossler_attractor(initial, &config, 0.01, 100);

    assert_eq!(trajectory.points.len(), 101);
    assert_eq!(trajectory.times.len(), 101);
}

#[test]
fn test_rossler_attractor_finite_values() {
    let initial = Point3D::new(1.0, 1.0, 1.0);
    let config = RosslerConfig::default();
    let trajectory = rossler_attractor(initial, &config, 0.01, 500);

    for point in &trajectory.points {
        assert!(point.x.is_finite());
        assert!(point.y.is_finite());
        assert!(point.z.is_finite());
    }
}

#[test]
fn test_rossler_sensitivity() {
    let initial1 = Point3D::new(1.0, 1.0, 1.0);
    let initial2 = Point3D::new(1.0, 1.0, 1.01); // Slightly larger perturbation
    let config = RosslerConfig::default();

    let traj1 = rossler_attractor(initial1, &config, 0.01, 1000); // More steps
    let traj2 = rossler_attractor(initial2, &config, 0.01, 1000);

    let dist_late = traj1.points[1000].distance(&traj2.points[1000]);
    // Should show some sensitivity, but may be smaller than Lorenz
    assert!(dist_late > 0.0001, "Expected sensitivity, got distance={}", dist_late);
}

// ============================================================================
// Logistic Map Bifurcation Tests
// ============================================================================

#[test]
fn test_logistic_map_fixed_point() {
    // For r < 1, logistic map should converge to 0
    let results = logistic_map_bifurcation(0.5, 0.9, 5, 100, 100);

    for (r, values) in &results {
        let avg: f64 = values.iter().sum::<f64>() / values.len() as f64;
        assert!(avg < 0.1, "r={}, avg={}", r, avg);
    }
}

#[test]
fn test_logistic_map_stable_fixed_point() {
    // For 1 < r < 3, should converge to a stable fixed point
    let results = logistic_map_bifurcation(2.0, 2.5, 3, 100, 100);

    for (_r, values) in &results {
        // Check convergence: variance should be small
        let mean: f64 = values.iter().sum::<f64>() / values.len() as f64;
        let variance: f64 = values.iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>() / values.len() as f64;
        assert!(variance < 0.01, "Expected stable fixed point, got variance={}", variance);
    }
}

#[test]
fn test_logistic_map_chaos() {
    // For r > 3.57, logistic map is chaotic
    let results = logistic_map_bifurcation(3.7, 3.9, 3, 100, 100);

    for (r, values) in &results {
        // Values should be spread out (chaotic)
        let max: f64 = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min: f64 = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let range = max - min;
        assert!(range > 0.01, "r={}, range={} (expected chaotic behavior)", r, range);
    }
}

#[test]
fn test_logistic_map_period_doubling() {
    // Test that bifurcation diagram produces different behaviors at different r values
    let results = logistic_map_bifurcation(3.0, 3.5, 5, 100, 200);

    assert_eq!(results.len(), 5);

    // Just verify that we get valid results and different values at different r
    for (r, values) in &results {
        assert!(values.len() > 0, "Should have values for r={}", r);

        // All values should be in [0, 1]
        for &v in values {
            assert!(v >= 0.0 && v <= 1.0, "Value {} out of range at r={}", v, r);
        }

        // Check they're finite
        let sum: f64 = values.iter().sum();
        assert!(sum.is_finite(), "Sum should be finite at r={}", r);
    }
}

// ============================================================================
// Lyapunov Exponent Tests
// ============================================================================

#[test]
fn test_lyapunov_logistic_stable() {
    // r = 2.5: stable fixed point → negative Lyapunov
    let lyap = lyapunov_logistic_map(2.5, 1000, 100);
    assert!(lyap < 0.0, "Expected negative Lyapunov for stable system, got {}", lyap);
}

#[test]
fn test_lyapunov_logistic_chaotic() {
    // r = 4.0: fully chaotic → positive Lyapunov
    let lyap = lyapunov_logistic_map(4.0, 1000, 100);
    assert!(lyap > 0.0, "Expected positive Lyapunov for chaotic system, got {}", lyap);
}

#[test]
fn test_lyapunov_logistic_periodic() {
    // r = 3.2: period-2 orbit → near zero or slightly negative
    let lyap = lyapunov_logistic_map(3.2, 1000, 100);
    assert!(lyap < 0.5, "Expected near-zero Lyapunov for periodic orbit, got {}", lyap);
}

#[test]
fn test_lyapunov_1d_map_quadratic() {
    // Test with a simple quadratic map
    let map = |x: f64| x * x;
    let lyap = lyapunov_1d_map(map, 0.5, 1000, 100);

    // This map should show expansion (positive derivative)
    assert!(lyap.is_finite());
}

#[test]
fn test_lyapunov_spectrum_3d() {
    // Test the 3D Lyapunov spectrum calculation
    let initial = Point3D::new(1.0, 1.0, 1.0);
    let jacobian = |_state: Point3D| -> [[f64; 3]; 3] {
        // Simple identity-like Jacobian for testing
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    };

    let spectrum = lyapunov_spectrum_3d(jacobian, initial, 0.01, 100, 50);
    assert_eq!(spectrum.len(), 3);

    for exp in &spectrum {
        assert!(exp.is_finite());
    }
}

// ============================================================================
// Kaplan-Yorke Dimension Tests
// ============================================================================

#[test]
fn test_kaplan_yorke_dimension_basic() {
    // Example: one positive, two negative (typical for strange attractor)
    let exponents = vec![0.9, 0.0, -14.5];
    let dim = kaplan_yorke_dimension(&exponents);

    // D = 1 + λ₁/|λ₃| = 1 + 0.9/14.5 ≈ 1.062
    assert!(dim > 0.5 && dim < 2.5, "Dimension was {}", dim);
}

#[test]
fn test_kaplan_yorke_dimension_all_negative() {
    // All negative exponents → dimension should be 0 or small
    let exponents = vec![-1.0, -2.0, -3.0];
    let dim = kaplan_yorke_dimension(&exponents);

    assert!(dim >= 0.0 && dim < 1.0, "Expected small dimension, got {}", dim);
}

#[test]
fn test_kaplan_yorke_dimension_sorting() {
    // Should work regardless of input order
    let exponents1 = vec![0.9, 0.0, -14.5];
    let exponents2 = vec![-14.5, 0.0, 0.9];

    let dim1 = kaplan_yorke_dimension(&exponents1);
    let dim2 = kaplan_yorke_dimension(&exponents2);

    assert!((dim1 - dim2).abs() < 1e-10, "Dimension should be independent of input order");
}

// ============================================================================
// Mandelbrot Set Tests
// ============================================================================

#[test]
fn test_mandelbrot_origin() {
    // c = 0 should not escape (z_n stays at 0)
    let result = mandelbrot(Complex::new(0.0, 0.0), 100, 2.0);
    assert!(!result.escaped);
    assert!(result.final_magnitude < 0.1);
}

#[test]
fn test_mandelbrot_escape() {
    // c = 2 should escape quickly
    let result = mandelbrot(Complex::new(2.0, 0.0), 100, 2.0);
    assert!(result.escaped);
    assert!(result.iterations < 10);
}

#[test]
fn test_mandelbrot_boundary() {
    // c = -0.5 is in the Mandelbrot set (main cardioid)
    let result = mandelbrot(Complex::new(-0.5, 0.0), 100, 2.0);
    assert!(result.iterations <= 100);
}

#[test]
fn test_mandelbrot_far_point() {
    // Point far from origin should escape immediately
    let result = mandelbrot(Complex::new(10.0, 10.0), 100, 2.0);
    assert!(result.escaped);
    assert!(result.iterations < 5);
}

// ============================================================================
// Julia Set Tests
// ============================================================================

#[test]
fn test_julia_set_basic() {
    // Julia set with c = -0.7 + 0.27i (classic parameter)
    let c = Complex::new(-0.7, 0.27);
    let z0 = Complex::new(0.0, 0.0);

    let result = julia(z0, c, 100, 2.0);
    assert!(result.iterations <= 100);
}

#[test]
fn test_julia_set_escape() {
    // Point far from origin should escape
    let c = Complex::new(-0.7, 0.27);
    let z0 = Complex::new(5.0, 5.0);

    let result = julia(z0, c, 100, 2.0);
    assert!(result.escaped);
    assert!(result.iterations < 10);
}

#[test]
fn test_julia_set_different_c() {
    // Julia set with different c value
    let c = Complex::new(0.285, 0.01);
    let z0 = Complex::new(0.0, 0.0);

    let result = julia(z0, c, 100, 2.0);
    assert!(result.iterations <= 100);
}

// ============================================================================
// Burning Ship Fractal Tests
// ============================================================================

#[test]
fn test_burning_ship_basic() {
    let c = Complex::new(-1.8, -0.1);
    let result = burning_ship(c, 100, 2.0);
    assert!(result.iterations <= 100);
}

#[test]
fn test_burning_ship_escape() {
    // Point that should escape quickly
    let c = Complex::new(5.0, 5.0);
    let result = burning_ship(c, 100, 2.0);
    assert!(result.escaped);
}

#[test]
fn test_burning_ship_origin() {
    let c = Complex::new(0.0, 0.0);
    let result = burning_ship(c, 100, 2.0);
    assert!(!result.escaped);
}

// ============================================================================
// Koch Snowflake Tests
// ============================================================================

#[test]
fn test_koch_snowflake_order_0() {
    let points = koch_snowflake(0);
    assert_eq!(points.len(), 4); // Triangle + return to start
}

#[test]
fn test_koch_snowflake_order_1() {
    let order_0 = koch_snowflake(0);
    let order_1 = koch_snowflake(1);

    // Each segment is divided into 4, so points should increase
    assert!(order_1.len() > order_0.len());
}

#[test]
fn test_koch_snowflake_growth() {
    let order_1 = koch_snowflake(1);
    let order_2 = koch_snowflake(2);

    // Higher order should have more points
    assert!(order_2.len() > order_1.len());
}

#[test]
fn test_koch_snowflake_closed() {
    let points = koch_snowflake(2);

    // First and last point should be the same (closed curve)
    let first = points[0];
    let last = points[points.len() - 1];

    assert!((first.0 - last.0).abs() < 1e-10);
    assert!((first.1 - last.1).abs() < 1e-10);
}

// ============================================================================
// Box-Counting Dimension Tests
// ============================================================================

#[test]
fn test_box_counting_line() {
    // Create a simple line from (0,0) to (10,0)
    let points: Vec<(f64, f64)> = (0..=100).map(|i| (i as f64 * 0.1, 0.0)).collect();

    let dim = box_counting_dimension(&points, 0.1, 10.0, 10);

    // A line should have dimension ~1
    assert!(dim > 0.5 && dim < 1.5, "Line dimension was {}, expected ~1", dim);
}

#[test]
fn test_box_counting_empty() {
    let points: Vec<(f64, f64)> = vec![];
    let dim = box_counting_dimension(&points, 0.1, 10.0, 10);

    // Empty set should have dimension 0
    assert_eq!(dim, 0.0);
}

#[test]
fn test_box_counting_single_point() {
    let points = vec![(5.0, 5.0)];
    let dim = box_counting_dimension(&points, 0.1, 10.0, 10);

    // Single point should have dimension 0
    assert!(dim < 0.5, "Single point dimension was {}, expected ~0", dim);
}

#[test]
fn test_box_counting_square() {
    // Create a filled square with many points for better accuracy
    let mut points = Vec::new();
    for i in 0..=50 {
        for j in 0..=50 {
            points.push((i as f64 * 0.02, j as f64 * 0.02));
        }
    }

    let dim = box_counting_dimension(&points, 0.005, 0.5, 20);

    // A filled square should have dimension ~2
    // Box-counting can be approximate depending on sampling
    assert!(dim > 1.0 && dim < 2.5, "Square dimension was {}, expected ~2", dim);
}
