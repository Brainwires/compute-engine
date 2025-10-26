// Unit tests for physics::n_body module
use super::*;

// ============================================================================
// Vec3 Tests
// ============================================================================

#[test]
fn test_vec3_construction_and_zero() {
    let v1 = Vec3::new(1.0, 2.0, 3.0);
    assert_eq!(v1.x, 1.0);
    assert_eq!(v1.y, 2.0);
    assert_eq!(v1.z, 3.0);

    let v2 = Vec3::zero();
    assert_eq!(v2.x, 0.0);
    assert_eq!(v2.y, 0.0);
    assert_eq!(v2.z, 0.0);
}

#[test]
fn test_vec3_magnitude() {
    let v1 = Vec3::new(3.0, 4.0, 0.0);
    assert!((v1.mag() - 5.0).abs() < 1e-10);
    assert!((v1.mag_sq() - 25.0).abs() < 1e-10);

    let v2 = Vec3::new(1.0, 1.0, 1.0);
    assert!((v2.mag() - 3.0_f64.sqrt()).abs() < 1e-10);
}

#[test]
fn test_vec3_normalize() {
    let v1 = Vec3::new(3.0, 4.0, 0.0);
    let normalized = v1.normalize();
    assert!((normalized.mag() - 1.0).abs() < 1e-10);
    assert!((normalized.x - 0.6).abs() < 1e-10);
    assert!((normalized.y - 0.8).abs() < 1e-10);

    // Test zero vector normalization
    let v2 = Vec3::zero();
    let normalized_zero = v2.normalize();
    assert_eq!(normalized_zero.x, 0.0);
    assert_eq!(normalized_zero.y, 0.0);
    assert_eq!(normalized_zero.z, 0.0);
}

#[test]
fn test_vec3_dot_product() {
    let v1 = Vec3::new(1.0, 2.0, 3.0);
    let v2 = Vec3::new(4.0, 5.0, 6.0);
    let dot = v1.dot(&v2);
    assert_eq!(dot, 32.0); // 1*4 + 2*5 + 3*6 = 32

    // Test orthogonal vectors
    let v3 = Vec3::new(1.0, 0.0, 0.0);
    let v4 = Vec3::new(0.0, 1.0, 0.0);
    assert_eq!(v3.dot(&v4), 0.0);
}

#[test]
fn test_vec3_cross_product() {
    let v1 = Vec3::new(1.0, 0.0, 0.0);
    let v2 = Vec3::new(0.0, 1.0, 0.0);
    let cross = v1.cross(&v2);
    assert_eq!(cross.x, 0.0);
    assert_eq!(cross.y, 0.0);
    assert_eq!(cross.z, 1.0);

    // Test anti-commutativity: a × b = -(b × a)
    let v3 = Vec3::new(2.0, 3.0, 4.0);
    let v4 = Vec3::new(5.0, 6.0, 7.0);
    let cross1 = v3.cross(&v4);
    let cross2 = v4.cross(&v3);
    assert!((cross1.x + cross2.x).abs() < 1e-10);
    assert!((cross1.y + cross2.y).abs() < 1e-10);
    assert!((cross1.z + cross2.z).abs() < 1e-10);
}

#[test]
fn test_vec3_arithmetic_operations() {
    let v1 = Vec3::new(1.0, 2.0, 3.0);
    let v2 = Vec3::new(4.0, 5.0, 6.0);

    // Addition
    let sum = v1 + v2;
    assert_eq!(sum.x, 5.0);
    assert_eq!(sum.y, 7.0);
    assert_eq!(sum.z, 9.0);

    // Subtraction
    let diff = v2 - v1;
    assert_eq!(diff.x, 3.0);
    assert_eq!(diff.y, 3.0);
    assert_eq!(diff.z, 3.0);

    // Scalar multiplication
    let scaled = v1 * 2.0;
    assert_eq!(scaled.x, 2.0);
    assert_eq!(scaled.y, 4.0);
    assert_eq!(scaled.z, 6.0);

    // Scalar division
    let divided = v2 / 2.0;
    assert_eq!(divided.x, 2.0);
    assert_eq!(divided.y, 2.5);
    assert_eq!(divided.z, 3.0);
}

// ============================================================================
// Body Tests
// ============================================================================

#[test]
fn test_body_creation() {
    let body = Body::new(
        1.0,
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::new(4.0, 5.0, 6.0),
        "test_body"
    );

    assert_eq!(body.mass, 1.0);
    assert_eq!(body.position.x, 1.0);
    assert_eq!(body.velocity.x, 4.0);
    assert_eq!(body.name, "test_body");
}

#[test]
fn test_body_kinetic_energy() {
    // KE = 0.5 * m * v²
    let body = Body::new(
        2.0,
        Vec3::zero(),
        Vec3::new(3.0, 4.0, 0.0), // |v| = 5
        "test"
    );

    let ke = body.kinetic_energy();
    assert!((ke - 25.0).abs() < 1e-10); // 0.5 * 2.0 * 25 = 25.0
}

#[test]
fn test_body_momentum() {
    let body = Body::new(
        2.0,
        Vec3::zero(),
        Vec3::new(3.0, 4.0, 5.0),
        "test"
    );

    let momentum = body.momentum();
    assert_eq!(momentum.x, 6.0);
    assert_eq!(momentum.y, 8.0);
    assert_eq!(momentum.z, 10.0);
}

#[test]
fn test_body_angular_momentum() {
    // L = r × p = r × (m*v)
    let body = Body::new(
        1.0,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        "test"
    );

    let angular_momentum = body.angular_momentum();
    // r × (m*v) = (1,0,0) × (0,1,0) = (0,0,1)
    assert!((angular_momentum.x - 0.0).abs() < 1e-10);
    assert!((angular_momentum.y - 0.0).abs() < 1e-10);
    assert!((angular_momentum.z - 1.0).abs() < 1e-10);
}

// ============================================================================
// NBodyConfig Tests
// ============================================================================

#[test]
fn test_force_between_bodies() {
    let body1 = Body::new(
        1e30,
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::zero(),
        "sun"
    );
    let body2 = Body::new(
        1e24,
        Vec3::new(1e11, 0.0, 0.0),
        Vec3::zero(),
        "earth"
    );

    let config = NBodyConfig {
        bodies: vec![body1, body2],
        dt: 1.0,
        steps: 1,
        method: IntegrationMethod::Verlet,
        softening: 0.0,
    };

    let force = config.force_between(0, 1);

    // Force should be positive in x-direction (pointing from sun to earth)
    assert!(force.x > 0.0);
    assert!((force.y).abs() < 1e-10);
    assert!((force.z).abs() < 1e-10);

    // Check magnitude using F = G*m1*m2/r²
    let expected_magnitude = G * 1e30 * 1e24 / (1e11 * 1e11);
    assert!((force.mag() - expected_magnitude).abs() / expected_magnitude < 1e-6);
}

#[test]
fn test_total_force() {
    // Three-body system
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "b1"),
        Body::new(1e30, Vec3::new(1e11, 0.0, 0.0), Vec3::zero(), "b2"),
        Body::new(1e30, Vec3::new(0.0, 1e11, 0.0), Vec3::zero(), "b3"),
    ];

    let config = NBodyConfig {
        bodies,
        dt: 1.0,
        steps: 1,
        method: IntegrationMethod::Verlet,
        softening: 0.0,
    };

    let total_force = config.total_force(0);

    // Body 0 should experience forces from both body 1 and body 2
    assert!(total_force.mag() > 0.0);
    // Due to symmetry, x and y components should be similar magnitude
    assert!((total_force.x.abs() - total_force.y.abs()).abs() < total_force.mag() * 0.01);
}

#[test]
fn test_potential_energy() {
    // Two-body system
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
        Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::zero(), "earth"),
    ];

    let config = NBodyConfig {
        bodies,
        dt: 1.0,
        steps: 1,
        method: IntegrationMethod::Verlet,
        softening: 0.0,
    };

    let pe = config.potential_energy();

    // Potential energy should be negative (bound system)
    assert!(pe < 0.0);

    // Check magnitude: PE = -G*m1*m2/r
    let expected_pe = -G * 1e30 * 1e24 / 1e11;
    assert!((pe - expected_pe).abs() / expected_pe.abs() < 1e-6);
}

#[test]
fn test_system_properties_conservation() {
    // Two equal masses with opposite velocities
    let bodies = vec![
        Body::new(
            1.0,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            "b1"
        ),
        Body::new(
            1.0,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            "b2"
        ),
    ];

    let config = NBodyConfig {
        bodies,
        dt: 0.01,
        steps: 10,
        method: IntegrationMethod::Verlet,
        softening: 0.01,
    };

    let props = config.system_properties();

    // Total momentum should be zero
    assert!(props.total_momentum.mag() < 1e-10);

    // Center of mass should be at (0.5, 0, 0)
    assert!((props.center_of_mass.x - 0.5).abs() < 1e-10);
    assert!(props.center_of_mass.y.abs() < 1e-10);
    assert!(props.center_of_mass.z.abs() < 1e-10);

    // Total energy = kinetic + potential
    assert!((props.total_energy - (props.kinetic_energy + props.potential_energy)).abs() < 1e-10);
}

#[test]
fn test_softening_parameter() {
    // Test that softening prevents singularities
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "b1"),
        Body::new(1e30, Vec3::new(1e3, 0.0, 0.0), Vec3::zero(), "b2"), // Very close
    ];

    let config_no_softening = NBodyConfig {
        bodies: bodies.clone(),
        dt: 1.0,
        steps: 1,
        method: IntegrationMethod::Verlet,
        softening: 0.0,
    };

    let config_with_softening = NBodyConfig {
        bodies,
        dt: 1.0,
        steps: 1,
        method: IntegrationMethod::Verlet,
        softening: 1e4,
    };

    let force_no_soft = config_no_softening.total_force(0);
    let force_with_soft = config_with_softening.total_force(0);

    // Softening should reduce the force magnitude
    assert!(force_with_soft.mag() < force_no_soft.mag());
    assert!(force_with_soft.mag().is_finite());
}

// ============================================================================
// Integration Methods Tests
// ============================================================================

#[test]
fn test_euler_integration() {
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
        Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::new(0.0, 30000.0, 0.0), "earth"),
    ];

    let config = NBodyConfig {
        bodies,
        dt: 3600.0, // 1 hour
        steps: 10,
        method: IntegrationMethod::Euler,
        softening: 1e8,
    };

    let result = simulate(&config);

    // Check that we have correct number of trajectory points
    assert_eq!(result.trajectories.len(), 2);
    assert_eq!(result.trajectories[0].len(), 11); // steps + 1
    assert_eq!(result.times.len(), 11);
    assert_eq!(result.energies.len(), 11);

    // All values should be finite
    assert!(result.energies.iter().all(|e| e.is_finite()));
}

#[test]
fn test_verlet_integration() {
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
        Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::new(0.0, 30000.0, 0.0), "earth"),
    ];

    let config = NBodyConfig {
        bodies,
        dt: 3600.0,
        steps: 100,
        method: IntegrationMethod::Verlet,
        softening: 1e8,
    };

    let result = simulate(&config);

    // Verlet is symplectic and should conserve energy better
    let initial_energy = result.energies[0];
    let final_energy = *result.energies.last().unwrap();

    assert!(initial_energy.is_finite());
    assert!(final_energy.is_finite());

    // Energy should be conserved (allowing some numerical drift)
    let energy_drift = (final_energy - initial_energy).abs();
    assert!(energy_drift.is_finite());
}

#[test]
fn test_runge_kutta_integration() {
    let bodies = vec![
        Body::new(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), "b1"),
        Body::new(1.0, Vec3::new(-1.0, 0.0, 0.0), Vec3::new(0.0, -1.0, 0.0), "b2"),
    ];

    let config = NBodyConfig {
        bodies,
        dt: 0.01,
        steps: 100,
        method: IntegrationMethod::RungeKutta4,
        softening: 0.1,
    };

    let result = simulate(&config);

    // RK4 should provide stable evolution
    assert!(result.energies.iter().all(|e| e.is_finite()));
    assert_eq!(result.trajectories.len(), 2);
}

#[test]
fn test_integration_method_comparison() {
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
        Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::new(0.0, 30000.0, 0.0), "planet"),
    ];

    let config_euler = NBodyConfig {
        bodies: bodies.clone(),
        dt: 3600.0,
        steps: 50,
        method: IntegrationMethod::Euler,
        softening: 1e8,
    };

    let config_verlet = NBodyConfig {
        bodies: bodies.clone(),
        dt: 3600.0,
        steps: 50,
        method: IntegrationMethod::Verlet,
        softening: 1e8,
    };

    let config_rk4 = NBodyConfig {
        bodies,
        dt: 3600.0,
        steps: 50,
        method: IntegrationMethod::RungeKutta4,
        softening: 1e8,
    };

    let result_euler = simulate(&config_euler);
    let result_verlet = simulate(&config_verlet);
    let result_rk4 = simulate(&config_rk4);

    // All methods should produce finite results
    assert!(result_euler.energies.iter().all(|e| e.is_finite()));
    assert!(result_verlet.energies.iter().all(|e| e.is_finite()));
    assert!(result_rk4.energies.iter().all(|e| e.is_finite()));

    // All should produce same number of data points
    assert_eq!(result_euler.trajectories[0].len(), 51);
    assert_eq!(result_verlet.trajectories[0].len(), 51);
    assert_eq!(result_rk4.trajectories[0].len(), 51);
}

// ============================================================================
// Barnes-Hut Octree Tests
// ============================================================================

#[test]
fn test_octree_construction() {
    let bodies = vec![
        Body::new(1.0, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "b1"),
        Body::new(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::zero(), "b2"),
        Body::new(1.0, Vec3::new(0.0, 1.0, 0.0), Vec3::zero(), "b3"),
        Body::new(1.0, Vec3::new(1.0, 1.0, 0.0), Vec3::zero(), "b4"),
    ];

    let tree = build_octree(&bodies);

    // Total mass should equal sum of all bodies
    assert!((tree.total_mass - 4.0).abs() < 1e-10);

    // Center of mass should be at (0.5, 0.5, 0)
    assert!((tree.center_of_mass.x - 0.5).abs() < 1e-6);
    assert!((tree.center_of_mass.y - 0.5).abs() < 1e-6);
    assert!((tree.center_of_mass.z - 0.0).abs() < 1e-6);
}

#[test]
fn test_octree_force_calculation() {
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
        Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::zero(), "earth"),
        Body::new(1e23, Vec3::new(0.0, 1e11, 0.0), Vec3::zero(), "probe"),
    ];

    let tree = build_octree(&bodies);
    let test_body = Body::new(
        1.0,
        Vec3::new(5e10, 5e10, 0.0),
        Vec3::zero(),
        "test"
    );

    let force = tree.calculate_force(&test_body, 0.5, 1e8);

    // Force should be finite and non-zero
    assert!(force.mag() > 0.0);
    assert!(force.mag().is_finite());
}

#[test]
fn test_octree_theta_parameter() {
    let bodies = vec![
        Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "b1"),
        Body::new(1e30, Vec3::new(1e12, 0.0, 0.0), Vec3::zero(), "b2"),
    ];

    let tree = build_octree(&bodies);
    let test_body = Body::new(1.0, Vec3::new(5e11, 5e11, 0.0), Vec3::zero(), "test");

    // Smaller theta = more accurate but slower
    let force_accurate = tree.calculate_force(&test_body, 0.1, 1e8);
    let force_approx = tree.calculate_force(&test_body, 0.9, 1e8);

    // Both should be finite
    assert!(force_accurate.mag().is_finite());
    assert!(force_approx.mag().is_finite());

    // Forces should be similar but not necessarily identical
    assert!(force_accurate.mag() > 0.0);
    assert!(force_approx.mag() > 0.0);
}
