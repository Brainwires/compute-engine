//! Comprehensive statistical physics test suite
//!
//! Tests for all statistical physics operations including:
//! - Partition Functions (3 operations)
//! - Statistical Distributions (3 operations)
//! - Thermodynamic Properties (2 operations)
//! - Phase Transitions (2 operations)

use computational_engine::engine::*;
use computational_engine::create_default_dispatcher;
use std::collections::HashMap;

// Physical constants for reference
const K_B: f64 = 1.380649e-23; // Boltzmann constant (J/K)

// ============================================================================
// PARTITION FUNCTION TESTS (3 operations)
// ============================================================================

#[test]
fn test_canonical_partition_function() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // 300 K
    params.insert("energy_levels".to_string(), serde_json::json!([0.0, 1.0e-20, 2.0e-20, 3.0e-20]));
    params.insert("degeneracies".to_string(), serde_json::json!([1, 2, 2, 1]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PartitionFunctionCanonical)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Canonical partition function should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let z = output.result.get("partition_function").and_then(|v| v.as_f64()).unwrap();
        assert!(z > 0.0, "Partition function should be positive");

        let u = output.result.get("internal_energy").and_then(|v| v.as_f64()).unwrap();
        assert!(u >= 0.0, "Internal energy should be non-negative");

        let s = output.result.get("entropy").and_then(|v| v.as_f64()).unwrap();
        assert!(s > 0.0, "Entropy should be positive");

        let heat_capacity = output.result.get("heat_capacity").and_then(|v| v.as_f64()).unwrap();
        assert!(heat_capacity >= 0.0, "Heat capacity should be non-negative");
    }
}

#[test]
fn test_partition_function_general() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("ensemble".to_string(), serde_json::json!("canonical"));
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    params.insert("energy_levels".to_string(), serde_json::json!([0.0, 1.0e-20, 2.0e-20]));
    params.insert("degeneracies".to_string(), serde_json::json!([1, 1, 1]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PartitionFunction)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "General partition function should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let ensemble_type = output.result.get("ensemble_type").and_then(|v| v.as_str()).unwrap();
        assert_eq!(ensemble_type, "canonical");

        let z = output.result.get("partition_function").and_then(|v| v.as_f64()).unwrap();
        assert!(z > 0.0, "Partition function should be positive");
    }
}

#[test]
fn test_grand_canonical_partition() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // K
    params.insert("volume".to_string(), serde_json::json!(1.0e-3)); // m³ (1 liter)
    params.insert("chemical_potential".to_string(), serde_json::json!(-1.0e-20)); // J
    params.insert("particle_type".to_string(), serde_json::json!("classical"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PartitionFunctionGrandCanonical)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Grand canonical partition should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let grand_potential = output.result.get("grand_potential").and_then(|v| v.as_f64()).unwrap();
        assert!(grand_potential.is_finite(), "Grand potential should be finite");

        let n_avg = output.result.get("average_particle_number").and_then(|v| v.as_f64()).unwrap();
        assert!(n_avg > 0.0, "Average particle number should be positive");

        let pressure = output.result.get("pressure").and_then(|v| v.as_f64()).unwrap();
        assert!(pressure > 0.0, "Pressure should be positive");

        let entropy = output.result.get("entropy").and_then(|v| v.as_f64()).unwrap();
        assert!(entropy >= 0.0, "Entropy should be non-negative");
    }
}

// ============================================================================
// STATISTICAL DISTRIBUTIONS TESTS (3 operations)
// ============================================================================

#[test]
fn test_maxwell_boltzmann_distribution() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // 300 K
    params.insert("mass".to_string(), serde_json::json!(4.65e-26)); // N₂ molecule mass (kg)
    params.insert("velocity".to_string(), serde_json::json!(500.0)); // m/s

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::MaxwellBoltzmann)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Maxwell-Boltzmann should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let prob_density = output.result.get("probability_density").and_then(|v| v.as_f64());
        assert!(prob_density.is_some(), "Should have probability density");
        assert!(prob_density.unwrap() > 0.0, "Probability density should be positive");

        let v_avg = output.result.get("average_speed").and_then(|v| v.as_f64()).unwrap();
        let v_rms = output.result.get("rms_speed").and_then(|v| v.as_f64()).unwrap();
        let v_p = output.result.get("most_probable_speed").and_then(|v| v.as_f64()).unwrap();

        // Physical ordering: v_p < v_avg < v_rms
        assert!(v_p < v_avg, "Most probable speed < average speed");
        assert!(v_avg < v_rms, "Average speed < RMS speed");
        assert!(v_rms > 400.0 && v_rms < 600.0, "RMS speed should be ~500 m/s for N₂ at 300K");
    }
}

#[test]
fn test_maxwell_boltzmann_without_velocity() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    params.insert("mass".to_string(), serde_json::json!(4.65e-26));
    // No velocity parameter - should return characteristic speeds only

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::MaxwellBoltzmann)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Maxwell-Boltzmann without velocity should succeed");

    if let Ok(ToolResponse::Compute(output)) = result {
        let prob_density = output.result.get("probability_density");
        assert!(prob_density.is_none(), "Should not have probability density without velocity");

        let v_avg = output.result.get("average_speed").and_then(|v| v.as_f64()).unwrap();
        assert!(v_avg > 0.0, "Should still return characteristic speeds");
    }
}

#[test]
fn test_fermi_dirac_distribution() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("energy".to_string(), serde_json::json!(1.6e-19)); // ~1 eV
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // K
    params.insert("chemical_potential".to_string(), serde_json::json!(1.5e-19)); // Fermi energy ~0.94 eV

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::FermiDirac)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Fermi-Dirac should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let occupation = output.result.get("occupation_probability").and_then(|v| v.as_f64()).unwrap();
        assert!(occupation >= 0.0 && occupation <= 1.0, "Occupation probability must be in [0,1]");

        let fermi_temp = output.result.get("fermi_temperature").and_then(|v| v.as_f64()).unwrap();
        assert!(fermi_temp > 0.0, "Fermi temperature should be positive");
        // Fermi temperature ~11,000 K for E_F ~ 1.5 eV
        assert!(fermi_temp > 10000.0, "Fermi temperature should be ~10,000+ K");

        let fermi_energy = output.result.get("fermi_energy").and_then(|v| v.as_f64()).unwrap();
        assert!((fermi_energy - 1.5e-19).abs() < 1e-21, "Fermi energy should match chemical potential");
    }
}

#[test]
fn test_bose_einstein_distribution() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("energy".to_string(), serde_json::json!(1.0e-20)); // J
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // K
    params.insert("chemical_potential".to_string(), serde_json::json!(0.5e-20)); // μ < E

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::BoseEinstein)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Bose-Einstein should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let occupation = output.result.get("occupation_number").and_then(|v| v.as_f64()).unwrap();
        assert!(occupation > 0.0, "Occupation number should be positive for bosons");
        assert!(occupation.is_finite(), "Occupation should be finite");
    }
}

#[test]
fn test_bose_einstein_requires_mu_less_than_energy() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("energy".to_string(), serde_json::json!(1.0e-20));
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    params.insert("chemical_potential".to_string(), serde_json::json!(1.5e-20)); // μ > E (invalid!)

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::BoseEinstein)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should reject μ ≥ E for bosons");
}

// ============================================================================
// THERMODYNAMIC PROPERTIES TESTS (2 operations)
// ============================================================================

#[test]
fn test_chemical_potential_classical() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0)); // K
    params.insert("pressure".to_string(), serde_json::json!(101325.0)); // Pa (1 atm)
    params.insert("particle_density".to_string(), serde_json::json!(2.5e25)); // particles/m³
    params.insert("particle_type".to_string(), serde_json::json!("classical"));
    params.insert("mass".to_string(), serde_json::json!(4.65e-26)); // N₂

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::ChemicalPotential)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Chemical potential should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let mu = output.result.get("chemical_potential").and_then(|v| v.as_f64()).unwrap();
        assert!(mu.is_finite(), "Chemical potential should be finite");

        let lambda = output.result.get("thermal_wavelength").and_then(|v| v.as_f64()).unwrap();
        assert!(lambda > 0.0, "Thermal wavelength should be positive");
        // Thermal wavelength ~1-10 nm for typical gases at room temperature
        assert!(lambda > 1.0e-11 && lambda < 1.0e-8, "Thermal wavelength should be ~nm scale");

        let fugacity = output.result.get("fugacity").and_then(|v| v.as_f64()).unwrap();
        assert!(fugacity > 0.0, "Fugacity should be positive");
    }
}

#[test]
fn test_chemical_potential_fermion() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    params.insert("pressure".to_string(), serde_json::json!(101325.0));
    params.insert("particle_density".to_string(), serde_json::json!(1.0e28)); // High density (degenerate)
    params.insert("particle_type".to_string(), serde_json::json!("fermion"));
    params.insert("mass".to_string(), serde_json::json!(9.109e-31)); // Electron mass

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::ChemicalPotential)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Fermion chemical potential should succeed");

    if let Ok(ToolResponse::Compute(output)) = result {
        let mu = output.result.get("chemical_potential").and_then(|v| v.as_f64()).unwrap();
        // For degenerate fermions, μ ≈ E_F > 0
        assert!(mu > 0.0, "Fermi energy should be positive for degenerate gas");
    }
}

#[test]
fn test_fugacity_coefficient() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    params.insert("pressure".to_string(), serde_json::json!(101325.0)); // 1 atm
    params.insert("particle_density".to_string(), serde_json::json!(2.5e25));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::FugacityCoefficient)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Fugacity coefficient should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let fugacity = output.result.get("fugacity").and_then(|v| v.as_f64()).unwrap();
        assert!(fugacity > 0.0, "Fugacity should be positive");

        let fugacity_coeff = output.result.get("fugacity_coefficient").and_then(|v| v.as_f64()).unwrap();
        assert!(fugacity_coeff > 0.0, "Fugacity coefficient should be positive");
        // For ideal gas, φ ≈ 1
        assert!(fugacity_coeff.is_finite(), "Fugacity coefficient should be finite");

        let mu = output.result.get("chemical_potential").and_then(|v| v.as_f64()).unwrap();
        assert!(mu.is_finite(), "Chemical potential should be finite");
    }
}

// ============================================================================
// PHASE TRANSITIONS TESTS (2 operations)
// ============================================================================

#[test]
fn test_phase_transition_ising_below_tc() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("model".to_string(), serde_json::json!("ising"));
    params.insert("temperature".to_string(), serde_json::json!(400.0)); // K
    params.insert("critical_temperature".to_string(), serde_json::json!(500.0)); // T_c

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PhaseTransition)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Phase transition should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let order_param = output.result.get("order_parameter").and_then(|v| v.as_f64()).unwrap();
        assert!(order_param > 0.0, "Order parameter should be positive below T_c");

        let phase = output.result.get("phase").and_then(|v| v.as_str()).unwrap();
        assert_eq!(phase, "ordered", "Should be in ordered phase below T_c");

        let correlation_length = output.result.get("correlation_length").and_then(|v| v.as_f64()).unwrap();
        assert!(correlation_length > 0.0, "Correlation length should be finite below T_c");
    }
}

#[test]
fn test_phase_transition_ising_above_tc() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("model".to_string(), serde_json::json!("ising"));
    params.insert("temperature".to_string(), serde_json::json!(600.0)); // T > T_c
    params.insert("critical_temperature".to_string(), serde_json::json!(500.0));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PhaseTransition)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Phase transition above T_c should succeed");

    if let Ok(ToolResponse::Compute(output)) = result {
        let order_param = output.result.get("order_parameter").and_then(|v| v.as_f64()).unwrap();
        assert_eq!(order_param, 0.0, "Order parameter should be zero above T_c");

        let phase = output.result.get("phase").and_then(|v| v.as_str()).unwrap();
        assert_eq!(phase, "disordered", "Should be in disordered phase above T_c");
    }
}

#[test]
fn test_critical_phenomena_ising_2d() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(495.0)); // Just below T_c
    params.insert("critical_temperature".to_string(), serde_json::json!(500.0));
    params.insert("model".to_string(), serde_json::json!("ising_2d"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::CriticalPhenomena)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Critical phenomena should succeed: {:?}", result);

    if let Ok(ToolResponse::Compute(output)) = result {
        let exponents = output.result.get("critical_exponents").and_then(|v| v.as_object()).unwrap();

        // Check 2D Ising exact exponents
        let beta = exponents.get("beta").and_then(|v| v.as_f64()).unwrap();
        assert!((beta - 0.125).abs() < 0.01, "2D Ising β = 1/8");

        let nu = exponents.get("nu").and_then(|v| v.as_f64()).unwrap();
        assert!((nu - 1.0).abs() < 0.01, "2D Ising ν = 1");

        let gamma = exponents.get("gamma").and_then(|v| v.as_f64()).unwrap();
        assert!((gamma - 1.75).abs() < 0.01, "2D Ising γ = 7/4");

        let reduced_t = output.result.get("reduced_temperature").and_then(|v| v.as_f64()).unwrap();
        assert!(reduced_t < 0.0, "Should be below T_c");
        assert!((reduced_t + 0.01).abs() < 0.001, "t = (T-T_c)/T_c = -0.01");
    }
}

#[test]
fn test_critical_phenomena_mean_field() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(510.0)); // Above T_c
    params.insert("critical_temperature".to_string(), serde_json::json!(500.0));
    params.insert("model".to_string(), serde_json::json!("mean_field"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::CriticalPhenomena)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Mean field critical phenomena should succeed");

    if let Ok(ToolResponse::Compute(output)) = result {
        let exponents = output.result.get("critical_exponents").and_then(|v| v.as_object()).unwrap();

        // Check mean field classical exponents
        let beta = exponents.get("beta").and_then(|v| v.as_f64()).unwrap();
        assert_eq!(beta, 0.5, "Mean field β = 1/2");

        let gamma = exponents.get("gamma").and_then(|v| v.as_f64()).unwrap();
        assert_eq!(gamma, 1.0, "Mean field γ = 1");

        let delta = exponents.get("delta").and_then(|v| v.as_f64()).unwrap();
        assert_eq!(delta, 3.0, "Mean field δ = 3");

        let order_param = output.result.get("order_parameter").and_then(|v| v.as_f64()).unwrap();
        assert_eq!(order_param, 0.0, "Order parameter zero above T_c");
    }
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_partition_function_requires_ensemble() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    // Missing ensemble parameter

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PartitionFunction)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should require ensemble parameter");
}

#[test]
fn test_maxwell_boltzmann_requires_mass() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    // Missing mass parameter

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::MaxwellBoltzmann)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should require mass parameter");
}

#[test]
fn test_canonical_partition_requires_energy_levels() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("temperature".to_string(), serde_json::json!(300.0));
    // Missing energy_levels and degeneracies

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Physics(PhysicsOp::StatisticalPhysics(StatPhysicsOp::PartitionFunctionCanonical)),
        data: serde_json::json!({}),
        parameters: params,
    });

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "Should require energy levels");
}
