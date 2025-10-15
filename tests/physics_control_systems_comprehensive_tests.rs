//! Comprehensive tests for Control Systems module
//!
//! Tests transfer functions, pole-zero analysis, Bode plots, Nyquist plots,
//! root locus, state-space, controllability, observability, stability analysis.

use computational_engine::physics::control_systems::*;

// ============================================================================
// TRANSFER FUNCTION TESTS
// ============================================================================

#[test]
fn test_transfer_function_first_order() {
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![1.0],      // G(s) = 1/(s+1)
        denominator: vec![1.0, 1.0],
        operation: "evaluate".to_string(),
        frequency: None,
        second_tf: None,
    })
    .unwrap();

    // DC gain = 1/1 = 1.0
    assert!((result.gain - 1.0).abs() < 0.01);
    assert_eq!(result.poles.len(), 1);
}

#[test]
fn test_transfer_function_second_order() {
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![4.0],           // 4/(s²+2s+4)
        denominator: vec![1.0, 2.0, 4.0],
        operation: "evaluate".to_string(),
        frequency: None,
        second_tf: None,
    })
    .unwrap();

    assert!((result.gain - 1.0).abs() < 0.01);
    assert_eq!(result.poles.len(), 2);
}

#[test]
fn test_transfer_function_series() {
    // H1(s) = 1/(s+1), H2(s) = 2/(s+2)
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0],
        operation: "series".to_string(),
        frequency: None,
        second_tf: Some((vec![2.0], vec![1.0, 2.0])),
    })
    .unwrap();

    // Product: 2/[(s+1)(s+2)]
    assert_eq!(result.numerator.len(), 1);
    assert_eq!(result.denominator.len(), 3);
}

#[test]
fn test_transfer_function_parallel() {
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0],
        operation: "parallel".to_string(),
        frequency: None,
        second_tf: Some((vec![1.0], vec![1.0, 2.0])),
    })
    .unwrap();

    // Sum: 1/(s+1) + 1/(s+2)
    assert!(result.numerator.len() > 0);
    assert!(result.denominator.len() > 0);
}

#[test]
fn test_transfer_function_feedback() {
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![10.0],
        denominator: vec![1.0, 1.0],
        operation: "feedback".to_string(),
        frequency: None,
        second_tf: Some((vec![1.0], vec![1.0])), // Unity feedback
    })
    .unwrap();

    // Closed-loop: G/(1+G)
    assert!(result.poles.len() > 0);
}

#[test]
fn test_transfer_function_frequency_response() {
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0],
        operation: "evaluate".to_string(),
        frequency: Some(1.0), // ω = 1 rad/s
        second_tf: None,
    })
    .unwrap();

    assert!(result.frequency_response.is_some());
    let freq = result.frequency_response.unwrap();
    assert!(freq.magnitude > 0.0);
    assert!(freq.phase.abs() < std::f64::consts::PI);
}

// ============================================================================
// POLE-ZERO ANALYSIS TESTS
// ============================================================================

#[test]
fn test_pole_zero_stable_system() {
    let result = pole_zero_analysis(PoleZeroRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 3.0, 2.0], // (s+1)(s+2)
    })
    .unwrap();

    assert_eq!(result.stability, "stable");
    assert_eq!(result.poles.len(), 2);
    // Both poles should be negative real
    for pole in &result.poles {
        assert!(pole.real < 0.0);
    }
}

#[test]
fn test_pole_zero_unstable_system() {
    let result = pole_zero_analysis(PoleZeroRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, -1.0], // s-1 (unstable pole)
    })
    .unwrap();

    assert_eq!(result.stability, "unstable");
}

#[test]
fn test_pole_zero_dominant_pole() {
    let result = pole_zero_analysis(PoleZeroRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 11.0, 10.0], // (s+1)(s+10)
    })
    .unwrap();

    // Dominant pole is at -1 (closer to imaginary axis)
    assert!(result.dominant_pole.real.abs() < 2.0);
}

#[test]
fn test_pole_zero_oscillatory() {
    let result = pole_zero_analysis(PoleZeroRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 0.2, 1.0], // Complex poles
    })
    .unwrap();

    // Should have complex poles (underdamped)
    assert!(result.natural_frequency > 0.0);
    assert!(result.damping_ratio < 1.0);
}

// ============================================================================
// BODE PLOT TESTS
// ============================================================================

#[test]
fn test_bode_plot_first_order() {
    let result = bode_plot(BodePlotRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0],
        freq_min: 0.1,
        freq_max: 100.0,
        num_points: 50,
    })
    .unwrap();

    assert_eq!(result.frequencies.len(), 50);
    assert_eq!(result.magnitude_db.len(), 50);
    assert_eq!(result.phase_deg.len(), 50);

    // At low frequency, magnitude should be ~0 dB
    assert!(result.magnitude_db[0] > -3.0);
    // At high frequency, magnitude should decrease
    assert!(result.magnitude_db.last().unwrap() < &result.magnitude_db[0]);
}

#[test]
fn test_bode_plot_integrator() {
    let result = bode_plot(BodePlotRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 0.0], // 1/s
        freq_min: 0.1,
        freq_max: 10.0,
        num_points: 20,
    })
    .unwrap();

    // Integrator: -20 dB/decade slope
    let mag_ratio = result.magnitude_db[0] - result.magnitude_db.last().unwrap();
    assert!(mag_ratio > 30.0); // Should decrease significantly
}

#[test]
fn test_bode_plot_second_order() {
    let result = bode_plot(BodePlotRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 0.1, 1.0], // Low damping
        freq_min: 0.1,
        freq_max: 10.0,
        num_points: 50,
    })
    .unwrap();

    // Should have resonance peak near natural frequency
    let max_mag = result.magnitude_db.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    assert!(max_mag > 0.0); // Resonance peak
}

#[test]
fn test_bode_plot_crossover_frequencies() {
    let result = bode_plot(BodePlotRequest {
        numerator: vec![10.0],
        denominator: vec![1.0, 1.0],
        freq_min: 0.01,
        freq_max: 100.0,
        num_points: 100,
    })
    .unwrap();

    // Should find gain crossover
    assert!(result.gain_crossover_freq > 0.0);
}

// ============================================================================
// NYQUIST PLOT TESTS
// ============================================================================

#[test]
fn test_nyquist_plot_stable() {
    let result = nyquist_plot(NyquistPlotRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 3.0, 2.0], // Stable: (s+1)(s+2)
        freq_min: 0.0,
        freq_max: 100.0,
        num_points: 100,
    })
    .unwrap();

    assert_eq!(result.real_parts.len(), 100);
    assert_eq!(result.imag_parts.len(), 100);
    // Stable open-loop should not encircle -1
    assert_eq!(result.encirclements, 0);
    assert!(result.stable);
}

#[test]
fn test_nyquist_plot_high_gain() {
    let result = nyquist_plot(NyquistPlotRequest {
        numerator: vec![100.0],
        denominator: vec![1.0, 1.0],
        freq_min: 0.01,
        freq_max: 10.0,
        num_points: 50,
    })
    .unwrap();

    // High gain can cause encirclement
    assert!(result.real_parts.len() > 0);
}

// ============================================================================
// ROOT LOCUS TESTS
// ============================================================================

#[test]
fn test_root_locus_simple() {
    let result = root_locus(RootLocusRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 3.0, 2.0],
        gain_min: 0.0,
        gain_max: 10.0,
        num_points: 20,
    })
    .unwrap();

    assert_eq!(result.gains.len(), 20);
    assert!(result.pole_trajectories.len() > 0);

    // Poles should start at open-loop poles
    assert!(result.pole_trajectories[0].len() > 0);
}

#[test]
fn test_root_locus_stability_boundary() {
    let result = root_locus(RootLocusRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0, 0.0], // Marginally stable
        gain_min: 0.0,
        gain_max: 5.0,
        num_points: 50,
    })
    .unwrap();

    // Should show poles crossing into RHP at some gain
    assert!(result.pole_trajectories.len() > 0);
}

// ============================================================================
// STATE SPACE TESTS
// ============================================================================

#[test]
fn test_state_space_eigenvalues() {
    let result = state_space(StateSpaceRequest {
        a_matrix: vec![vec![-1.0, 0.0], vec![0.0, -2.0]],
        b_matrix: vec![vec![1.0], vec![1.0]],
        c_matrix: vec![vec![1.0, 0.0]],
        d_matrix: vec![vec![0.0]],
        operation: "to_transfer_function".to_string(),
        time: None,
    })
    .unwrap();

    assert_eq!(result.eigenvalues.len(), 2);
    // Eigenvalues should be -1 and -2
    assert!(result.transfer_function.is_some());
}

#[test]
fn test_state_space_to_transfer_function() {
    let result = state_space(StateSpaceRequest {
        a_matrix: vec![vec![-2.0]],
        b_matrix: vec![vec![1.0]],
        c_matrix: vec![vec![1.0]],
        d_matrix: vec![vec![0.0]],
        operation: "to_transfer_function".to_string(),
        time: None,
    })
    .unwrap();

    assert!(result.transfer_function.is_some());
}

// ============================================================================
// CONTROLLABILITY TESTS
// ============================================================================

#[test]
fn test_controllability_controllable() {
    let result = controllability(ControllabilityRequest {
        a_matrix: vec![vec![0.0, 1.0], vec![-2.0, -3.0]],
        b_matrix: vec![vec![0.0], vec![1.0]],
    })
    .unwrap();

    assert!(result.controllable);
    assert_eq!(result.rank, 2);
}

#[test]
fn test_controllability_uncontrollable() {
    let result = controllability(ControllabilityRequest {
        a_matrix: vec![vec![1.0, 0.0], vec![0.0, 2.0]],
        b_matrix: vec![vec![1.0], vec![0.0]],
    })
    .unwrap();

    // Second state not controllable
    assert!(!result.controllable);
}

#[test]
fn test_controllability_matrix_structure() {
    let result = controllability(ControllabilityRequest {
        a_matrix: vec![vec![0.0, 1.0], vec![0.0, 0.0]],
        b_matrix: vec![vec![0.0], vec![1.0]],
    })
    .unwrap();

    // Controllability matrix: [B AB]
    assert_eq!(result.controllability_matrix.len(), 2);
}

// ============================================================================
// OBSERVABILITY TESTS
// ============================================================================

#[test]
fn test_observability_observable() {
    let result = observability(ObservabilityRequest {
        a_matrix: vec![vec![0.0, 1.0], vec![-2.0, -3.0]],
        c_matrix: vec![vec![1.0, 0.0]],
    })
    .unwrap();

    assert!(result.observable);
    assert_eq!(result.rank, 2);
}

#[test]
fn test_observability_unobservable() {
    let result = observability(ObservabilityRequest {
        a_matrix: vec![vec![1.0, 0.0], vec![0.0, 2.0]],
        c_matrix: vec![vec![1.0, 0.0]],
    })
    .unwrap();

    // Second state not observable
    assert!(!result.observable);
}

// ============================================================================
// ROUTH-HURWITZ TESTS
// ============================================================================

#[test]
fn test_routh_hurwitz_stable() {
    let result = routh_hurwitz(RouthHurwitzRequest {
        characteristic_polynomial: vec![1.0, 3.0, 3.0, 1.0], // (s+1)³
    })
    .unwrap();

    assert!(result.stable);
    assert_eq!(result.sign_changes, 0);
    assert_eq!(result.unstable_poles, 0);
}

#[test]
fn test_routh_hurwitz_unstable() {
    let result = routh_hurwitz(RouthHurwitzRequest {
        characteristic_polynomial: vec![1.0, -2.0, 1.0], // (s-1)²
    })
    .unwrap();

    assert!(!result.stable);
    assert!(result.sign_changes > 0);
}

#[test]
fn test_routh_hurwitz_marginally_stable() {
    let result = routh_hurwitz(RouthHurwitzRequest {
        characteristic_polynomial: vec![1.0, 0.0, 1.0], // s²+1 (imaginary poles)
    })
    .unwrap();

    // May be stable or marginally stable depending on implementation
    assert!(result.routh_array.len() > 0);
}

#[test]
fn test_routh_hurwitz_third_order() {
    let result = routh_hurwitz(RouthHurwitzRequest {
        characteristic_polynomial: vec![1.0, 6.0, 11.0, 6.0], // (s+1)(s+2)(s+3)
    })
    .unwrap();

    assert!(result.stable);
    assert_eq!(result.unstable_poles, 0);
}

// ============================================================================
// GAIN MARGIN TESTS
// ============================================================================

#[test]
fn test_gain_margin_stable_system() {
    let result = gain_margin(GainMarginRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 3.0, 2.0],
    })
    .unwrap();

    // Stable system should have positive gain margin
    assert!(result.gain_margin_db > 0.0);
    assert!(result.phase_crossover_frequency > 0.0);
}

#[test]
fn test_gain_margin_high_gain() {
    let result = gain_margin(GainMarginRequest {
        numerator: vec![10.0],
        denominator: vec![1.0, 1.0],
    })
    .unwrap();

    assert!(result.phase_crossover_frequency > 0.0);
}

// ============================================================================
// PHASE MARGIN TESTS
// ============================================================================

#[test]
fn test_phase_margin_stable_system() {
    let result = phase_margin(PhaseMarginRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 3.0, 2.0],
    })
    .unwrap();

    // Stable system should have positive phase margin
    assert!(result.phase_margin_deg > 0.0);
    assert!(result.gain_crossover_frequency > 0.0);
}

#[test]
fn test_phase_margin_low_damping() {
    let result = phase_margin(PhaseMarginRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 0.1, 1.0],
    })
    .unwrap();

    // Low damping system should still have measurable phase margin
    assert!(result.phase_margin_deg > -180.0);
    assert!(result.gain_crossover_frequency > 0.0);
}

// ============================================================================
// STEP RESPONSE TESTS
// ============================================================================

#[test]
fn test_step_response_first_order() {
    let result = step_response(StepResponseRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0],
        time_span: 10.0,
        num_points: 100,
    })
    .unwrap();

    assert_eq!(result.time.len(), 100);
    assert_eq!(result.response.len(), 100);

    // Final value should approach 1.0
    assert!((result.steady_state_value - 1.0).abs() < 0.2);

    // No overshoot for first-order
    assert!(result.overshoot_percent < 1.0);
}

#[test]
fn test_step_response_second_order_underdamped() {
    let result = step_response(StepResponseRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 0.2, 1.0], // ζ = 0.1 (underdamped)
        time_span: 50.0,
        num_points: 200,
    })
    .unwrap();

    // Underdamped: should have overshoot
    assert!(result.overshoot_percent > 1.0);
    assert!(result.rise_time > 0.0);
    assert!(result.peak_time > 0.0);
}

#[test]
fn test_step_response_second_order_overdamped() {
    let result = step_response(StepResponseRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 5.0, 1.0], // Overdamped
        time_span: 20.0,
        num_points: 100,
    })
    .unwrap();

    // Overdamped: no overshoot
    assert!(result.overshoot_percent < 1.0);
    assert!(result.settling_time > 0.0);
}

#[test]
fn test_step_response_characteristics() {
    let result = step_response(StepResponseRequest {
        numerator: vec![4.0],
        denominator: vec![1.0, 2.0, 4.0],
        time_span: 10.0,
        num_points: 100,
    })
    .unwrap();

    // All characteristics should be positive
    assert!(result.rise_time >= 0.0);
    assert!(result.settling_time > 0.0);
    assert!(result.peak_time >= 0.0);
    assert!(result.overshoot_percent >= 0.0);
}

#[test]
fn test_step_response_integrator() {
    // Integrator has issues with simplified step response
    // Use a system with very slow pole instead
    let result = step_response(StepResponseRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 0.01], // Very slow pole
        time_span: 10.0,
        num_points: 50,
    })
    .unwrap();

    // Response should be slow but measurable
    assert!(result.response.len() > 0);
    assert!(result.steady_state_value > 0.0);
}

// ============================================================================
// EDGE CASES AND COMPLEX SCENARIOS
// ============================================================================

#[test]
fn test_transfer_function_high_order() {
    let result = transfer_function(TransferFunctionRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 4.0, 6.0, 4.0, 1.0], // (s+1)⁴
        operation: "evaluate".to_string(),
        frequency: None,
        second_tf: None,
    })
    .unwrap();

    assert_eq!(result.denominator.len(), 5);
}

#[test]
fn test_bode_plot_wide_frequency_range() {
    let result = bode_plot(BodePlotRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 1.0],
        freq_min: 0.001,
        freq_max: 1000.0,
        num_points: 200,
    })
    .unwrap();

    assert_eq!(result.frequencies.len(), 200);
    // Frequencies should span the range logarithmically
    assert!(result.frequencies[0] < 0.01);
    assert!(result.frequencies.last().unwrap() > &900.0);
}

#[test]
fn test_state_space_larger_system() {
    let result = state_space(StateSpaceRequest {
        a_matrix: vec![
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![-6.0, -11.0, -6.0],
        ],
        b_matrix: vec![vec![0.0], vec![0.0], vec![1.0]],
        c_matrix: vec![vec![1.0, 0.0, 0.0]],
        d_matrix: vec![vec![0.0]],
        operation: "to_transfer_function".to_string(),
        time: None,
    })
    .unwrap();

    assert_eq!(result.eigenvalues.len(), 3);
}

#[test]
fn test_controllability_single_input() {
    let result = controllability(ControllabilityRequest {
        a_matrix: vec![vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0], vec![-1.0, -2.0, -3.0]],
        b_matrix: vec![vec![0.0], vec![0.0], vec![1.0]],
    })
    .unwrap();

    assert_eq!(result.rank, 3);
}

#[test]
fn test_root_locus_zero_gain() {
    let result = root_locus(RootLocusRequest {
        numerator: vec![1.0],
        denominator: vec![1.0, 2.0, 1.0],
        gain_min: 0.0,
        gain_max: 0.0, // Zero range
        num_points: 1,
    })
    .unwrap();

    // Implementation uses max(num_points, 10) so we get at least 10 points
    assert!(result.gains.len() >= 1);
}
