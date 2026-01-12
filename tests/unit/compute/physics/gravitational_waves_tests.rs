// Unit tests for physics::gravitational_waves module
use super::*;
use std::f64::consts::PI;

// ============================================================================
// Binary System Tests
// ============================================================================

#[test]
fn test_binary_black_hole_system() {
    let system = BinarySystem {
        mass1: 36.0,
        mass2: 29.0,
        distance: 410.0, // Mpc (like GW150914)
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    assert_eq!(system.mass1, 36.0);
    assert_eq!(system.mass2, 29.0);
    assert_eq!(system.system_type, BinaryType::BinaryBlackHole);
}

#[test]
fn test_binary_neutron_star_system() {
    let system = BinarySystem {
        mass1: 1.4,
        mass2: 1.4,
        distance: 40.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryNeutronStar,
    };

    assert_eq!(system.system_type, BinaryType::BinaryNeutronStar);
    assert!(system.mass1 < 3.0); // Typical NS mass
    assert!(system.mass2 < 3.0);
}

#[test]
fn test_total_mass() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 20.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let total = system.total_mass();
    assert_eq!(total, 50.0);
}

#[test]
fn test_reduced_mass() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 20.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let mu = system.reduced_mass();
    // μ = m1*m2/(m1+m2) = 30*20/50 = 12
    assert!((mu - 12.0).abs() < 1e-10);
}

#[test]
fn test_chirp_mass_gw150914() {
    let system = BinarySystem {
        mass1: 36.0,
        mass2: 29.0,
        distance: 410.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let m_chirp = system.chirp_mass();
    // GW150914 chirp mass ~ 30 solar masses
    assert!(m_chirp > 28.0 && m_chirp < 32.0);
}

#[test]
fn test_symmetric_mass_ratio() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0, // Equal masses
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let eta = system.symmetric_mass_ratio();
    // For equal masses: η = 0.25
    assert!((eta - 0.25).abs() < 1e-10);
}

#[test]
fn test_symmetric_mass_ratio_unequal() {
    let system = BinarySystem {
        mass1: 40.0,
        mass2: 20.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let eta = system.symmetric_mass_ratio();
    // η = m1*m2/(m1+m2)^2 = 40*20/3600 = 0.222...
    assert!(eta > 0.2 && eta < 0.25);
}

#[test]
fn test_schwarzschild_radius() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let r_s = system.schwarzschild_radius();
    // 60 solar masses: r_s ~ 177 km
    assert!(r_s > 1.5e5 && r_s < 2e5);
}

#[test]
fn test_isco_frequency() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let f_isco = system.isco_frequency();
    // ISCO frequency is very high - should be positive and finite
    assert!(f_isco > 0.0 && f_isco.is_finite());
}

#[test]
fn test_time_to_coalescence() {
    let system = BinarySystem {
        mass1: 1.4,
        mass2: 1.4,
        distance: 40.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryNeutronStar,
    };

    let t_coal = system.time_to_coalescence(100.0); // At 100 Hz
    // Should be positive and finite
    assert!(t_coal > 0.0 && t_coal.is_finite());
}

// ============================================================================
// Waveform Generation Tests
// ============================================================================

#[test]
fn test_taylor_f2_waveform_generation() {
    let system = BinarySystem {
        mass1: 10.0,
        mass2: 10.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 20.0, 100.0, 1024.0);

    assert!(!waveform.is_empty());
    assert_eq!(waveform.times.len(), waveform.h_plus.len());
    assert_eq!(waveform.times.len(), waveform.h_cross.len());
    assert_eq!(waveform.times.len(), waveform.frequency.len());
}

#[test]
fn test_waveform_frequency_evolution() {
    let system = BinarySystem {
        mass1: 10.0,
        mass2: 10.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 30.0, 100.0, 1024.0);

    // Frequency should increase monotonically (chirp)
    for i in 1..waveform.frequency.len() {
        assert!(waveform.frequency[i] >= waveform.frequency[i - 1]);
    }
}

#[test]
fn test_waveform_peak_strain() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 20.0, 200.0, 1024.0);
    let peak = waveform.peak_strain();

    assert!(peak > 0.0);
    assert!(peak.is_finite());
}

#[test]
fn test_ringdown_waveform() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = ringdown_waveform(&system, 0.1, 2048.0);

    assert!(!waveform.is_empty());
    assert_eq!(waveform.len(), 204); // 0.1s * 2048Hz ≈ 204 samples (integer truncation)
}

#[test]
fn test_ringdown_exponential_decay() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = ringdown_waveform(&system, 0.1, 2048.0);

    // Amplitude should decay exponentially
    let first_amp = waveform.h_plus[0].abs();
    let last_amp = waveform.h_plus[waveform.len() - 1].abs();
    assert!(last_amp < first_amp);
}

#[test]
fn test_imr_waveform() {
    let system = BinarySystem {
        mass1: 20.0,
        mass2: 20.0,
        distance: 200.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = imr_waveform(&system, 1024.0);

    assert!(!waveform.is_empty());
    // IMR should have both inspiral and ringdown
    assert!(waveform.len() > 100);
}

// ============================================================================
// Detector Response Tests
// ============================================================================

#[test]
fn test_ligo_hanford_detector() {
    let hanford = LIGODetector::hanford();

    assert_eq!(hanford.name, "LIGO Hanford");
    assert_eq!(hanford.noise_model, NoiseModel::AdvancedLIGO);
    assert!(hanford.latitude.abs() > 0.0); // Non-zero latitude
}

#[test]
fn test_ligo_livingston_detector() {
    let livingston = LIGODetector::livingston();

    assert_eq!(livingston.name, "LIGO Livingston");
    assert_eq!(livingston.noise_model, NoiseModel::AdvancedLIGO);
    assert!(livingston.longitude.abs() > 0.0);
}

#[test]
fn test_aligo_noise_curve() {
    let detector = LIGODetector::hanford();

    // Test at minimum noise frequency (~215 Hz)
    let noise_215 = detector.noise_asd(215.0);
    assert!(noise_215 > 0.0 && noise_215 < 1e-10); // Realistic noise level ~1e-12

    // Noise should be higher at low frequencies (seismic)
    let noise_20 = detector.noise_asd(20.0);
    assert!(noise_20 > noise_215);

    // Noise should be higher at high frequencies (shot noise)
    let noise_2000 = detector.noise_asd(2000.0);
    assert!(noise_2000 > noise_215);
}

#[test]
fn test_noise_model_comparison() {
    let mut detector = LIGODetector::hanford();

    let aligo_noise = detector.noise_asd(100.0);

    detector.noise_model = NoiseModel::InitialLIGO;
    let iligo_noise = detector.noise_asd(100.0);

    detector.noise_model = NoiseModel::EinsteinTelescope;
    let et_noise = detector.noise_asd(100.0);

    // ET should be most sensitive, iLIGO least sensitive
    assert!(et_noise < aligo_noise);
    assert!(aligo_noise < iligo_noise);
}

#[test]
fn test_antenna_pattern() {
    let detector = LIGODetector::hanford();
    let (f_plus, f_cross) = detector.antenna_pattern(0.0, 0.0, 0.0);

    // Simplified model returns (1.0, 1.0)
    assert_eq!(f_plus, 1.0);
    assert_eq!(f_cross, 1.0);
}

#[test]
fn test_project_waveform() {
    let system = BinarySystem {
        mass1: 10.0,
        mass2: 10.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 30.0, 100.0, 1024.0);
    let detector = LIGODetector::hanford();

    let strain = detector.project_waveform(&waveform);

    assert_eq!(strain.len(), waveform.len());
    assert!(strain.iter().all(|&s| s.is_finite()));
}

#[test]
fn test_detector_network() {
    let network = DetectorNetwork::ligo_hl();

    assert_eq!(network.detectors.len(), 2);
    assert_eq!(network.detectors[0].name, "LIGO Hanford");
    assert_eq!(network.detectors[1].name, "LIGO Livingston");
}

#[test]
fn test_network_sensitivity() {
    let network = DetectorNetwork::ligo_hl();

    let network_asd = network.network_sensitivity(100.0);
    let h_asd = network.detectors[0].noise_asd(100.0);
    let l_asd = network.detectors[1].noise_asd(100.0);

    // Network should be more sensitive than individual detectors
    assert!(network_asd < h_asd);
    assert!(network_asd < l_asd);
}

// ============================================================================
// Signal-to-Noise Ratio Tests
// ============================================================================

#[test]
fn test_snr_calculation() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let detector = LIGODetector::hanford();

    let snr = calculate_snr(&waveform, &detector);

    assert!(snr > 0.0);
    assert!(snr.is_finite());
}

#[test]
fn test_snr_scales_with_distance() {
    let mut system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform1 = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let detector = LIGODetector::hanford();
    let snr1 = calculate_snr(&waveform1, &detector);

    // Double the distance
    system.distance = 200.0;
    let waveform2 = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let snr2 = calculate_snr(&waveform2, &detector);

    // SNR should decrease with distance (approximately as 1/d)
    assert!(snr2 < snr1);
}

#[test]
fn test_optimal_snr() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let detector = LIGODetector::hanford();

    let snr = calculate_snr(&waveform, &detector);
    let opt_snr = optimal_snr(&waveform, &detector);

    // Optimal SNR should be higher
    assert!(opt_snr >= snr);
}

#[test]
fn test_matched_filter_snr() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let signal = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let template = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let detector = LIGODetector::hanford();

    let mf_snr = matched_filter_snr(&signal, &template, &detector);

    assert!(mf_snr.is_finite());
}

#[test]
fn test_detection_threshold() {
    let threshold = detection_threshold();
    assert_eq!(threshold, 8.0);
}

#[test]
fn test_false_alarm_rate() {
    let far_low = false_alarm_rate(5.0);
    let far_high = false_alarm_rate(10.0);

    // Higher SNR should have lower false alarm rate
    assert!(far_high < far_low);
    assert!(far_high > 0.0);
}

#[test]
fn test_horizon_distance() {
    let system = BinarySystem {
        mass1: 30.0,
        mass2: 30.0,
        distance: 1.0, // 1 Mpc reference
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryBlackHole,
    };

    let waveform = taylor_f2_waveform(&system, 30.0, 200.0, 1024.0);
    let detector = LIGODetector::hanford();

    let horizon = horizon_distance(&waveform, &detector, 8.0);

    assert!(horizon > 0.0);
    assert!(horizon.is_finite());
}

// ============================================================================
// Physical Constants Tests
// ============================================================================

#[test]
fn test_physical_constants() {
    assert_eq!(C, 299792458.0); // Speed of light
    assert_eq!(G, 6.67430e-11); // Gravitational constant
    assert_eq!(M_SUN, 1.989e30); // Solar mass
    assert_eq!(PC, 3.086e16); // Parsec
}

// ============================================================================
// Edge Cases and Error Handling Tests
// ============================================================================

#[test]
fn test_empty_waveform() {
    let waveform = Waveform {
        times: vec![],
        h_plus: vec![],
        h_cross: vec![],
        frequency: vec![],
    };

    assert!(waveform.is_empty());
    assert_eq!(waveform.len(), 0);
    assert_eq!(waveform.peak_strain(), 0.0);
}

#[test]
fn test_snr_empty_waveform() {
    let waveform = Waveform {
        times: vec![],
        h_plus: vec![],
        h_cross: vec![],
        frequency: vec![],
    };

    let detector = LIGODetector::hanford();
    let snr = calculate_snr(&waveform, &detector);

    assert_eq!(snr, 0.0);
}

#[test]
fn test_extreme_mass_ratio() {
    let system = BinarySystem {
        mass1: 100.0,
        mass2: 1.0, // Extreme mass ratio
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BlackHoleNeutronStar,
    };

    let eta = system.symmetric_mass_ratio();
    // Should be much less than 0.25
    assert!(eta < 0.1);

    let m_chirp = system.chirp_mass();
    assert!(m_chirp > 0.0 && m_chirp.is_finite());
}

#[test]
fn test_high_frequency_limit() {
    let system = BinarySystem {
        mass1: 1.4,
        mass2: 1.4,
        distance: 100.0,
        inclination: 0.0,
        coalescence_phase: 0.0,
        system_type: BinaryType::BinaryNeutronStar,
    };

    // Neutron stars can reach higher frequencies
    let f_isco = system.isco_frequency();
    assert!(f_isco > 1000.0); // NS ISCO > 1 kHz
}
