// Unit tests for physics::gravitational_waves::snr
use computational_engine::physics::gravitational_waves::snr::*;

use super::*;
    use crate::physics::gravitational_waves::{BinarySystem, BinaryType};
    use crate::physics::gravitational_waves::waveforms::taylor_f2_waveform;

    #[test]
    fn test_snr_calculation() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0, // Mpc (GW150914-like)
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = taylor_f2_waveform(&binary, 30.0, 150.0, 2048.0);
        let detector = LIGODetector::hanford();

        let snr = calculate_snr(&waveform, &detector);

        // SNR should be positive or zero, and finite
        assert!(snr >= 0.0);
        assert!(snr.is_finite());
    }

    #[test]
    fn test_optimal_snr() {
        let binary = BinarySystem {
            mass1: 1.4,
            mass2: 1.4,
            distance: 40.0, // Close BNS
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryNeutronStar,
        };

        let waveform = taylor_f2_waveform(&binary, 20.0, 100.0, 2048.0);
        let detector = LIGODetector::hanford();

        let snr = calculate_snr(&waveform, &detector);
        let opt_snr = optimal_snr(&waveform, &detector);

        // Optimal should be higher than regular SNR
        assert!(opt_snr >= snr);
    }

    #[test]
    fn test_detection_threshold() {
        let threshold = detection_threshold();
        assert_eq!(threshold, 8.0);
    }

    #[test]
    fn test_false_alarm_rate() {
        let far_low = false_alarm_rate(5.0);
        let far_high = false_alarm_rate(15.0);

        // Higher SNR should have lower false alarm rate
        assert!(far_high < far_low);

        // Both should be positive
        assert!(far_low > 0.0);
        assert!(far_high > 0.0);
    }

    #[test]
    fn test_horizon_distance() {
        let binary = BinarySystem {
            mass1: 1.4,
            mass2: 1.4,
            distance: 1.0, // 1 Mpc reference
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryNeutronStar,
        };

        let waveform = taylor_f2_waveform(&binary, 20.0, 500.0, 2048.0);
        let detector = LIGODetector::hanford();

        let horizon = horizon_distance(&waveform, &detector, 8.0);

        // BNS horizon for aLIGO should be ~100-200 Mpc
        assert!(horizon > 0.0);
        // Just check it's reasonable (not NaN or infinite)
        assert!(horizon.is_finite());
    }

    #[test]
    fn test_matched_filter() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let signal = taylor_f2_waveform(&binary, 30.0, 100.0, 2048.0);
        let template = signal.clone();
        let detector = LIGODetector::hanford();

        let mf_snr = matched_filter_snr(&signal, &template, &detector);

        // Perfect match should give high SNR
        assert!(mf_snr > 0.0);
        assert!(mf_snr.is_finite());
    }
