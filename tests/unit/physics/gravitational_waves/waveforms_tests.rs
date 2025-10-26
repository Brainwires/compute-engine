// Unit tests for physics::gravitational_waves::waveforms
use computational_engine::physics::gravitational_waves::waveforms::*;

use super::*;
    use crate::physics::gravitational_waves::BinaryType;

    #[test]
    fn test_taylor_f2_waveform() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = taylor_f2_waveform(&binary, 20.0, 100.0, 2048.0);

        assert!(waveform.len() > 0);
        assert_eq!(waveform.times.len(), waveform.h_plus.len());
        assert_eq!(waveform.times.len(), waveform.h_cross.len());

        // Check that frequency increases monotonically
        for i in 1..waveform.frequency.len().min(100) {
            assert!(waveform.frequency[i] >= waveform.frequency[i - 1]);
        }
    }

    #[test]
    fn test_ringdown_waveform() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = ringdown_waveform(&binary, 0.05, 2048.0);

        assert!(waveform.len() > 0);

        // Ringdown should decay exponentially
        let early_amp = waveform.h_plus[10].abs();
        let late_amp = waveform.h_plus[waveform.len() - 10].abs();
        assert!(late_amp < early_amp);
    }

    #[test]
    fn test_imr_waveform() {
        let binary = BinarySystem {
            mass1: 36.0,
            mass2: 29.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = imr_waveform(&binary, 2048.0);

        assert!(waveform.len() > 100);

        // Peak strain should be finite and non-zero
        let peak = waveform.peak_strain();
        assert!(peak > 0.0 && peak.is_finite());
    }

    #[test]
    fn test_waveform_strain_amplitude() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 100.0, // Closer -> larger strain
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let waveform = taylor_f2_waveform(&binary, 30.0, 50.0, 1024.0);

        // Strain should be finite
        for &h in &waveform.h_plus {
            assert!(h.is_finite());
        }

        // At least some non-zero values
        let has_signal = waveform.h_plus.iter().any(|&h| h.abs() > 1e-30);
        assert!(has_signal);
    }
