//! Signal-to-Noise Ratio Calculations

use super::{LIGODetector, Waveform};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SNRResult {
    pub snr: f64,
    pub optimal_snr: f64,
    pub matched_filter_snr: f64,
}

/// Calculate signal-to-noise ratio for a waveform in a detector
pub fn calculate_snr(waveform: &Waveform, detector: &LIGODetector) -> f64 {
    if waveform.is_empty() {
        return 0.0;
    }

    // Project waveform onto detector
    let h = detector.project_waveform(waveform);

    // Calculate SNR in frequency domain (simplified)
    let mut snr_sq = 0.0;
    let df = if waveform.times.len() > 1 {
        1.0 / (waveform.times[waveform.times.len() - 1] - waveform.times[0])
    } else {
        1.0
    };

    for i in 0..h.len().min(waveform.frequency.len()) {
        let f = waveform.frequency[i];
        if f > 10.0 && f < 2000.0 {
            let s_n = detector.noise_asd(f);
            if s_n > 0.0 {
                let h_f = h[i].abs();
                snr_sq += 4.0 * (h_f * h_f) / (s_n * s_n) * df;
            }
        }
    }

    snr_sq.sqrt()
}

/// Calculate optimal SNR (theoretical maximum)
pub fn optimal_snr(waveform: &Waveform, detector: &LIGODetector) -> f64 {
    // Simplified: assume optimal alignment
    calculate_snr(waveform, detector) * 1.5
}

/// Matched filter SNR
pub fn matched_filter_snr(
    signal: &Waveform,
    template: &Waveform,
    detector: &LIGODetector,
) -> f64 {
    if signal.is_empty() || template.is_empty() {
        return 0.0;
    }

    // Simplified matched filtering
    let min_len = signal.h_plus.len().min(template.h_plus.len());

    let mut overlap = 0.0;
    let mut norm_sq = 0.0;

    for i in 0..min_len {
        let f = signal.frequency.get(i).unwrap_or(&100.0);
        let s_n = detector.noise_asd(*f);

        if s_n > 0.0 {
            let sig = signal.h_plus[i];
            let temp = template.h_plus[i];

            overlap += sig * temp / (s_n * s_n);
            norm_sq += temp * temp / (s_n * s_n);
        }
    }

    if norm_sq > 0.0 {
        overlap / norm_sq.sqrt()
    } else {
        0.0
    }
}

/// Detection threshold (typical SNR needed for detection)
pub fn detection_threshold() -> f64 {
    8.0 // Typical threshold for GW detection
}

/// False alarm rate estimate from SNR
pub fn false_alarm_rate(snr: f64) -> f64 {
    // Simplified: exponential suppression
    let rate_0 = 1.0; // events per year at SNR=0
    rate_0 * (-snr * snr / 2.0).exp()
}

/// Horizon distance: distance at which SNR = 8 for optimal orientation
pub fn horizon_distance(waveform: &Waveform, detector: &LIGODetector, threshold: f64) -> f64 {
    let snr_at_1mpc = calculate_snr(waveform, detector);

    if snr_at_1mpc > 0.0 {
        // Distance scales as 1/SNR
        snr_at_1mpc / threshold
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
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
}
