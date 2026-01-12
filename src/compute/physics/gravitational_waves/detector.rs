//! LIGO Detector Response and Noise Curves

use super::{Waveform, PI};
use serde::{Deserialize, Serialize};

/// LIGO detector configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LIGODetector {
    pub name: String,
    pub latitude: f64,      // radians
    pub longitude: f64,     // radians
    pub arm_azimuth: f64,   // radians
    pub noise_model: NoiseModel,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum NoiseModel {
    AdvancedLIGO,    // aLIGO design sensitivity
    InitialLIGO,      // Initial LIGO
    EinsteinTelescope, // ET design
}

impl LIGODetector {
    /// LIGO Hanford
    pub fn hanford() -> Self {
        Self {
            name: "LIGO Hanford".to_string(),
            latitude: 46.4551 * PI / 180.0,
            longitude: -119.4078 * PI / 180.0,
            arm_azimuth: 0.0,
            noise_model: NoiseModel::AdvancedLIGO,
        }
    }

    /// LIGO Livingston
    pub fn livingston() -> Self {
        Self {
            name: "LIGO Livingston".to_string(),
            latitude: 30.5629 * PI / 180.0,
            longitude: -90.7742 * PI / 180.0,
            arm_azimuth: 0.0,
            noise_model: NoiseModel::AdvancedLIGO,
        }
    }

    /// Noise amplitude spectral density (strain/√Hz)
    pub fn noise_asd(&self, frequency: f64) -> f64 {
        match self.noise_model {
            NoiseModel::AdvancedLIGO => aligo_noise_asd(frequency),
            NoiseModel::InitialLIGO => initial_ligo_noise_asd(frequency),
            NoiseModel::EinsteinTelescope => einstein_telescope_noise_asd(frequency),
        }
    }

    /// Antenna pattern functions (F+, Fx)
    /// Simplified: assumes optimal orientation
    pub fn antenna_pattern(&self, _theta: f64, _phi: f64, _psi: f64) -> (f64, f64) {
        // Simplified: return optimal response
        (1.0, 1.0)
    }

    /// Project waveform onto detector
    pub fn project_waveform(&self, waveform: &Waveform) -> Vec<f64> {
        let (f_plus, f_cross) = self.antenna_pattern(0.0, 0.0, 0.0);

        waveform.h_plus.iter().zip(&waveform.h_cross)
            .map(|(&hp, &hc)| f_plus * hp + f_cross * hc)
            .collect()
    }
}

/// Advanced LIGO noise curve (analytical approximation)
fn aligo_noise_asd(f: f64) -> f64 {
    if f < 10.0 || f > 5000.0 {
        return 1e-20; // Outside sensitive band
    }

    // Simplified analytical fit
    let f0 = 215.0; // Minimum noise frequency
    let seismic = 1e-23 * (f0 / f).powf(4.0);
    let shot = 1e-23 * (f / f0);
    let thermal = 1e-23;

    (seismic + thermal + shot).sqrt()
}

/// Initial LIGO noise curve
fn initial_ligo_noise_asd(f: f64) -> f64 {
    if f < 40.0 || f > 2000.0 {
        return 1e-20;
    }

    let f0 = 150.0;
    let seismic = 1e-22 * (f0 / f).powf(4.0);
    let shot = 1e-22 * (f / f0);

    (seismic + shot).sqrt()
}

/// Einstein Telescope (3rd generation) noise curve
fn einstein_telescope_noise_asd(f: f64) -> f64 {
    if f < 1.0 || f > 10000.0 {
        return 1e-26;
    }

    let f0 = 100.0;
    let seismic = 1e-25 * (f0 / f).powf(5.0);
    let shot = 1e-25 * (f / f0).powf(0.5);

    (seismic + shot).sqrt()
}

/// Detector network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectorNetwork {
    pub detectors: Vec<LIGODetector>,
}

impl DetectorNetwork {
    /// LIGO Hanford-Livingston network
    pub fn ligo_hl() -> Self {
        Self {
            detectors: vec![LIGODetector::hanford(), LIGODetector::livingston()],
        }
    }

    /// Combined network sensitivity
    pub fn network_sensitivity(&self, frequency: f64) -> f64 {
        // Quadrature sum of sensitivities
        let sum_sq: f64 = self.detectors.iter()
            .map(|det| {
                let asd = det.noise_asd(frequency);
                1.0 / (asd * asd)
            })
            .sum();

        1.0 / sum_sq.sqrt()
    }
}

