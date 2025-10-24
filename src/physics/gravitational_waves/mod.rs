//! Gravitational Waves and LIGO Detection
//!
//! Simulates gravitational wave signals from:
//! - Binary black hole mergers (BBH)
//! - Binary neutron star mergers (BNS)
//! - Black hole - neutron star mergers (BHNS)
//!
//! Includes:
//! - Post-Newtonian waveform approximations
//! - LIGO detector response functions
//! - Signal-to-noise ratio calculations
//! - Matched filtering
//!
//! **References:**
//! - Abbott et al. (2016) - GW150914 observation
//! - Blanchet (2014) - Post-Newtonian formalism

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub mod waveforms;
pub mod detector;
pub mod snr;

pub use waveforms::*;
pub use detector::*;
pub use snr::*;

/// Physical constants
pub const C: f64 = 299792458.0; // m/s
pub const G: f64 = 6.67430e-11; // m³/(kg·s²)
pub const M_SUN: f64 = 1.989e30; // kg
pub const PC: f64 = 3.086e16; // parsec in meters

/// Binary system configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinarySystem {
    pub mass1: f64,        // Primary mass (solar masses)
    pub mass2: f64,        // Secondary mass (solar masses)
    pub distance: f64,     // Luminosity distance (Mpc)
    pub inclination: f64,  // Orbital inclination (radians)
    pub coalescence_phase: f64, // Phase at coalescence (radians)
    pub system_type: BinaryType,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum BinaryType {
    BinaryBlackHole,
    BinaryNeutronStar,
    BlackHoleNeutronStar,
}

impl BinarySystem {
    /// Total mass in solar masses
    pub fn total_mass(&self) -> f64 {
        self.mass1 + self.mass2
    }

    /// Reduced mass μ = m₁m₂/(m₁+m₂) in solar masses
    pub fn reduced_mass(&self) -> f64 {
        (self.mass1 * self.mass2) / self.total_mass()
    }

    /// Chirp mass M = (m₁m₂)^(3/5) / (m₁+m₂)^(1/5)
    pub fn chirp_mass(&self) -> f64 {
        let eta = self.symmetric_mass_ratio();
        self.total_mass() * eta.powf(3.0 / 5.0)
    }

    /// Symmetric mass ratio η = μ/M = m₁m₂/(m₁+m₂)²
    pub fn symmetric_mass_ratio(&self) -> f64 {
        (self.mass1 * self.mass2) / (self.total_mass().powi(2))
    }

    /// Schwarzschild radius of total mass
    pub fn schwarzschild_radius(&self) -> f64 {
        2.0 * G * (self.total_mass() * M_SUN) / (C * C)
    }

    /// ISCO frequency for non-spinning binary
    pub fn isco_frequency(&self) -> f64 {
        let r_s = self.schwarzschild_radius();
        let r_isco = 6.0 * r_s;
        C.powi(3) / (2.0 * PI * (G * self.total_mass() * M_SUN).sqrt() * r_isco.sqrt())
    }

    /// Time to coalescence from frequency (PN approximation)
    pub fn time_to_coalescence(&self, frequency: f64) -> f64 {
        let m_chirp = self.chirp_mass() * M_SUN;
        let coeff = 5.0 / (256.0 * (PI * frequency).powf(8.0 / 3.0));
        coeff * (G * m_chirp / C.powi(3)).powf(-5.0 / 3.0)
    }
}

/// Gravitational wave strain timeseries
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Waveform {
    pub times: Vec<f64>,      // Time samples (seconds)
    pub h_plus: Vec<f64>,     // Plus polarization
    pub h_cross: Vec<f64>,    // Cross polarization
    pub frequency: Vec<f64>,  // Instantaneous frequency
}

impl Waveform {
    pub fn len(&self) -> usize {
        self.times.len()
    }

    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// Peak strain amplitude
    pub fn peak_strain(&self) -> f64 {
        let max_plus = self.h_plus.iter().map(|h| h.abs()).fold(0.0, f64::max);
        let max_cross = self.h_cross.iter().map(|h| h.abs()).fold(0.0, f64::max);
        max_plus.max(max_cross)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_mass_parameters() {
        let binary = BinarySystem {
            mass1: 36.0,
            mass2: 29.0,
            distance: 410.0, // Mpc (GW150914-like)
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let total = binary.total_mass();
        assert_eq!(total, 65.0);

        let chirp = binary.chirp_mass();
        assert!(chirp > 25.0 && chirp < 35.0); // ~30 M_sun for GW150914

        let eta = binary.symmetric_mass_ratio();
        assert!(eta > 0.2 && eta < 0.26); // Should be ~0.247
    }

    #[test]
    fn test_isco_frequency() {
        let binary = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 100.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let f_isco = binary.isco_frequency();

        // ISCO frequency should be positive and finite
        assert!(f_isco > 0.0);
        assert!(f_isco.is_finite());
    }

    #[test]
    fn test_time_to_coalescence() {
        let binary = BinarySystem {
            mass1: 1.4,
            mass2: 1.4,
            distance: 100.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryNeutronStar,
        };

        // Time to coalescence from 10 Hz
        let t_coal = binary.time_to_coalescence(10.0);

        // For BNS, should be minutes to hours at 10 Hz
        assert!(t_coal > 60.0); // More than 1 minute
    }

    #[test]
    fn test_binary_types() {
        let bbh = BinarySystem {
            mass1: 30.0,
            mass2: 30.0,
            distance: 410.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryBlackHole,
        };

        let bns = BinarySystem {
            mass1: 1.4,
            mass2: 1.4,
            distance: 40.0,
            inclination: 0.0,
            coalescence_phase: 0.0,
            system_type: BinaryType::BinaryNeutronStar,
        };

        assert_eq!(bbh.system_type, BinaryType::BinaryBlackHole);
        assert_eq!(bns.system_type, BinaryType::BinaryNeutronStar);
    }
}
