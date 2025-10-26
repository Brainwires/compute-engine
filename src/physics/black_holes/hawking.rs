//! Hawking Radiation and Quantum Effects

use super::BlackHoleConfig;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HawkingRadiation {
    pub temperature: f64,        // Kelvin
    pub luminosity: f64,          // Watts
    pub peak_wavelength: f64,     // meters
    pub evaporation_time: f64,    // seconds
}

/// Calculate Hawking radiation properties
pub fn hawking_radiation(config: &BlackHoleConfig) -> HawkingRadiation {
    let temp = config.hawking_temperature();

    // Stefan-Boltzmann law: L = σ A T⁴
    const STEFAN_BOLTZMANN: f64 = 5.670374419e-8; // W/(m²·K⁴)
    let area = config.horizon_area();
    let luminosity = STEFAN_BOLTZMANN * area * temp.powi(4);

    // Wien's law: λ_peak = b/T
    const WIEN_CONSTANT: f64 = 2.897771955e-3; // m·K
    let peak_wavelength = WIEN_CONSTANT / temp;

    let evaporation_time = config.evaporation_time();

    HawkingRadiation {
        temperature: temp,
        luminosity,
        peak_wavelength,
        evaporation_time,
    }
}

/// Mass loss rate: dM/dt = -ℏc⁴/(15360πG²M²)
pub fn mass_loss_rate(config: &BlackHoleConfig) -> f64 {
    const HBAR: f64 = 1.054571817e-34;
    let c4 = super::C2 * super::C2;
    let g2 = super::G * super::G;
    let m2 = config.mass * config.mass;

    -(HBAR * c4) / (15360.0 * std::f64::consts::PI * g2 * m2)
}

