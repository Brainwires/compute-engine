//! Statistical Physics Module
//!
//! Implements statistical mechanics and thermodynamics:
//! - Partition Functions (canonical, grand canonical, microcanonical)
//! - Statistical Distributions (Maxwell-Boltzmann, Fermi-Dirac, Bose-Einstein)
//! - Thermodynamic Properties (chemical potential, fugacity)
//! - Phase Transitions and Critical Phenomena

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Physical constants
const K_B: f64 = 1.380649e-23; // Boltzmann constant (J/K)
const H: f64 = 6.62607015e-34; // Planck constant (J·s)
const H_BAR: f64 = 1.054571817e-34; // Reduced Planck constant (J·s)

// ============================================================================
// PARTITION FUNCTIONS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct PartitionFunctionRequest {
    pub ensemble: String, // "canonical", "grand_canonical", "microcanonical"
    pub temperature: f64, // Kelvin
    pub volume: Option<f64>, // m³
    pub num_particles: Option<usize>,
    pub energy_levels: Option<Vec<f64>>, // Energy levels (J)
    pub degeneracies: Option<Vec<usize>>, // Degeneracy of each level
    pub chemical_potential: Option<f64>, // For grand canonical (J)
}

#[derive(Debug, Serialize)]
pub struct PartitionFunctionResult {
    pub partition_function: f64,
    pub ensemble_type: String,
    pub helmholtz_free_energy: Option<f64>, // F = -kT ln(Z)
    pub internal_energy: Option<f64>, // U = -∂ln(Z)/∂β
    pub entropy: Option<f64>, // S = k ln(Z) + U/T
    pub pressure: Option<f64>, // P = kT ∂ln(Z)/∂V
}

#[derive(Debug, Deserialize)]
pub struct CanonicalPartitionRequest {
    pub temperature: f64,
    pub energy_levels: Vec<f64>,
    pub degeneracies: Vec<usize>,
}

#[derive(Debug, Serialize)]
pub struct CanonicalPartitionResult {
    pub partition_function: f64,
    pub helmholtz_free_energy: f64,
    pub internal_energy: f64,
    pub entropy: f64,
    pub heat_capacity: f64,
}

#[derive(Debug, Deserialize)]
pub struct GrandCanonicalRequest {
    pub temperature: f64,
    pub volume: f64,
    pub chemical_potential: f64,
    pub particle_type: String, // "fermion", "boson", "classical"
}

#[derive(Debug, Serialize)]
pub struct GrandCanonicalResult {
    pub grand_potential: f64, // Ω = -kT ln(Ξ)
    pub average_particle_number: f64,
    pub pressure: f64,
    pub entropy: f64,
}

// ============================================================================
// STATISTICAL DISTRIBUTIONS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct MaxwellBoltzmannRequest {
    pub temperature: f64, // Kelvin
    pub mass: f64, // Particle mass (kg)
    pub velocity: Option<f64>, // Speed (m/s) - if None, return average
}

#[derive(Debug, Serialize)]
pub struct MaxwellBoltzmannResult {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub probability_density: Option<f64>, // f(v) at given velocity
    pub average_speed: f64, // <v> = √(8kT/πm)
    pub rms_speed: f64, // v_rms = √(3kT/m)
    pub most_probable_speed: f64, // v_p = √(2kT/m)
}

#[derive(Debug, Deserialize)]
pub struct FermiDiracRequest {
    pub energy: f64, // Energy level (J)
    pub temperature: f64, // Kelvin
    pub chemical_potential: f64, // Fermi energy (J)
}

#[derive(Debug, Serialize)]
pub struct FermiDiracResult {
    pub occupation_probability: f64, // f(E) = 1/(exp((E-μ)/kT) + 1)
    pub fermi_temperature: f64, // T_F = E_F/k
    pub fermi_energy: f64, // Chemical potential
}

#[derive(Debug, Deserialize)]
pub struct BoseEinsteinRequest {
    pub energy: f64, // Energy level (J)
    pub temperature: f64, // Kelvin
    pub chemical_potential: f64, // Must be < energy (J)
}

#[derive(Debug, Serialize)]
pub struct BoseEinsteinResult {
    pub occupation_number: f64, // n(E) = 1/(exp((E-μ)/kT) - 1)
    pub bose_condensation_temp: Option<f64>, // Critical temperature
}

// ============================================================================
// THERMODYNAMIC PROPERTIES
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct ChemicalPotentialRequest {
    pub temperature: f64,
    pub pressure: f64,
    pub particle_density: f64, // particles/m³
    pub particle_type: String, // "fermion", "boson", "classical"
    pub mass: f64, // Particle mass (kg)
}

#[derive(Debug, Serialize)]
pub struct ChemicalPotentialResult {
    pub chemical_potential: f64, // μ (J)
    pub thermal_wavelength: f64, // λ = h/√(2πmkT)
    pub fugacity: f64, // z = exp(μ/kT)
}

#[derive(Debug, Deserialize)]
pub struct FugacityRequest {
    pub temperature: f64,
    pub pressure: f64,
    pub particle_density: f64,
}

#[derive(Debug, Serialize)]
pub struct FugacityResult {
    pub fugacity: f64, // z = exp(μ/kT)
    pub fugacity_coefficient: f64, // φ = f/P
    pub chemical_potential: f64,
}

// ============================================================================
// PHASE TRANSITIONS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct PhaseTransitionRequest {
    pub model: String, // "ising", "van_der_waals", "landau"
    pub temperature: f64,
    pub critical_temperature: f64,
    pub parameters: std::collections::HashMap<String, f64>,
}

#[derive(Debug, Serialize)]
pub struct PhaseTransitionResult {
    pub order_parameter: f64,
    pub correlation_length: f64,
    pub susceptibility: f64,
    pub phase: String, // "ordered", "disordered", "critical"
}

#[derive(Debug, Deserialize)]
pub struct CriticalPhenomenaRequest {
    pub temperature: f64,
    pub critical_temperature: f64,
    pub model: String, // "ising_2d", "ising_3d", "mean_field"
}

#[derive(Debug, Serialize)]
pub struct CriticalPhenomenaResult {
    pub critical_exponents: std::collections::HashMap<String, f64>,
    pub reduced_temperature: f64, // t = (T - T_c)/T_c
    pub order_parameter: f64,
    pub correlation_length: f64,
}

// ============================================================================
// IMPLEMENTATIONS
// ============================================================================

/// General partition function calculator
pub fn partition_function(request: PartitionFunctionRequest) -> Result<PartitionFunctionResult, String> {
    let t = request.temperature;
    let beta = 1.0 / (K_B * t);

    match request.ensemble.as_str() {
        "canonical" => {
            let energy_levels = request.energy_levels.ok_or("Energy levels required for canonical ensemble")?;
            let degeneracies = request.degeneracies.unwrap_or_else(|| vec![1; energy_levels.len()]);

            let mut z = 0.0;
            let mut u = 0.0; // Internal energy

            for (i, &energy) in energy_levels.iter().enumerate() {
                let g = degeneracies.get(i).unwrap_or(&1);
                let boltzmann_factor = (*g as f64) * (-beta * energy).exp();
                z += boltzmann_factor;
                u += energy * boltzmann_factor;
            }

            u /= z; // <E> = Σ E_i P_i

            let f = -K_B * t * z.ln(); // Helmholtz free energy
            let s = K_B * z.ln() + u / t; // Entropy

            Ok(PartitionFunctionResult {
                partition_function: z,
                ensemble_type: "canonical".to_string(),
                helmholtz_free_energy: Some(f),
                internal_energy: Some(u),
                entropy: Some(s),
                pressure: None,
            })
        },

        "grand_canonical" => {
            Err("Use grand_canonical_partition for grand canonical ensemble".to_string())
        },

        "microcanonical" => {
            let num_states = request.energy_levels.map(|e| e.len()).unwrap_or(1);
            let s = K_B * (num_states as f64).ln();

            Ok(PartitionFunctionResult {
                partition_function: num_states as f64,
                ensemble_type: "microcanonical".to_string(),
                helmholtz_free_energy: None,
                internal_energy: None,
                entropy: Some(s),
                pressure: None,
            })
        },

        _ => Err(format!("Unknown ensemble type: {}", request.ensemble)),
    }
}

/// Canonical partition function (detailed)
pub fn canonical_partition(request: CanonicalPartitionRequest) -> Result<CanonicalPartitionResult, String> {
    let t = request.temperature;
    let beta = 1.0 / (K_B * t);

    let mut z = 0.0;
    let mut u = 0.0;
    let mut u_squared = 0.0;

    for (i, &energy) in request.energy_levels.iter().enumerate() {
        let g = request.degeneracies.get(i).unwrap_or(&1);
        let boltzmann_factor = (*g as f64) * (-beta * energy).exp();
        z += boltzmann_factor;
        u += energy * boltzmann_factor;
        u_squared += energy * energy * boltzmann_factor;
    }

    u /= z;
    u_squared /= z;

    let variance = u_squared - u * u;
    let heat_capacity = variance / (K_B * t * t);

    let f = -K_B * t * z.ln();
    let s = K_B * z.ln() + u / t;

    Ok(CanonicalPartitionResult {
        partition_function: z,
        helmholtz_free_energy: f,
        internal_energy: u,
        entropy: s,
        heat_capacity,
    })
}

/// Grand canonical partition function
pub fn grand_canonical_partition(request: GrandCanonicalRequest) -> Result<GrandCanonicalResult, String> {
    let t = request.temperature;
    let v = request.volume;
    let mu = request.chemical_potential;
    let beta = 1.0 / (K_B * t);

    // Thermal wavelength
    let m = 1.0e-26; // Typical atomic mass (kg) - would be parameter in real implementation
    let lambda = H / (2.0 * PI * m * K_B * t).sqrt();

    match request.particle_type.as_str() {
        "classical" => {
            // Classical ideal gas
            let fugacity = (beta * mu).exp();
            let n_avg = fugacity * v / (lambda.powi(3));
            let omega = -K_B * t * fugacity * v / (lambda.powi(3));
            let pressure = K_B * t * n_avg / v;
            let entropy = n_avg * K_B * (5.0/2.0 + (v / (n_avg * lambda.powi(3))).ln());

            Ok(GrandCanonicalResult {
                grand_potential: omega,
                average_particle_number: n_avg,
                pressure,
                entropy,
            })
        },

        "fermion" | "boson" => {
            // Simplified quantum gas
            let fugacity = (beta * mu).exp();
            let n_avg = v / lambda.powi(3) * fugacity; // Approximate
            let pressure = K_B * t * n_avg / v;
            let entropy = n_avg * K_B * 3.0/2.0;

            Ok(GrandCanonicalResult {
                grand_potential: -pressure * v,
                average_particle_number: n_avg,
                pressure,
                entropy,
            })
        },

        _ => Err(format!("Unknown particle type: {}", request.particle_type)),
    }
}

/// Maxwell-Boltzmann velocity distribution
pub fn maxwell_boltzmann(request: MaxwellBoltzmannRequest) -> Result<MaxwellBoltzmannResult, String> {
    let t = request.temperature;
    let m = request.mass;

    // Characteristic speeds
    let v_avg = (8.0 * K_B * t / (PI * m)).sqrt(); // Average speed
    let v_rms = (3.0 * K_B * t / m).sqrt(); // RMS speed
    let v_p = (2.0 * K_B * t / m).sqrt(); // Most probable speed

    let probability_density = if let Some(v) = request.velocity {
        // f(v) = 4π(m/2πkT)^(3/2) v² exp(-mv²/2kT)
        let normalization = 4.0 * PI * (m / (2.0 * PI * K_B * t)).powf(1.5);
        let prob = normalization * v * v * (- m * v * v / (2.0 * K_B * t)).exp();
        Some(prob)
    } else {
        None
    };

    Ok(MaxwellBoltzmannResult {
        probability_density,
        average_speed: v_avg,
        rms_speed: v_rms,
        most_probable_speed: v_p,
    })
}

/// Fermi-Dirac distribution
pub fn fermi_dirac(request: FermiDiracRequest) -> Result<FermiDiracResult, String> {
    let e = request.energy;
    let t = request.temperature;
    let mu = request.chemical_potential;

    // f(E) = 1 / (exp((E-μ)/kT) + 1)
    let occupation = 1.0 / (((e - mu) / (K_B * t)).exp() + 1.0);

    let fermi_temp = mu / K_B;

    Ok(FermiDiracResult {
        occupation_probability: occupation,
        fermi_temperature: fermi_temp,
        fermi_energy: mu,
    })
}

/// Bose-Einstein distribution
pub fn bose_einstein(request: BoseEinsteinRequest) -> Result<BoseEinsteinResult, String> {
    let e = request.energy;
    let t = request.temperature;
    let mu = request.chemical_potential;

    if mu >= e {
        return Err("Chemical potential must be less than energy for bosons".to_string());
    }

    // n(E) = 1 / (exp((E-μ)/kT) - 1)
    let occupation = 1.0 / (((e - mu) / (K_B * t)).exp() - 1.0);

    // Bose condensation temperature (simplified)
    let t_c = if mu < 0.0 {
        Some((-mu / K_B) * 0.527) // Approximate for 3D
    } else {
        None
    };

    Ok(BoseEinsteinResult {
        occupation_number: occupation,
        bose_condensation_temp: t_c,
    })
}

/// Chemical potential calculation
pub fn chemical_potential(request: ChemicalPotentialRequest) -> Result<ChemicalPotentialResult, String> {
    let t = request.temperature;
    let n = request.particle_density;
    let m = request.mass;

    // Thermal de Broglie wavelength
    let lambda = H / (2.0 * PI * m * K_B * t).sqrt();

    match request.particle_type.as_str() {
        "classical" => {
            // μ = kT ln(nλ³)
            let mu = K_B * t * (n * lambda.powi(3)).ln();
            let fugacity = (mu / (K_B * t)).exp();

            Ok(ChemicalPotentialResult {
                chemical_potential: mu,
                thermal_wavelength: lambda,
                fugacity,
            })
        },

        "fermion" => {
            // Approximate Fermi energy for degenerate gas
            let e_f = (H_BAR * H_BAR / (2.0 * m)) * (3.0 * PI * PI * n).powf(2.0/3.0);
            let fugacity = (e_f / (K_B * t)).exp();

            Ok(ChemicalPotentialResult {
                chemical_potential: e_f,
                thermal_wavelength: lambda,
                fugacity,
            })
        },

        "boson" => {
            // For bosons, μ < 0 typically
            let mu = -K_B * t * 0.1; // Simplified
            let fugacity = (mu / (K_B * t)).exp();

            Ok(ChemicalPotentialResult {
                chemical_potential: mu,
                thermal_wavelength: lambda,
                fugacity,
            })
        },

        _ => Err(format!("Unknown particle type: {}", request.particle_type)),
    }
}

/// Fugacity coefficient
pub fn fugacity_coefficient(request: FugacityRequest) -> Result<FugacityResult, String> {
    let t = request.temperature;
    let p = request.pressure;
    let n = request.particle_density;

    // For ideal gas: μ = kT ln(P/P₀)
    let p0 = 101325.0; // Standard pressure (Pa)
    let mu = K_B * t * (p / p0).ln();

    let fugacity = (mu / (K_B * t)).exp();
    let fugacity_coeff = fugacity / p;

    Ok(FugacityResult {
        fugacity,
        fugacity_coefficient: fugacity_coeff,
        chemical_potential: mu,
    })
}

/// Phase transition analysis
pub fn phase_transition(request: PhaseTransitionRequest) -> Result<PhaseTransitionResult, String> {
    let t = request.temperature;
    let t_c = request.critical_temperature;
    let reduced_t = (t - t_c) / t_c;

    match request.model.as_str() {
        "ising" => {
            // Ising model (mean field approximation)
            let order_param = if t < t_c {
                (1.0 - t / t_c).powf(0.5) // β = 0.5 for mean field
            } else {
                0.0
            };

            let correlation_length = if reduced_t.abs() > 1e-10 {
                1.0 / reduced_t.abs().powf(0.5) // ν = 0.5 for mean field
            } else {
                f64::INFINITY
            };

            let susceptibility = if reduced_t.abs() > 1e-10 {
                1.0 / reduced_t.abs() // γ = 1.0 for mean field
            } else {
                f64::INFINITY
            };

            let phase = if t < t_c {
                "ordered"
            } else if (t - t_c).abs() < 0.01 * t_c {
                "critical"
            } else {
                "disordered"
            }.to_string();

            Ok(PhaseTransitionResult {
                order_parameter: order_param,
                correlation_length,
                susceptibility,
                phase,
            })
        },

        "van_der_waals" => {
            // Van der Waals gas phase transition
            let pressure = request.parameters.get("pressure").ok_or("pressure required")?;
            let a = request.parameters.get("a").unwrap_or(&0.1); // Attraction parameter
            let b = request.parameters.get("b").unwrap_or(&0.01); // Excluded volume

            // Reduced variables
            let t_r = t / t_c;
            let p_c = 1.0 / (27.0 * b * b); // Critical pressure (simplified)
            let v_c = 3.0 * b; // Critical volume

            // Order parameter: density difference between liquid and gas
            let order_param = if t < t_c {
                (1.0 - t_r).powf(0.5)
            } else {
                0.0
            };

            let correlation_length = if (t - t_c).abs() > 1e-10 {
                1.0 / ((t - t_c) / t_c).abs().powf(0.5)
            } else {
                f64::INFINITY
            };

            let susceptibility = if (t - t_c).abs() > 1e-10 {
                1.0 / ((t - t_c) / t_c).abs()
            } else {
                f64::INFINITY
            };

            let phase = if t < t_c && *pressure < p_c {
                "two_phase"
            } else if t < t_c {
                "liquid"
            } else {
                "gas"
            }.to_string();

            Ok(PhaseTransitionResult {
                order_parameter: order_param,
                correlation_length,
                susceptibility,
                phase,
            })
        },

        "potts" => {
            // Potts model (generalization of Ising)
            let q = request.parameters.get("q").unwrap_or(&3.0).round() as usize; // Number of states

            // Mean field critical exponents
            let order_param = if t < t_c {
                ((q as f64 - 1.0) / q as f64 * (1.0 - t / t_c)).powf(0.5)
            } else {
                0.0
            };

            let correlation_length = if (t - t_c).abs() > 1e-10 {
                1.0 / ((t - t_c) / t_c).abs().powf(0.5)
            } else {
                f64::INFINITY
            };

            let susceptibility = if (t - t_c).abs() > 1e-10 {
                1.0 / ((t - t_c) / t_c).abs()
            } else {
                f64::INFINITY
            };

            let phase = if t < t_c {
                "ordered"
            } else if (t - t_c).abs() < 0.01 * t_c {
                "critical"
            } else {
                "disordered"
            }.to_string();

            Ok(PhaseTransitionResult {
                order_parameter: order_param,
                correlation_length,
                susceptibility,
                phase,
            })
        },

        "xy" => {
            // XY model (continuous spins in 2D)
            // Exhibits Kosterlitz-Thouless transition

            // For XY model, β_KT ≈ 0.5 (not mean field)
            let order_param = if t < t_c {
                (1.0 - t / t_c).powf(0.231) // β ≈ 1/8 for 2D XY
            } else {
                0.0
            };

            // Kosterlitz-Thouless transition has exponential correlation length
            let correlation_length = if t > t_c && (t - t_c) / t_c < 0.5 {
                ((t - t_c) / t_c).exp()
            } else if (t - t_c).abs() < 1e-10 {
                f64::INFINITY
            } else {
                1.0 / ((t - t_c) / t_c).abs().powf(0.67)
            };

            let susceptibility = if (t - t_c).abs() > 1e-10 {
                1.0 / ((t - t_c) / t_c).abs().powf(1.3)
            } else {
                f64::INFINITY
            };

            let phase = if t < t_c {
                "quasi_long_range_order"
            } else if (t - t_c).abs() < 0.01 * t_c {
                "kt_transition"
            } else {
                "disordered"
            }.to_string();

            Ok(PhaseTransitionResult {
                order_parameter: order_param,
                correlation_length,
                susceptibility,
                phase,
            })
        },

        "heisenberg" => {
            // Heisenberg model (3D continuous spins)
            // Critical exponents for 3D O(3) universality class
            let beta_exp = 0.365; // Order parameter exponent
            let nu_exp = 0.705; // Correlation length exponent
            let gamma_exp = 1.386; // Susceptibility exponent

            let order_param = if t < t_c {
                (1.0 - t / t_c).powf(beta_exp)
            } else {
                0.0
            };

            let correlation_length = if (t - t_c).abs() > 1e-10 {
                ((t - t_c) / t_c).abs().powf(-nu_exp)
            } else {
                f64::INFINITY
            };

            let susceptibility = if (t - t_c).abs() > 1e-10 {
                ((t - t_c) / t_c).abs().powf(-gamma_exp)
            } else {
                f64::INFINITY
            };

            let phase = if t < t_c {
                "ferromagnetic"
            } else if (t - t_c).abs() < 0.01 * t_c {
                "critical"
            } else {
                "paramagnetic"
            }.to_string();

            Ok(PhaseTransitionResult {
                order_parameter: order_param,
                correlation_length,
                susceptibility,
                phase,
            })
        },

        _ => Err(format!("Phase transition model {} not implemented", request.model)),
    }
}

/// Critical phenomena and exponents
pub fn critical_phenomena(request: CriticalPhenomenaRequest) -> Result<CriticalPhenomenaResult, String> {
    let t = request.temperature;
    let t_c = request.critical_temperature;
    let reduced_t = (t - t_c) / t_c;

    let mut exponents = std::collections::HashMap::new();

    match request.model.as_str() {
        "ising_2d" => {
            // Exact Onsager solution exponents
            exponents.insert("alpha".to_string(), 0.0); // Specific heat
            exponents.insert("beta".to_string(), 1.0/8.0); // Order parameter
            exponents.insert("gamma".to_string(), 7.0/4.0); // Susceptibility
            exponents.insert("delta".to_string(), 15.0); // Critical isotherm
            exponents.insert("nu".to_string(), 1.0); // Correlation length
            exponents.insert("eta".to_string(), 1.0/4.0); // Correlation function
        },

        "ising_3d" => {
            // Approximate 3D Ising exponents
            exponents.insert("alpha".to_string(), 0.110);
            exponents.insert("beta".to_string(), 0.326);
            exponents.insert("gamma".to_string(), 1.237);
            exponents.insert("delta".to_string(), 4.789);
            exponents.insert("nu".to_string(), 0.630);
            exponents.insert("eta".to_string(), 0.036);
        },

        "mean_field" => {
            // Classical mean field exponents
            exponents.insert("alpha".to_string(), 0.0);
            exponents.insert("beta".to_string(), 0.5);
            exponents.insert("gamma".to_string(), 1.0);
            exponents.insert("delta".to_string(), 3.0);
            exponents.insert("nu".to_string(), 0.5);
            exponents.insert("eta".to_string(), 0.0);
        },

        _ => return Err(format!("Unknown critical model: {}", request.model)),
    }

    let beta_exp = exponents.get("beta").unwrap_or(&0.5);
    let nu_exp = exponents.get("nu").unwrap_or(&0.5);

    let order_param = if reduced_t < 0.0 {
        (-reduced_t).powf(*beta_exp)
    } else {
        0.0
    };

    let correlation_length = if reduced_t.abs() > 1e-10 {
        reduced_t.abs().powf(-nu_exp)
    } else {
        f64::INFINITY
    };

    Ok(CriticalPhenomenaResult {
        critical_exponents: exponents,
        reduced_temperature: reduced_t,
        order_parameter: order_param,
        correlation_length,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maxwell_boltzmann() {
        let result = maxwell_boltzmann(MaxwellBoltzmannRequest {
            temperature: 300.0,
            mass: 4.65e-26, // Nitrogen molecule
            velocity: Some(500.0),
        }).unwrap();

        assert!(result.probability_density.is_some());
        assert!(result.average_speed > 0.0);
        assert!(result.rms_speed > result.average_speed);
    }

    #[test]
    fn test_fermi_dirac() {
        let result = fermi_dirac(FermiDiracRequest {
            energy: 1.0e-19, // 1 eV in Joules
            temperature: 300.0,
            chemical_potential: 0.8e-19,
        }).unwrap();

        assert!(result.occupation_probability > 0.0 && result.occupation_probability < 1.0);
    }

    #[test]
    fn test_canonical_partition() {
        let result = canonical_partition(CanonicalPartitionRequest {
            temperature: 300.0,
            energy_levels: vec![0.0, 1.0e-20, 2.0e-20],
            degeneracies: vec![1, 2, 1],
        }).unwrap();

        assert!(result.partition_function > 0.0);
        assert!(result.internal_energy >= 0.0);
    }
}
