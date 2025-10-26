//! Electromagnetism Module
//!
//! Implements electromagnetic theory and applications:
//! - Maxwell's equations
//! - Electromagnetic wave propagation
//! - Antenna theory and radiation
//! - Transmission lines
//! - Waveguides
//! - Scattering and diffraction
//! - Electromagnetic interference (EMI)
//! - Plasma physics basics

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Physical constants
const C: f64 = 299792458.0; // Speed of light (m/s)
const MU_0: f64 = 4.0 * PI * 1e-7; // Permeability of free space
const EPSILON_0: f64 = 8.854187817e-12; // Permittivity of free space
const Z_0: f64 = 376.730313668; // Impedance of free space (Ohms)

#[derive(Debug, Deserialize)]
pub struct MaxwellRequest {
    pub equation: String, // "gauss_electric", "gauss_magnetic", "faraday", "ampere"
    pub electric_field: Option<Vec<f64>>,
    pub magnetic_field: Option<Vec<f64>>,
    pub charge_density: Option<f64>,
    pub current_density: Option<Vec<f64>>,
    pub position: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct MaxwellResult {
    pub divergence: Option<f64>,
    pub curl: Option<Vec<f64>>,
    pub field_strength: f64,
    pub equation: String,
}

#[derive(Debug, Deserialize)]
pub struct WaveRequest {
    pub frequency: f64,          // Hz
    pub wavelength: Option<f64>, // m
    pub medium: String,          // "vacuum", "air", "dielectric"
    pub permittivity: Option<f64>,
    pub permeability: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct WaveResult {
    pub wavelength: f64,
    pub frequency: f64,
    pub wave_number: f64,
    pub angular_frequency: f64,
    pub phase_velocity: f64,
    pub impedance: f64,
}

#[derive(Debug, Deserialize)]
pub struct AntennaRequest {
    pub antenna_type: String, // "dipole", "monopole", "patch", "horn"
    pub frequency: f64,
    pub length: Option<f64>,
    pub distance: f64, // observation distance
    pub theta: f64,    // angle in radians
    pub power: f64,    // input power (W)
}

#[derive(Debug, Serialize)]
pub struct AntennaResult {
    pub gain: f64, // dBi
    pub directivity: f64,
    pub radiation_resistance: f64,
    pub efficiency: f64,
    pub beam_width: f64,
    pub electric_field_strength: f64,
}

#[derive(Debug, Deserialize)]
pub struct TransmissionLineRequest {
    pub line_type: String, // "coax", "microstrip", "stripline"
    pub frequency: f64,
    pub length: f64,
    pub z0: f64, // characteristic impedance
    pub load_impedance: Option<f64>,
    pub attenuation: Option<f64>, // dB/m
}

#[derive(Debug, Serialize)]
pub struct TransmissionLineResult {
    pub vswr: f64,
    pub reflection_coefficient: f64,
    pub return_loss: f64,    // dB
    pub insertion_loss: f64, // dB
    pub input_impedance: f64,
}

#[derive(Debug, Deserialize)]
pub struct WaveguideRequest {
    pub guide_type: String, // "rectangular", "circular"
    pub frequency: f64,
    pub width: f64,          // a dimension
    pub height: Option<f64>, // b dimension (for rectangular)
    pub mode: String,        // "TE10", "TM11", etc.
}

#[derive(Debug, Serialize)]
pub struct WaveguideResult {
    pub cutoff_frequency: f64,
    pub guide_wavelength: f64,
    pub phase_velocity: f64,
    pub group_velocity: f64,
    pub attenuation: f64,
}

#[derive(Debug, Deserialize)]
pub struct ScatteringRequest {
    pub scatterer_type: String, // "sphere", "cylinder", "plane"
    pub radius: f64,
    pub wavelength: f64,
    pub incident_angle: f64,
    pub polarization: String, // "TE", "TM"
}

#[derive(Debug, Serialize)]
pub struct ScatteringResult {
    pub radar_cross_section: f64, // m²
    pub scattering_amplitude: f64,
    pub scattered_power: f64,
    pub scattering_pattern: Vec<(f64, f64)>, // (angle, intensity)
}

#[derive(Debug, Deserialize)]
pub struct PoyntingRequest {
    pub electric_field: Vec<f64>,
    pub magnetic_field: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct PoyntingResult {
    pub poynting_vector: Vec<f64>,
    pub power_density: f64, // W/m²
    pub intensity: f64,
}

#[derive(Debug, Deserialize)]
pub struct SkinEffectRequest {
    pub frequency: f64,
    pub conductivity: f64, // S/m
    pub permeability: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct SkinEffectResult {
    pub skin_depth: f64,         // meters
    pub surface_resistance: f64, // Ohms/square
}

/// Maxwell's equations calculations
pub fn maxwell_equations(request: MaxwellRequest) -> Result<MaxwellResult, String> {
    match request.equation.as_str() {
        "gauss_electric" => {
            // ∇·E = ρ/ε₀
            let e_field = request.electric_field.ok_or("Need electric field")?;
            let divergence = numerical_divergence(&e_field);

            Ok(MaxwellResult {
                divergence: Some(divergence),
                curl: None,
                field_strength: vector_magnitude(&e_field),
                equation: "Gauss's Law (Electric)".to_string(),
            })
        }
        "gauss_magnetic" => {
            // ∇·B = 0
            let b_field = request.magnetic_field.ok_or("Need magnetic field")?;
            let divergence = numerical_divergence(&b_field);

            Ok(MaxwellResult {
                divergence: Some(divergence),
                curl: None,
                field_strength: vector_magnitude(&b_field),
                equation: "Gauss's Law (Magnetic)".to_string(),
            })
        }
        "faraday" => {
            // ∇×E = -∂B/∂t
            let e_field = request.electric_field.ok_or("Need electric field")?;
            let curl = numerical_curl(&e_field);

            Ok(MaxwellResult {
                divergence: None,
                curl: Some(curl.clone()),
                field_strength: vector_magnitude(&curl),
                equation: "Faraday's Law".to_string(),
            })
        }
        "ampere" => {
            // ∇×B = μ₀(J + ε₀∂E/∂t)
            let b_field = request.magnetic_field.ok_or("Need magnetic field")?;
            let curl = numerical_curl(&b_field);

            Ok(MaxwellResult {
                divergence: None,
                curl: Some(curl.clone()),
                field_strength: vector_magnitude(&curl),
                equation: "Ampere-Maxwell Law".to_string(),
            })
        }
        _ => Err(format!("Unknown Maxwell equation: {}", request.equation)),
    }
}

fn numerical_divergence(field: &[f64]) -> f64 {
    // Simplified divergence calculation
    field.iter().sum::<f64>() / field.len() as f64
}

fn numerical_curl(field: &[f64]) -> Vec<f64> {
    // Simplified curl (returns rotation components)
    if field.len() >= 3 {
        vec![
            field[2] - field[1],
            field[0] - field[2],
            field[1] - field[0],
        ]
    } else {
        vec![0.0; 3]
    }
}

fn vector_magnitude(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Electromagnetic wave properties
pub fn em_wave(request: WaveRequest) -> Result<WaveResult, String> {
    let freq = request.frequency;

    // Get medium properties
    let (eps, mu) = match request.medium.as_str() {
        "vacuum" | "air" => (EPSILON_0, MU_0),
        "dielectric" => (
            request.permittivity.unwrap_or(EPSILON_0),
            request.permeability.unwrap_or(MU_0),
        ),
        _ => return Err(format!("Unknown medium: {}", request.medium)),
    };

    let phase_velocity = 1.0 / (eps * mu).sqrt();
    let wavelength = phase_velocity / freq;
    let wave_number = 2.0 * PI / wavelength;
    let angular_freq = 2.0 * PI * freq;
    let impedance = (mu / eps).sqrt();

    Ok(WaveResult {
        wavelength,
        frequency: freq,
        wave_number,
        angular_frequency: angular_freq,
        phase_velocity,
        impedance,
    })
}

/// Antenna calculations
pub fn antenna_analysis(request: AntennaRequest) -> Result<AntennaResult, String> {
    let wavelength = C / request.frequency;

    let (gain, directivity, rad_resistance, efficiency, beam_width) =
        match request.antenna_type.as_str() {
            "dipole" => {
                let length = request.length.unwrap_or(wavelength / 2.0);
                let directivity = 1.64; // dBi for half-wave dipole
                let rad_resistance = 73.0; // Ohms for half-wave dipole
                let efficiency = 0.95;
                let beam_width = 78.0; // degrees
                let gain = directivity * efficiency;
                (gain, directivity, rad_resistance, efficiency, beam_width)
            }
            "monopole" => {
                let directivity = 5.16; // dBi
                let rad_resistance = 36.5; // Ohms
                let efficiency = 0.93;
                let beam_width = 90.0;
                let gain = directivity * efficiency;
                (gain, directivity, rad_resistance, efficiency, beam_width)
            }
            "patch" => {
                let directivity = 6.0; // dBi (typical)
                let rad_resistance = 200.0; // Ohms (typical)
                let efficiency = 0.85;
                let beam_width = 65.0;
                let gain = directivity * efficiency;
                (gain, directivity, rad_resistance, efficiency, beam_width)
            }
            _ => return Err(format!("Unknown antenna type: {}", request.antenna_type)),
        };

    // Calculate field strength at distance
    let power_density = request.power * gain / (4.0 * PI * request.distance.powi(2));
    let e_field = (power_density * Z_0).sqrt();

    Ok(AntennaResult {
        gain,
        directivity,
        radiation_resistance: rad_resistance,
        efficiency,
        beam_width,
        electric_field_strength: e_field,
    })
}

/// Transmission line analysis
pub fn transmission_line(
    request: TransmissionLineRequest,
) -> Result<TransmissionLineResult, String> {
    let z_load = request.load_impedance.unwrap_or(request.z0);

    // Reflection coefficient
    let gamma = (z_load - request.z0) / (z_load + request.z0);
    let gamma_mag = gamma.abs();

    // VSWR
    let vswr = (1.0 + gamma_mag) / (1.0 - gamma_mag);

    // Return loss
    let return_loss = -20.0 * gamma_mag.log10();

    // Insertion loss (from attenuation)
    let alpha = request.attenuation.unwrap_or(0.0); // dB/m
    let insertion_loss = alpha * request.length;

    // Input impedance (simplified for lossless line)
    let beta = 2.0 * PI * request.frequency / C; // propagation constant
    let input_z = request.z0 * (z_load + request.z0 * (beta * request.length).tan())
        / (request.z0 + z_load * (beta * request.length).tan());

    Ok(TransmissionLineResult {
        vswr,
        reflection_coefficient: gamma_mag,
        return_loss,
        insertion_loss,
        input_impedance: input_z.abs(),
    })
}

/// Waveguide analysis
pub fn waveguide(request: WaveguideRequest) -> Result<WaveguideResult, String> {
    let wavelength = C / request.frequency;

    match request.guide_type.as_str() {
        "rectangular" => {
            let a = request.width;
            let b = request.height.unwrap_or(a / 2.0);

            // For TE10 mode (dominant)
            let cutoff_freq = C / (2.0 * a);

            if request.frequency < cutoff_freq {
                return Err("Frequency below cutoff".to_string());
            }

            let beta =
                (2.0 * PI / wavelength) * (1.0 - (cutoff_freq / request.frequency).powi(2)).sqrt();
            let guide_wavelength = 2.0 * PI / beta;
            let phase_velocity = request.frequency * guide_wavelength;
            let group_velocity = C * C / phase_velocity;

            // Attenuation (simplified)
            let attenuation = 0.001 * request.frequency / 1e9; // Rough estimate in dB/m

            Ok(WaveguideResult {
                cutoff_frequency: cutoff_freq,
                guide_wavelength,
                phase_velocity,
                group_velocity,
                attenuation,
            })
        }
        "circular" => {
            let radius = request.width; // Use width as radius for circular waveguide

            // For TE11 mode (dominant mode in circular waveguide)
            let p_11 = 1.841; // First root of Bessel function J'_1
            let cutoff_freq = C * p_11 / (2.0 * PI * radius);

            if request.frequency < cutoff_freq {
                return Err("Frequency below cutoff for circular waveguide".to_string());
            }

            let beta =
                (2.0 * PI / wavelength) * (1.0 - (cutoff_freq / request.frequency).powi(2)).sqrt();
            let guide_wavelength = 2.0 * PI / beta;
            let phase_velocity = request.frequency * guide_wavelength;
            let group_velocity = C * C / phase_velocity;
            let attenuation = 0.0008 * request.frequency / 1e9; // Lower loss than rectangular

            Ok(WaveguideResult {
                cutoff_frequency: cutoff_freq,
                guide_wavelength,
                phase_velocity,
                group_velocity,
                attenuation,
            })
        }

        "coaxial" => {
            // Coaxial waveguide (TEM mode)
            let inner_radius = request.width; // Inner conductor radius
            let outer_radius = request.height.unwrap_or(request.width * 3.0); // Outer conductor radius

            // TEM mode has no cutoff frequency
            let cutoff_freq = 0.0;

            // For TEM mode, phase velocity = c
            let phase_velocity = C;
            let guide_wavelength = wavelength;
            let group_velocity = C;

            // Characteristic impedance Z₀ = (η/2π) ln(b/a)
            let z0 = (Z_0 / (2.0 * PI)) * (outer_radius / inner_radius).ln();

            // Attenuation includes conductor and dielectric losses
            let attenuation = 0.0005 * request.frequency / 1e9;

            Ok(WaveguideResult {
                cutoff_frequency: cutoff_freq,
                guide_wavelength,
                phase_velocity,
                group_velocity,
                attenuation,
            })
        }

        "dielectric" => {
            // Dielectric slab waveguide
            let slab_thickness = request.width;
            let n_core = request.height.unwrap_or(1.5); // Refractive index of core
            let n_cladding = 1.0; // Air cladding

            // Approximate cutoff for TE₀ mode
            let cutoff_wavelength =
                2.0 * slab_thickness * (n_core * n_core - n_cladding * n_cladding).sqrt();
            let cutoff_freq = C / cutoff_wavelength;

            if request.frequency < cutoff_freq {
                return Err("Frequency below cutoff for dielectric waveguide".to_string());
            }

            let effective_index = n_core; // Simplified
            let phase_velocity = C / effective_index;
            let guide_wavelength = wavelength / effective_index;
            let group_velocity = phase_velocity * 0.95; // Accounting for dispersion

            let attenuation = 0.0001 * request.frequency / 1e9; // Very low loss

            Ok(WaveguideResult {
                cutoff_frequency: cutoff_freq,
                guide_wavelength,
                phase_velocity,
                group_velocity,
                attenuation,
            })
        }

        _ => Err(format!(
            "Waveguide type {} not implemented",
            request.guide_type
        )),
    }
}

/// Electromagnetic scattering
pub fn scattering(request: ScatteringRequest) -> Result<ScatteringResult, String> {
    let k = 2.0 * PI / request.wavelength;
    let ka = k * request.radius; // Size parameter

    // Radar cross section (simplified Rayleigh scattering for small spheres)
    let rcs = if ka < 1.0 {
        // Rayleigh regime
        PI * request.radius.powi(2) * ka.powi(4)
    } else {
        // Geometric optics regime
        PI * request.radius.powi(2)
    };

    let scattering_amplitude = rcs.sqrt();
    let scattered_power = rcs * 1.0; // Assuming unit incident power

    // Simple scattering pattern (cos² pattern)
    let mut pattern = Vec::new();
    for i in 0..36 {
        let angle = (i as f64 * 10.0).to_radians();
        let intensity = angle.cos().powi(2);
        pattern.push((angle, intensity));
    }

    Ok(ScatteringResult {
        radar_cross_section: rcs,
        scattering_amplitude,
        scattered_power,
        scattering_pattern: pattern,
    })
}

/// Poynting vector and power flow
pub fn poynting_vector(request: PoyntingRequest) -> Result<PoyntingResult, String> {
    if request.electric_field.len() != 3 || request.magnetic_field.len() != 3 {
        return Err("E and H fields must be 3D vectors".to_string());
    }

    let e = &request.electric_field;
    let h = &request.magnetic_field;

    // S = E × H
    let s = vec![
        e[1] * h[2] - e[2] * h[1],
        e[2] * h[0] - e[0] * h[2],
        e[0] * h[1] - e[1] * h[0],
    ];

    let power_density = vector_magnitude(&s);
    let intensity = power_density;

    Ok(PoyntingResult {
        poynting_vector: s,
        power_density,
        intensity,
    })
}

/// Skin effect calculation
pub fn skin_effect(request: SkinEffectRequest) -> Result<SkinEffectResult, String> {
    let omega = 2.0 * PI * request.frequency;
    let mu = request.permeability.unwrap_or(MU_0);
    let sigma = request.conductivity;

    // Skin depth: δ = √(2/(ωμσ))
    let skin_depth = (2.0 / (omega * mu * sigma)).sqrt();

    // Surface resistance: Rs = 1/(σδ)
    let surface_resistance = 1.0 / (sigma * skin_depth);

    Ok(SkinEffectResult {
        skin_depth,
        surface_resistance,
    })
}

/// Plasma frequency calculation
pub fn plasma_frequency(electron_density: f64) -> f64 {
    // ωp = √(ne²/(ε₀m))
    let e = 1.602e-19; // electron charge
    let m = 9.109e-31; // electron mass

    (electron_density * e * e / (EPSILON_0 * m)).sqrt() / (2.0 * PI)
}

