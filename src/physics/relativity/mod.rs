//! Relativity Module
//!
//! Implements special and general relativity calculations:
//! - Special Relativity: Lorentz transformations, time dilation, length contraction
//! - Relativistic kinematics: energy, momentum, velocity addition
//! - General Relativity: Schwarzschild metric, gravitational effects
//! - Black hole physics: event horizon, Hawking radiation

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Physical constants
const C: f64 = 299792458.0; // Speed of light (m/s)
const G: f64 = 6.67430e-11; // Gravitational constant (m³/kg/s²)
const H_BAR: f64 = 1.054571817e-34; // Reduced Planck constant (J·s)
const K_B: f64 = 1.380649e-23; // Boltzmann constant (J/K)

// ============================================================================
// SPECIAL RELATIVITY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct LorentzTransformRequest {
    pub velocity: f64,      // Relative velocity (m/s)
    pub position: Vec<f64>, // [x, y, z] (m)
    pub time: f64,          // Time (s)
}

#[derive(Debug, Serialize)]
pub struct LorentzTransformResult {
    pub gamma: f64,               // Lorentz factor
    pub position_prime: Vec<f64>, // Transformed position
    pub time_prime: f64,          // Transformed time
    pub beta: f64,                // v/c
}

#[derive(Debug, Deserialize)]
pub struct TimeDilationRequest {
    pub proper_time: f64, // Time in rest frame (s)
    pub velocity: f64,    // Velocity (m/s)
}

#[derive(Debug, Serialize)]
pub struct TimeDilationResult {
    pub gamma: f64,
    pub dilated_time: f64, // Time in moving frame (s)
    pub time_difference: f64,
}

#[derive(Debug, Deserialize)]
pub struct LengthContractionRequest {
    pub proper_length: f64, // Length in rest frame (m)
    pub velocity: f64,      // Velocity (m/s)
}

#[derive(Debug, Serialize)]
pub struct LengthContractionResult {
    pub gamma: f64,
    pub contracted_length: f64, // Length in moving frame (m)
    pub contraction_factor: f64,
}

#[derive(Debug, Deserialize)]
pub struct RelativisticEnergyRequest {
    pub mass: f64,     // Rest mass (kg)
    pub velocity: f64, // Velocity (m/s)
}

#[derive(Debug, Serialize)]
pub struct RelativisticEnergyResult {
    pub gamma: f64,
    pub rest_energy: f64,    // E₀ = mc² (J)
    pub kinetic_energy: f64, // K = (γ-1)mc² (J)
    pub total_energy: f64,   // E = γmc² (J)
    pub momentum: f64,       // p = γmv (kg·m/s)
}

#[derive(Debug, Deserialize)]
pub struct VelocityAdditionRequest {
    pub velocity1: f64, // First velocity (m/s)
    pub velocity2: f64, // Second velocity (m/s)
}

#[derive(Debug, Serialize)]
pub struct VelocityAdditionResult {
    pub relativistic_sum: f64, // Relativistic velocity addition
    pub classical_sum: f64,    // Classical (Galilean) sum
    pub difference: f64,
}

// ============================================================================
// GENERAL RELATIVITY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SchwarzschildRequest {
    pub mass: f64,   // Mass of object (kg)
    pub radius: f64, // Radial distance from center (m)
}

#[derive(Debug, Serialize)]
pub struct SchwarzschildResult {
    pub schwarzschild_radius: f64,   // Event horizon radius (m)
    pub metric_g00: f64,             // Time component of metric
    pub metric_grr: f64,             // Radial component of metric
    pub escape_velocity: f64,        // Escape velocity at radius (m/s)
    pub gravitational_redshift: f64, // z = sqrt(g00) - 1
}

#[derive(Debug, Deserialize)]
pub struct GravitationalTimeDilationRequest {
    pub mass: f64,        // Mass of gravitating body (kg)
    pub radius: f64,      // Distance from center (m)
    pub proper_time: f64, // Time measured at radius (s)
}

#[derive(Debug, Serialize)]
pub struct GravitationalTimeDilationResult {
    pub schwarzschild_radius: f64,
    pub time_dilation_factor: f64, // sqrt(1 - rs/r)
    pub coordinate_time: f64,      // Time at infinity
    pub time_difference: f64,
}

#[derive(Debug, Deserialize)]
pub struct OrbitalPrecessionRequest {
    pub mass: f64,            // Central mass (kg)
    pub semi_major_axis: f64, // Semi-major axis (m)
    pub eccentricity: f64,    // Orbital eccentricity
}

#[derive(Debug, Serialize)]
pub struct OrbitalPrecessionResult {
    pub precession_per_orbit: f64,   // Radians per orbit
    pub precession_per_century: f64, // Arcseconds per century
    pub orbital_period: f64,         // Seconds
}

#[derive(Debug, Deserialize)]
pub struct GravitationalLensingRequest {
    pub lens_mass: f64,        // Mass of lensing object (kg)
    pub impact_parameter: f64, // Closest approach distance (m)
}

#[derive(Debug, Serialize)]
pub struct GravitationalLensingResult {
    pub deflection_angle: f64, // Radians
    pub einstein_radius: f64,  // Characteristic lensing scale (m)
}

// ============================================================================
// BLACK HOLE PHYSICS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct BlackHoleRequest {
    pub mass: f64, // Black hole mass (kg)
}

#[derive(Debug, Serialize)]
pub struct BlackHoleResult {
    pub schwarzschild_radius: f64, // Event horizon (m)
    pub surface_gravity: f64,      // m/s² at horizon
    pub hawking_temperature: f64,  // Kelvin
    pub hawking_luminosity: f64,   // Watts
    pub evaporation_time: f64,     // Seconds
    pub entropy: f64,              // Bekenstein-Hawking entropy (J/K)
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/// Calculate Lorentz factor γ = 1/√(1 - v²/c²)
fn lorentz_factor(velocity: f64) -> Result<f64, String> {
    if velocity.abs() >= C {
        return Err("Velocity must be less than speed of light".to_string());
    }
    let beta = velocity / C;
    let gamma = 1.0 / (1.0 - beta * beta).sqrt();
    Ok(gamma)
}

// ============================================================================
// SPECIAL RELATIVITY IMPLEMENTATIONS
// ============================================================================

/// Lorentz transformation (boost in x-direction)
pub fn lorentz_transform(
    request: LorentzTransformRequest,
) -> Result<LorentzTransformResult, String> {
    let gamma = lorentz_factor(request.velocity)?;
    let beta = request.velocity / C;

    if request.position.len() != 3 {
        return Err("Position must be 3D vector [x, y, z]".to_string());
    }

    let x = request.position[0];
    let y = request.position[1];
    let z = request.position[2];
    let t = request.time;

    // Lorentz transformation (boost in x-direction)
    let x_prime = gamma * (x - beta * C * t);
    let y_prime = y;
    let z_prime = z;
    let t_prime = gamma * (t - beta * x / C);

    Ok(LorentzTransformResult {
        gamma,
        position_prime: vec![x_prime, y_prime, z_prime],
        time_prime: t_prime,
        beta,
    })
}

/// Time dilation: Δt' = γΔt
pub fn time_dilation(request: TimeDilationRequest) -> Result<TimeDilationResult, String> {
    let gamma = lorentz_factor(request.velocity)?;
    let dilated_time = gamma * request.proper_time;

    Ok(TimeDilationResult {
        gamma,
        dilated_time,
        time_difference: dilated_time - request.proper_time,
    })
}

/// Length contraction: L = L₀/γ
pub fn length_contraction(
    request: LengthContractionRequest,
) -> Result<LengthContractionResult, String> {
    let gamma = lorentz_factor(request.velocity)?;
    let contracted_length = request.proper_length / gamma;

    Ok(LengthContractionResult {
        gamma,
        contracted_length,
        contraction_factor: 1.0 / gamma,
    })
}

/// Relativistic energy and momentum
pub fn relativistic_energy(
    request: RelativisticEnergyRequest,
) -> Result<RelativisticEnergyResult, String> {
    let gamma = lorentz_factor(request.velocity)?;
    let rest_energy = request.mass * C * C;
    let total_energy = gamma * rest_energy;
    let kinetic_energy = (gamma - 1.0) * rest_energy;
    let momentum = gamma * request.mass * request.velocity;

    Ok(RelativisticEnergyResult {
        gamma,
        rest_energy,
        kinetic_energy,
        total_energy,
        momentum,
    })
}

/// Relativistic velocity addition: v = (v₁ + v₂)/(1 + v₁v₂/c²)
pub fn velocity_addition(
    request: VelocityAdditionRequest,
) -> Result<VelocityAdditionResult, String> {
    let v1 = request.velocity1;
    let v2 = request.velocity2;

    if v1.abs() >= C || v2.abs() >= C {
        return Err("Velocities must be less than speed of light".to_string());
    }

    let relativistic_sum = (v1 + v2) / (1.0 + v1 * v2 / (C * C));
    let classical_sum = v1 + v2;

    Ok(VelocityAdditionResult {
        relativistic_sum,
        classical_sum,
        difference: classical_sum - relativistic_sum,
    })
}

// ============================================================================
// GENERAL RELATIVITY IMPLEMENTATIONS
// ============================================================================

/// Schwarzschild metric and properties
pub fn schwarzschild_metric(request: SchwarzschildRequest) -> Result<SchwarzschildResult, String> {
    let rs = 2.0 * G * request.mass / (C * C); // Schwarzschild radius

    if request.radius <= rs {
        return Err("Radius must be outside event horizon".to_string());
    }

    let r = request.radius;

    // Metric components (Schwarzschild coordinates)
    let g00 = -(1.0 - rs / r); // Time component
    let grr = 1.0 / (1.0 - rs / r); // Radial component

    // Escape velocity at radius r
    let escape_velocity = C * (rs / r).sqrt();

    // Gravitational redshift: z = 1/sqrt(1 - rs/r) - 1
    let gravitational_redshift = 1.0 / (1.0 - rs / r).sqrt() - 1.0;

    Ok(SchwarzschildResult {
        schwarzschild_radius: rs,
        metric_g00: g00,
        metric_grr: grr,
        escape_velocity,
        gravitational_redshift,
    })
}

/// Gravitational time dilation
pub fn gravitational_time_dilation(
    request: GravitationalTimeDilationRequest,
) -> Result<GravitationalTimeDilationResult, String> {
    let rs = 2.0 * G * request.mass / (C * C);

    if request.radius <= rs {
        return Err("Radius must be outside event horizon".to_string());
    }

    let time_dilation_factor = (1.0 - rs / request.radius).sqrt();
    let coordinate_time = request.proper_time / time_dilation_factor;

    Ok(GravitationalTimeDilationResult {
        schwarzschild_radius: rs,
        time_dilation_factor,
        coordinate_time,
        time_difference: coordinate_time - request.proper_time,
    })
}

/// Orbital precession (perihelion advance)
pub fn orbital_precession(
    request: OrbitalPrecessionRequest,
) -> Result<OrbitalPrecessionResult, String> {
    let a = request.semi_major_axis;
    let e = request.eccentricity;
    let m = request.mass;

    // Orbital period (Kepler's third law)
    let period = 2.0 * PI * (a.powi(3) / (G * m)).sqrt();

    // Precession per orbit (Einstein's formula)
    let precession_per_orbit = 6.0 * PI * G * m / (a * C * C * (1.0 - e * e));

    // Convert to arcseconds per century
    let orbits_per_century = 100.0 * 365.25 * 24.0 * 3600.0 / period;
    let precession_per_century = precession_per_orbit * orbits_per_century * 206265.0; // radians to arcseconds

    Ok(OrbitalPrecessionResult {
        precession_per_orbit,
        precession_per_century,
        orbital_period: period,
    })
}

/// Gravitational lensing deflection angle
pub fn gravitational_lensing(
    request: GravitationalLensingRequest,
) -> Result<GravitationalLensingResult, String> {
    let rs = 2.0 * G * request.lens_mass / (C * C);
    let b = request.impact_parameter;

    // Einstein deflection angle: θ = 4GM/(bc²)
    let deflection_angle = 2.0 * rs / b;

    // Einstein radius (characteristic scale)
    let einstein_radius = (2.0 * rs * b).sqrt();

    Ok(GravitationalLensingResult {
        deflection_angle,
        einstein_radius,
    })
}

// ============================================================================
// BLACK HOLE PHYSICS IMPLEMENTATIONS
// ============================================================================

/// Black hole thermodynamics and Hawking radiation
pub fn black_hole_properties(request: BlackHoleRequest) -> Result<BlackHoleResult, String> {
    let m = request.mass;
    let rs = 2.0 * G * m / (C * C); // Schwarzschild radius

    // Surface gravity at event horizon
    let surface_gravity = C.powi(4) / (4.0 * G * m);

    // Hawking temperature: T = ℏc³/(8πGMk_B)
    let hawking_temp = H_BAR * C.powi(3) / (8.0 * PI * G * m * K_B);

    // Stefan-Boltzmann law for luminosity
    let sigma_sb = 5.670374419e-8; // Stefan-Boltzmann constant
    let area = 4.0 * PI * rs.powi(2);
    let hawking_luminosity = sigma_sb * area * hawking_temp.powi(4);

    // Evaporation time: t ~ M³
    let evaporation_time = 2.1e67 * (m / 1e30).powi(3); // Approximate formula

    // Bekenstein-Hawking entropy: S = kc³A/(4ℏG)
    let entropy = K_B * C.powi(3) * area / (4.0 * H_BAR * G);

    Ok(BlackHoleResult {
        schwarzschild_radius: rs,
        surface_gravity,
        hawking_temperature: hawking_temp,
        hawking_luminosity,
        evaporation_time,
        entropy,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lorentz_factor() {
        // At 0.6c, γ = 1.25
        let gamma = lorentz_factor(0.6 * C).unwrap();
        assert!((gamma - 1.25).abs() < 0.01);
    }

    #[test]
    fn test_time_dilation() {
        let result = time_dilation(TimeDilationRequest {
            proper_time: 1.0,
            velocity: 0.6 * C,
        })
        .unwrap();

        assert!((result.gamma - 1.25).abs() < 0.01);
        assert!((result.dilated_time - 1.25).abs() < 0.01);
    }

    #[test]
    fn test_length_contraction() {
        let result = length_contraction(LengthContractionRequest {
            proper_length: 10.0,
            velocity: 0.6 * C,
        })
        .unwrap();

        assert!((result.contracted_length - 8.0).abs() < 0.1);
    }

    #[test]
    fn test_relativistic_energy() {
        let result = relativistic_energy(RelativisticEnergyRequest {
            mass: 1.0, // 1 kg
            velocity: 0.6 * C,
        })
        .unwrap();

        let expected_rest = C * C;
        assert!((result.rest_energy - expected_rest).abs() < 1e10);
        assert!(result.total_energy > result.rest_energy);
    }

    #[test]
    fn test_velocity_addition() {
        let result = velocity_addition(VelocityAdditionRequest {
            velocity1: 0.6 * C,
            velocity2: 0.6 * C,
        })
        .unwrap();

        // Should be less than c
        assert!(result.relativistic_sum < C);
        // Should be ~0.882c
        assert!((result.relativistic_sum / C - 0.882).abs() < 0.01);
    }

    #[test]
    fn test_schwarzschild_radius() {
        // Solar mass
        let solar_mass = 1.989e30; // kg
        let result = schwarzschild_metric(SchwarzschildRequest {
            mass: solar_mass,
            radius: 1e10, // 10,000 km
        })
        .unwrap();

        // Schwarzschild radius of sun ~2.95 km
        assert!((result.schwarzschild_radius - 2950.0).abs() < 100.0);
    }

    #[test]
    fn test_black_hole_properties() {
        let solar_mass = 1.989e30; // kg
        let result = black_hole_properties(BlackHoleRequest { mass: solar_mass }).unwrap();

        assert!(result.schwarzschild_radius > 0.0);
        assert!(result.hawking_temperature > 0.0);
        assert!(result.entropy > 0.0);
    }
}
