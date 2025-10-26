/**
 * Geophysics Module
 *
 * Implements Earth science formulas:
 * - Seismology (magnitude scales, wave velocities)
 * - Atmospheric physics (hydrostatic equation, humidity)
 * - Geochemistry (radioactive dating)
 * - Planetary science (escape velocity, Roche limit)
 */
use serde::{Deserialize, Serialize};
use std::f64::consts::E;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeophysicsInput {
    pub category: GeophysicsCategory,
    pub parameters: GeophysicsParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeophysicsCategory {
    Seismology,
    Atmosphere,
    RadiometricDating,
    PlanetaryScience,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeophysicsParams {
    // Seismology
    pub amplitude: Option<f64>,      // μm
    pub distance: Option<f64>,       // km
    pub seismic_moment: Option<f64>, // N·m
    pub energy: Option<f64>,         // J
    pub magnitude: Option<f64>,      // Richter or Mw

    // Atmospheric
    pub pressure: Option<f64>,          // Pa or hPa
    pub temperature: Option<f64>,       // K or °C
    pub altitude: Option<f64>,          // m
    pub density: Option<f64>,           // kg/m³
    pub relative_humidity: Option<f64>, // %
    pub vapor_pressure: Option<f64>,    // Pa

    // Radioactive dating
    pub parent_isotope: Option<f64>,    // amount
    pub daughter_isotope: Option<f64>,  // amount
    pub half_life: Option<f64>,         // years
    pub isotope_system: Option<String>, // "C14", "U238Pb206", "K40Ar40"

    // Planetary
    pub mass_primary: Option<f64>,      // kg
    pub mass_secondary: Option<f64>,    // kg
    pub radius_primary: Option<f64>,    // m
    pub density_primary: Option<f64>,   // kg/m³
    pub density_secondary: Option<f64>, // kg/m³
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeophysicsResult {
    pub value: f64,
    pub unit: String,
    pub formula_used: String,
    pub uncertainty: Option<f64>,
    pub interpretation: String,
    pub additional_data: Option<std::collections::HashMap<String, f64>>,
}

const G: f64 = 6.67430e-11; // Gravitational constant (m³/(kg·s²))
const G_ACCEL: f64 = 9.80665; // Standard gravity (m/s²)
const R_GAS: f64 = 287.05; // Specific gas constant for dry air (J/(kg·K))

pub fn calculate_geophysics(input: GeophysicsInput) -> Result<GeophysicsResult, String> {
    match input.category {
        GeophysicsCategory::Seismology => calculate_seismology(&input.parameters),
        GeophysicsCategory::Atmosphere => calculate_atmosphere(&input.parameters),
        GeophysicsCategory::RadiometricDating => calculate_dating(&input.parameters),
        GeophysicsCategory::PlanetaryScience => calculate_planetary(&input.parameters),
    }
}

/// Earthquake magnitude calculations
fn calculate_seismology(params: &GeophysicsParams) -> Result<GeophysicsResult, String> {
    let mut additional = std::collections::HashMap::new();

    if let Some(m0) = params.seismic_moment {
        // Moment magnitude: Mw = (2/3)·log10(M0) - 6.07
        // M0 in N·m
        let mw = (2.0 / 3.0) * m0.log10() - 6.07;

        // Energy release: log10(E) = 1.5·Mw + 4.8
        let log_energy = 1.5 * mw + 4.8;
        let energy = 10_f64.powf(log_energy);

        additional.insert("moment_magnitude_mw".to_string(), mw);
        additional.insert("energy_joules".to_string(), energy);
        additional.insert("energy_tnt_equivalent_tons".to_string(), energy / 4.184e9);

        let severity = match mw {
            m if m < 3.0 => "Micro (rarely felt)",
            m if m < 4.0 => "Minor (often felt, no damage)",
            m if m < 5.0 => "Light (some damage)",
            m if m < 6.0 => "Moderate (considerable damage)",
            m if m < 7.0 => "Strong (major damage)",
            m if m < 8.0 => "Major (serious damage over large areas)",
            _ => "Great (devastating damage)",
        };

        return Ok(GeophysicsResult {
            value: mw,
            unit: "Mw".to_string(),
            formula_used: "Moment magnitude: Mw = (2/3)·log10(M0) - 6.07".to_string(),
            uncertainty: Some(0.1),
            interpretation: format!("{} earthquake", severity),
            additional_data: Some(additional),
        });
    }

    if let Some(energy) = params.energy {
        // Richter scale from energy: M = (2/3)·log10(E) - 2.9
        let magnitude = (2.0 / 3.0) * energy.log10() - 2.9;

        additional.insert("richter_magnitude".to_string(), magnitude);
        additional.insert("tnt_equivalent_tons".to_string(), energy / 4.184e9);

        return Ok(GeophysicsResult {
            value: magnitude,
            unit: "ML (Richter)".to_string(),
            formula_used: "Richter: M = (2/3)·log10(E) - 2.9".to_string(),
            uncertainty: Some(0.2),
            interpretation: format!("Local magnitude {:.1}", magnitude),
            additional_data: Some(additional),
        });
    }

    Err("Seismic moment or energy required".to_string())
}

/// Atmospheric physics calculations
fn calculate_atmosphere(params: &GeophysicsParams) -> Result<GeophysicsResult, String> {
    let mut additional = std::collections::HashMap::new();

    // Hydrostatic equation: dP/dz = -ρg
    if let (Some(p0), Some(t), Some(h)) = (params.pressure, params.temperature, params.altitude) {
        // Barometric formula: P(h) = P0·exp(-g·M·h/(R·T))
        // Simplified for small h and isothermal atmosphere
        let temp_k = if t < 200.0 { t + 273.15 } else { t }; // Convert if in Celsius

        // Scale height: H = RT/g
        let scale_height = (R_GAS * temp_k) / G_ACCEL;

        let pressure_at_h = p0 * E.powf(-h / scale_height);

        additional.insert("scale_height_m".to_string(), scale_height);
        additional.insert(
            "pressure_drop_percent".to_string(),
            (p0 - pressure_at_h) / p0 * 100.0,
        );

        return Ok(GeophysicsResult {
            value: pressure_at_h,
            unit: "Pa".to_string(),
            formula_used: "Hydrostatic: P(h) = P0·exp(-h/H)".to_string(),
            uncertainty: Some(pressure_at_h * 0.02),
            interpretation: format!("Pressure at {:.0}m altitude: {:.0} Pa", h, pressure_at_h),
            additional_data: Some(additional),
        });
    }

    // Clausius-Clapeyron for humidity
    if let (Some(t), Some(rh)) = (params.temperature, params.relative_humidity) {
        let temp_k = if t < 200.0 { t + 273.15 } else { t };

        // Saturation vapor pressure (Magnus formula): es = 6.112·exp(17.67·T/(T+243.5))
        // T in Celsius
        let temp_c = temp_k - 273.15;
        let e_sat = 6.112 * E.powf(17.67 * temp_c / (temp_c + 243.5)) * 100.0; // Convert to Pa

        let e_actual = (rh / 100.0) * e_sat;

        // Dew point approximation
        let dew_point =
            243.5 * (e_actual / 100.0 / 6.112).ln() / (17.67 - (e_actual / 100.0 / 6.112).ln());

        additional.insert("saturation_vapor_pressure_pa".to_string(), e_sat);
        additional.insert("actual_vapor_pressure_pa".to_string(), e_actual);
        additional.insert("dew_point_celsius".to_string(), dew_point);

        return Ok(GeophysicsResult {
            value: dew_point,
            unit: "°C".to_string(),
            formula_used: "Clausius-Clapeyron (Magnus approximation)".to_string(),
            uncertainty: Some(0.5),
            interpretation: format!("Dew point: {:.1}°C at {:.0}% RH", dew_point, rh),
            additional_data: Some(additional),
        });
    }

    Err("Insufficient atmospheric parameters".to_string())
}

/// Radioactive dating
fn calculate_dating(params: &GeophysicsParams) -> Result<GeophysicsResult, String> {
    let system = params
        .isotope_system
        .as_ref()
        .ok_or("Isotope system required (C14, U238Pb206, K40Ar40)")?;

    let half_life = match system.as_str() {
        "C14" => 5730.0,        // Carbon-14 (years)
        "U238Pb206" => 4.468e9, // Uranium-238 (years)
        "K40Ar40" => 1.25e9,    // Potassium-40 (years)
        "Rb87Sr87" => 4.88e10,  // Rubidium-87 (years)
        _ => params
            .half_life
            .ok_or("Unknown isotope system or half_life not provided")?,
    };

    let decay_constant = 0.693147 / half_life; // λ = ln(2)/t½

    let mut additional = std::collections::HashMap::new();

    if let (Some(parent), Some(daughter)) = (params.parent_isotope, params.daughter_isotope) {
        // Age calculation: t = (1/λ)·ln(1 + D/P)
        // where D = daughter atoms, P = parent atoms
        let age = (1.0 / decay_constant) * (1.0 + daughter / parent).ln();

        // Calculate uncertainty (simplified, assumes 2% analytical uncertainty)
        let uncertainty = age * 0.02;

        additional.insert("half_life_years".to_string(), half_life);
        additional.insert("decay_constant".to_string(), decay_constant);
        additional.insert("daughter_parent_ratio".to_string(), daughter / parent);

        let time_scale = if age < 1e6 {
            format!("{:.0} years", age)
        } else if age < 1e9 {
            format!("{:.2} million years", age / 1e6)
        } else {
            format!("{:.3} billion years", age / 1e9)
        };

        return Ok(GeophysicsResult {
            value: age,
            unit: "years".to_string(),
            formula_used: format!("{} dating: t = (1/λ)·ln(1 + D/P)", system),
            uncertainty: Some(uncertainty),
            interpretation: format!("Age: {}", time_scale),
            additional_data: Some(additional),
        });
    }

    Err("Parent and daughter isotope amounts required".to_string())
}

/// Planetary science calculations
fn calculate_planetary(params: &GeophysicsParams) -> Result<GeophysicsResult, String> {
    let mut additional = std::collections::HashMap::new();

    // Escape velocity: v = √(2GM/r)
    if let (Some(mass), Some(radius)) = (params.mass_primary, params.radius_primary) {
        let v_escape = (2.0 * G * mass / radius).sqrt();

        additional.insert("escape_velocity_m_per_s".to_string(), v_escape);
        additional.insert("escape_velocity_km_per_s".to_string(), v_escape / 1000.0);

        // Calculate orbital velocity at surface: v_orb = √(GM/r)
        let v_orbital = (G * mass / radius).sqrt();
        additional.insert("orbital_velocity_surface".to_string(), v_orbital);

        return Ok(GeophysicsResult {
            value: v_escape / 1000.0,
            unit: "km/s".to_string(),
            formula_used: "Escape velocity: v = √(2GM/r)".to_string(),
            uncertainty: Some(v_escape / 1000.0 * 0.01),
            interpretation: format!("Escape velocity: {:.2} km/s", v_escape / 1000.0),
            additional_data: Some(additional),
        });
    }

    // Roche limit: d = 2.46·R·(ρM/ρm)^(1/3)
    if let (Some(r_primary), Some(rho_primary), Some(rho_secondary)) = (
        params.radius_primary,
        params.density_primary,
        params.density_secondary,
    ) {
        let roche_rigid = 2.46 * r_primary * (rho_primary / rho_secondary).powf(1.0 / 3.0);
        let roche_fluid = 2.44 * r_primary * (rho_primary / rho_secondary).powf(1.0 / 3.0);

        additional.insert("roche_limit_rigid".to_string(), roche_rigid);
        additional.insert("roche_limit_fluid".to_string(), roche_fluid);
        additional.insert(
            "roche_limit_primary_radii".to_string(),
            roche_fluid / r_primary,
        );

        return Ok(GeophysicsResult {
            value: roche_fluid / 1000.0,
            unit: "km".to_string(),
            formula_used: "Roche limit: d = 2.44·R·(ρM/ρm)^(1/3)".to_string(),
            uncertainty: Some(roche_fluid / 1000.0 * 0.05),
            interpretation: format!(
                "Roche limit: {:.0} km ({:.2} primary radii)",
                roche_fluid / 1000.0,
                roche_fluid / r_primary
            ),
            additional_data: Some(additional),
        });
    }

    Err("Insufficient planetary parameters".to_string())
}


impl Default for GeophysicsParams {
    fn default() -> Self {
        Self {
            amplitude: None,
            distance: None,
            seismic_moment: None,
            energy: None,
            magnitude: None,
            pressure: None,
            temperature: None,
            altitude: None,
            density: None,
            relative_humidity: None,
            vapor_pressure: None,
            parent_isotope: None,
            daughter_isotope: None,
            half_life: None,
            isotope_system: None,
            mass_primary: None,
            mass_secondary: None,
            radius_primary: None,
            density_primary: None,
            density_secondary: None,
        }
    }
}

#[cfg(test)]
#[path = "../../tests/unit/geophysics_tests.rs"]
mod tests;
