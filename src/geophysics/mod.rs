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

#[cfg(test)]
mod tests {
    use super::*;

    // ===== SEISMOLOGY TESTS =====

    #[test]
    fn test_seismology_moment_magnitude_large_quake() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e20), // Large earthquake
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value > 6.0 && result.value < 8.0);
        assert_eq!(result.unit, "Mw");
        assert!(result.additional_data.is_some());
        let additional = result.additional_data.unwrap();
        assert!(additional.contains_key("moment_magnitude_mw"));
        assert!(additional.contains_key("energy_joules"));
        assert!(additional.contains_key("energy_tnt_equivalent_tons"));
    }

    #[test]
    fn test_seismology_moment_magnitude_micro() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e12), // Micro earthquake
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value < 3.0);
        assert!(result.interpretation.contains("Micro"));
    }

    #[test]
    fn test_seismology_moment_magnitude_moderate() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e18), // Moderate earthquake
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value >= 5.0 && result.value < 6.0);
        assert!(result.interpretation.contains("Moderate"));
    }

    #[test]
    fn test_seismology_moment_magnitude_great() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e23), // Great earthquake (like 2004 Indian Ocean)
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value >= 8.0);
        assert!(result.interpretation.contains("Great"));
    }

    #[test]
    fn test_seismology_richter_from_energy() {
        let params = GeophysicsParams {
            energy: Some(1e15), // ~7.1 magnitude
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        assert!(result.value > 6.5 && result.value < 7.5);
        assert_eq!(result.unit, "ML (Richter)");
        let additional = result.additional_data.unwrap();
        assert!(additional.contains_key("richter_magnitude"));
        assert!(additional.contains_key("tnt_equivalent_tons"));
    }

    #[test]
    fn test_seismology_energy_conversion() {
        let params = GeophysicsParams {
            seismic_moment: Some(1e20),
            ..Default::default()
        };

        let result = calculate_seismology(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let energy_joules = additional.get("energy_joules").unwrap();
        let tnt_tons = additional.get("energy_tnt_equivalent_tons").unwrap();

        // Verify conversion: 1 ton TNT = 4.184e9 J
        assert!((*energy_joules / 4.184e9 - tnt_tons).abs() < 1.0);
    }

    #[test]
    fn test_seismology_missing_parameters() {
        let params = GeophysicsParams::default();
        let result = calculate_seismology(&params);
        assert!(result.is_err());
    }

    // ===== ATMOSPHERE TESTS =====

    #[test]
    fn test_atmosphere_pressure_at_altitude() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),  // Sea level pressure (Pa)
            temperature: Some(288.15), // 15°C in Kelvin
            altitude: Some(1000.0),    // 1 km
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Pressure should decrease with altitude
        assert!(result.value < 101325.0);
        assert!(result.value > 85000.0); // Reasonable range
        assert_eq!(result.unit, "Pa");
        assert!(result.interpretation.contains("1000m altitude"));
    }

    #[test]
    fn test_atmosphere_scale_height() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(288.15),
            altitude: Some(8500.0), // Approximately one scale height
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let scale_height = additional.get("scale_height_m").unwrap();

        // Scale height for Earth's atmosphere ~8400m
        assert!((*scale_height - 8400.0).abs() < 500.0);

        // After one scale height, pressure should drop to ~1/e of initial
        let pressure_drop = additional.get("pressure_drop_percent").unwrap();
        assert!((*pressure_drop - 63.2).abs() < 5.0); // 1-1/e ≈ 63.2%
    }

    #[test]
    fn test_atmosphere_high_altitude() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(288.15),
            altitude: Some(10000.0), // 10 km (cruising altitude)
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Pressure at 10 km should be ~30% of sea level (~31000 Pa)
        assert!(result.value < 32000.0);
        assert!(result.value > 29000.0);
    }

    #[test]
    fn test_atmosphere_temperature_celsius_conversion() {
        let params_celsius = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(15.0), // Celsius
            altitude: Some(1000.0),
            ..Default::default()
        };

        let params_kelvin = GeophysicsParams {
            pressure: Some(101325.0),
            temperature: Some(288.15), // Kelvin equivalent
            altitude: Some(1000.0),
            ..Default::default()
        };

        let result_c = calculate_atmosphere(&params_celsius).unwrap();
        let result_k = calculate_atmosphere(&params_kelvin).unwrap();

        // Should give same result
        assert!((result_c.value - result_k.value).abs() < 10.0);
    }

    #[test]
    fn test_atmosphere_dew_point() {
        let params = GeophysicsParams {
            temperature: Some(20.0),       // 20°C
            relative_humidity: Some(60.0), // 60% RH
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Dew point should be below air temperature
        assert!(result.value < 20.0);
        assert!(result.value > 10.0); // Reasonable range for 60% RH
        assert_eq!(result.unit, "°C");
        assert!(result.interpretation.contains("Dew point"));
        assert!(result.interpretation.contains("60% RH"));
    }

    #[test]
    fn test_atmosphere_saturation_vapor_pressure() {
        let params = GeophysicsParams {
            temperature: Some(25.0),        // 25°C
            relative_humidity: Some(100.0), // Saturated
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let e_sat = additional.get("saturation_vapor_pressure_pa").unwrap();
        let e_actual = additional.get("actual_vapor_pressure_pa").unwrap();

        // At 100% RH, actual = saturation
        assert!((e_sat - e_actual).abs() < 1.0);

        // Saturation vapor pressure at 25°C ~3167 Pa
        assert!((*e_sat - 3167.0).abs() < 100.0);
    }

    #[test]
    fn test_atmosphere_dew_point_low_humidity() {
        let params = GeophysicsParams {
            temperature: Some(30.0),       // 30°C
            relative_humidity: Some(20.0), // Low humidity
            ..Default::default()
        };

        let result = calculate_atmosphere(&params).unwrap();
        // Low humidity means large difference between temp and dew point
        assert!(30.0 - result.value > 15.0);
    }

    #[test]
    fn test_atmosphere_missing_parameters() {
        let params = GeophysicsParams {
            pressure: Some(101325.0),
            // Missing temperature and altitude
            ..Default::default()
        };

        let result = calculate_atmosphere(&params);
        assert!(result.is_err());
    }

    // ===== RADIOMETRIC DATING TESTS =====

    #[test]
    fn test_carbon_dating_one_half_life() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // When P=D, age should be ~1 half-life
        assert!((result.value - 5730.0).abs() < 100.0);
        assert_eq!(result.unit, "years");
        assert!(result.interpretation.contains("years"));
    }

    #[test]
    fn test_carbon_dating_recent() {
        let params = GeophysicsParams {
            parent_isotope: Some(90.0),
            daughter_isotope: Some(10.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Young sample (less than 1 half-life)
        assert!(result.value < 5730.0);
        assert!(result.value > 0.0);
    }

    #[test]
    fn test_carbon_dating_old() {
        let params = GeophysicsParams {
            parent_isotope: Some(10.0),
            daughter_isotope: Some(90.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Old sample (multiple half-lives)
        assert!(result.value > 10000.0);
        assert!(result.value < 40000.0); // C14 limit
    }

    #[test]
    fn test_uranium_dating() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("U238Pb206".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // One half-life of U-238 is 4.468 billion years
        assert!((result.value - 4.468e9).abs() < 1e8);
        assert!(result.interpretation.contains("billion years"));
    }

    #[test]
    fn test_potassium_dating() {
        let params = GeophysicsParams {
            parent_isotope: Some(75.0),
            daughter_isotope: Some(25.0),
            isotope_system: Some("K40Ar40".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Younger than 1 half-life (1.25 billion years)
        assert!(result.value < 1.25e9);
        assert!(result.value > 0.0);
    }

    #[test]
    fn test_rubidium_dating() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("Rb87Sr87".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Rb-87 has very long half-life (48.8 billion years)
        assert!(result.value > 4e10);
        assert!(result.value < 5e10);
        assert!(result.interpretation.contains("billion years"));
    }

    #[test]
    fn test_dating_custom_half_life() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("Custom".to_string()),
            half_life: Some(10000.0),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        // Should use provided half-life
        assert!((result.value - 10000.0).abs() < 200.0);
    }

    #[test]
    fn test_dating_decay_constant() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let decay_const = additional.get("decay_constant").unwrap();
        let half_life = additional.get("half_life_years").unwrap();

        // λ = ln(2) / t½
        assert!((decay_const * half_life - 0.693147).abs() < 1e-5);
    }

    #[test]
    fn test_dating_daughter_parent_ratio() {
        let params = GeophysicsParams {
            parent_isotope: Some(25.0),
            daughter_isotope: Some(75.0), // 3:1 D/P ratio
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let ratio = additional.get("daughter_parent_ratio").unwrap();

        assert!((*ratio - 3.0).abs() < 0.01);
        // High D/P ratio means old sample
        assert!(result.value > 11000.0); // ~2 half-lives
    }

    #[test]
    fn test_dating_missing_isotope_system() {
        let params = GeophysicsParams {
            parent_isotope: Some(50.0),
            daughter_isotope: Some(50.0),
            ..Default::default()
        };

        let result = calculate_dating(&params);
        assert!(result.is_err());
    }

    #[test]
    fn test_dating_missing_isotope_amounts() {
        let params = GeophysicsParams {
            isotope_system: Some("C14".to_string()),
            ..Default::default()
        };

        let result = calculate_dating(&params);
        assert!(result.is_err());
    }

    // ===== PLANETARY SCIENCE TESTS =====

    #[test]
    fn test_escape_velocity_earth() {
        let params = GeophysicsParams {
            mass_primary: Some(5.972e24),  // Earth mass (kg)
            radius_primary: Some(6.371e6), // Earth radius (m)
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Earth's escape velocity ~11.2 km/s
        assert!((result.value - 11.2).abs() < 0.2);
        assert_eq!(result.unit, "km/s");
        assert!(result.interpretation.contains("Escape velocity"));
    }

    #[test]
    fn test_escape_velocity_mars() {
        let params = GeophysicsParams {
            mass_primary: Some(6.4171e23),  // Mars mass
            radius_primary: Some(3.3895e6), // Mars radius
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Mars escape velocity ~5.0 km/s
        assert!((result.value - 5.0).abs() < 0.3);
    }

    #[test]
    fn test_escape_velocity_jupiter() {
        let params = GeophysicsParams {
            mass_primary: Some(1.898e27),   // Jupiter mass
            radius_primary: Some(6.9911e7), // Jupiter radius
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Jupiter escape velocity ~59.5 km/s
        assert!((result.value - 59.5).abs() < 2.0);
    }

    #[test]
    fn test_escape_velocity_moon() {
        let params = GeophysicsParams {
            mass_primary: Some(7.342e22),  // Moon mass
            radius_primary: Some(1.737e6), // Moon radius
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Moon escape velocity ~2.38 km/s
        assert!((result.value - 2.38).abs() < 0.1);
    }

    #[test]
    fn test_orbital_velocity_surface() {
        let params = GeophysicsParams {
            mass_primary: Some(5.972e24),
            radius_primary: Some(6.371e6),
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let v_orbital = additional.get("orbital_velocity_surface").unwrap();

        // Earth surface orbital velocity ~7.9 km/s
        assert!((v_orbital / 1000.0 - 7.9).abs() < 0.2);

        // Escape velocity should be √2 times orbital velocity
        let v_escape = result.value * 1000.0; // Convert back to m/s
        assert!((v_escape / v_orbital - 2_f64.sqrt()).abs() < 0.01);
    }

    #[test]
    fn test_roche_limit_earth_moon() {
        let params = GeophysicsParams {
            radius_primary: Some(6.371e6),   // Earth radius
            density_primary: Some(5514.0),   // Earth density (kg/m³)
            density_secondary: Some(3344.0), // Moon density
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        // Earth-Moon Roche limit ~18,365 km (for fluid body)
        assert!((result.value - 18365.0).abs() < 500.0);
        assert_eq!(result.unit, "km");
    }

    #[test]
    fn test_roche_limit_saturn_rings() {
        let params = GeophysicsParams {
            radius_primary: Some(5.8232e7),  // Saturn radius
            density_primary: Some(687.0),    // Saturn density
            density_secondary: Some(1000.0), // Ice density
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let roche_radii = additional.get("roche_limit_primary_radii").unwrap();

        // Saturn's Roche limit is ~2.15 Saturn radii (density ratio < 1)
        assert!((*roche_radii - 2.15).abs() < 0.1);

        // Saturn's rings are inside Roche limit
        assert!(result.value > 100000.0); // >100,000 km
    }

    #[test]
    fn test_roche_limit_rigid_vs_fluid() {
        let params = GeophysicsParams {
            radius_primary: Some(6.371e6),
            density_primary: Some(5514.0),
            density_secondary: Some(3344.0),
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let rigid = additional.get("roche_limit_rigid").unwrap();
        let fluid = additional.get("roche_limit_fluid").unwrap();

        // Rigid body Roche limit (2.46) is slightly larger than fluid (2.44)
        assert!(rigid > fluid);
        assert!((rigid / fluid - 1.008).abs() < 0.01);
    }

    #[test]
    fn test_roche_limit_equal_densities() {
        let params = GeophysicsParams {
            radius_primary: Some(6.371e6),
            density_primary: Some(5000.0),
            density_secondary: Some(5000.0), // Same density
            ..Default::default()
        };

        let result = calculate_planetary(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let roche_radii = additional.get("roche_limit_primary_radii").unwrap();

        // When densities are equal, Roche limit = 2.44 radii
        assert!((*roche_radii - 2.44).abs() < 0.01);
    }

    #[test]
    fn test_planetary_missing_parameters() {
        let params = GeophysicsParams {
            mass_primary: Some(5.972e24),
            // Missing radius
            ..Default::default()
        };

        let result = calculate_planetary(&params);
        assert!(result.is_err());
    }

    // ===== INTEGRATION TESTS =====

    #[test]
    fn test_full_workflow_seismology() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::Seismology,
            parameters: GeophysicsParams {
                seismic_moment: Some(5e19),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert!(result.uncertainty.is_some());
        assert!(!result.interpretation.is_empty());
    }

    #[test]
    fn test_full_workflow_atmosphere() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::Atmosphere,
            parameters: GeophysicsParams {
                pressure: Some(101325.0),
                temperature: Some(288.15),
                altitude: Some(5000.0),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "Pa");
    }

    #[test]
    fn test_full_workflow_dating() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::RadiometricDating,
            parameters: GeophysicsParams {
                parent_isotope: Some(60.0),
                daughter_isotope: Some(40.0),
                isotope_system: Some("U238Pb206".to_string()),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "years");
    }

    #[test]
    fn test_full_workflow_planetary() {
        let input = GeophysicsInput {
            category: GeophysicsCategory::PlanetaryScience,
            parameters: GeophysicsParams {
                mass_primary: Some(5.972e24),
                radius_primary: Some(6.371e6),
                ..Default::default()
            },
        };

        let result = calculate_geophysics(input).unwrap();
        assert!(result.value > 0.0);
        assert_eq!(result.unit, "km/s");
    }
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
