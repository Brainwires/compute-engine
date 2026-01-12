/**
 * Thermodynamics Module
 *
 * Implements fundamental thermodynamic operations:
 * - Heat transfer: Conduction, Convection, Radiation
 * - Thermal resistance networks
 * - Entropy calculations
 * - (Future: enthalpy, cycles, phase transitions)
 */
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermodynamicsInput {
    pub operation: ThermodynamicsOperation,
    pub parameters: ThermodynamicsParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThermodynamicsOperation {
    Conduction,
    Convection,
    Radiation,
    ThermalResistance,
    Entropy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermodynamicsParams {
    // Heat Transfer - Conduction
    pub thermal_conductivity: Option<f64>, // k (W/(m·K))
    pub area: Option<f64>,                 // A (m²)
    pub temperature_gradient: Option<f64>, // dT/dx (K/m)
    pub thickness: Option<f64>,            // L (m)
    pub temp_hot: Option<f64>,             // T1 (K or °C)
    pub temp_cold: Option<f64>,            // T2 (K or °C)

    // Heat Transfer - Convection
    pub heat_transfer_coefficient: Option<f64>, // h (W/(m²·K))
    pub surface_temp: Option<f64>,              // Ts (K or °C)
    pub fluid_temp: Option<f64>,                // T∞ (K or °C)

    // Heat Transfer - Radiation
    pub emissivity: Option<f64>,     // ε (0-1)
    pub surface_temp_1: Option<f64>, // T1 (K)
    pub surface_temp_2: Option<f64>, // T2 (K)

    // Thermal resistance
    pub resistances: Option<Vec<f64>>, // R values (K/W)
    pub configuration: Option<String>, // "series" or "parallel"

    // Entropy calculations
    pub heat_transfer: Option<f64>,   // Q (J)
    pub temperature: Option<f64>,     // T (K)
    pub mass: Option<f64>,            // m (kg)
    pub specific_heat: Option<f64>,   // c (J/(kg·K))
    pub initial_temp: Option<f64>,    // T1 (K)
    pub final_temp: Option<f64>,      // T2 (K)
    pub num_microstates: Option<f64>, // Ω (for Boltzmann entropy)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermodynamicsResult {
    pub value: f64,
    pub unit: String,
    pub formula_used: String,
    pub interpretation: String,
    pub additional_info: Option<serde_json::Value>,
}

const STEFAN_BOLTZMANN: f64 = 5.67e-8; // W/(m²·K⁴)
const BOLTZMANN_CONSTANT: f64 = 1.380649e-23; // J/K

pub fn calculate_thermodynamics(
    input: ThermodynamicsInput,
) -> Result<ThermodynamicsResult, String> {
    match input.operation {
        ThermodynamicsOperation::Conduction => calculate_conduction(&input.parameters),
        ThermodynamicsOperation::Convection => calculate_convection(&input.parameters),
        ThermodynamicsOperation::Radiation => calculate_radiation(&input.parameters),
        ThermodynamicsOperation::ThermalResistance => calculate_resistance(&input.parameters),
        ThermodynamicsOperation::Entropy => calculate_entropy(&input.parameters),
    }
}

/// Fourier's law: q = -k·A·(dT/dx)
fn calculate_conduction(params: &ThermodynamicsParams) -> Result<ThermodynamicsResult, String> {
    let k = params
        .thermal_conductivity
        .ok_or("Thermal conductivity required")?;
    let area = params.area.ok_or("Area required")?;

    let heat_rate = if let Some(grad) = params.temperature_gradient {
        // Direct gradient provided
        -k * area * grad
    } else {
        // Calculate from temperatures and thickness
        let t_hot = params.temp_hot.ok_or("Hot temperature required")?;
        let t_cold = params.temp_cold.ok_or("Cold temperature required")?;
        let thickness = params.thickness.ok_or("Thickness required")?;

        if thickness <= 0.0 {
            return Err("Thickness must be positive".to_string());
        }

        k * area * (t_hot - t_cold) / thickness
    };

    let thermal_r = params.thickness.map(|l| l / (k * area));

    Ok(ThermodynamicsResult {
        value: heat_rate.abs(),
        unit: "W".to_string(),
        formula_used: "Fourier's law: q = k·A·ΔT/L".to_string(),
        interpretation: format!("Conduction heat transfer: {:.2} W", heat_rate.abs()),
        additional_info: Some(serde_json::json!({
            "thermal_resistance": thermal_r
        })),
    })
}

/// Newton's law of cooling: q = h·A·(Ts - T∞)
fn calculate_convection(params: &ThermodynamicsParams) -> Result<ThermodynamicsResult, String> {
    let h = params
        .heat_transfer_coefficient
        .ok_or("Heat transfer coefficient required")?;
    let area = params.area.ok_or("Area required")?;
    let t_surface = params.surface_temp.ok_or("Surface temperature required")?;
    let t_fluid = params.fluid_temp.ok_or("Fluid temperature required")?;

    let heat_rate = h * area * (t_surface - t_fluid);
    let thermal_r = 1.0 / (h * area);

    let flow_type = if h < 10.0 {
        "Natural convection (low h)"
    } else if h < 100.0 {
        "Forced convection (moderate h)"
    } else {
        "Phase change or high velocity (high h)"
    };

    Ok(ThermodynamicsResult {
        value: heat_rate.abs(),
        unit: "W".to_string(),
        formula_used: "Newton's law: q = h·A·(Ts - T∞)".to_string(),
        interpretation: format!("{}: {:.2} W", flow_type, heat_rate.abs()),
        additional_info: Some(serde_json::json!({
            "thermal_resistance": thermal_r
        })),
    })
}

/// Stefan-Boltzmann law: q = σ·ε·A·(T₁⁴ - T₂⁴)
fn calculate_radiation(params: &ThermodynamicsParams) -> Result<ThermodynamicsResult, String> {
    let emissivity = params.emissivity.ok_or("Emissivity required")?;
    let area = params.area.ok_or("Area required")?;
    let t1 = params
        .surface_temp_1
        .ok_or("Surface 1 temperature (K) required")?;
    let t2 = params.surface_temp_2.unwrap_or(0.0); // Default to radiation to space

    if emissivity < 0.0 || emissivity > 1.0 {
        return Err("Emissivity must be between 0 and 1".to_string());
    }

    if t1 <= 0.0 {
        return Err("Temperature must be in Kelvin (> 0)".to_string());
    }

    let heat_rate = STEFAN_BOLTZMANN * emissivity * area * (t1.powi(4) - t2.powi(4));

    let surface_type = if emissivity > 0.9 {
        "Black body (ε ≈ 1)"
    } else if emissivity > 0.5 {
        "Gray body"
    } else {
        "Reflective surface (low ε)"
    };

    Ok(ThermodynamicsResult {
        value: heat_rate.abs(),
        unit: "W".to_string(),
        formula_used: "Stefan-Boltzmann: q = σ·ε·A·(T₁⁴ - T₂⁴)".to_string(),
        interpretation: format!("{}: {:.2} W", surface_type, heat_rate.abs()),
        additional_info: Some(serde_json::json!({
            "surface_type": surface_type
        })),
    })
}

/// Thermal resistance networks
fn calculate_resistance(params: &ThermodynamicsParams) -> Result<ThermodynamicsResult, String> {
    let resistances = params.resistances.as_ref().ok_or("Resistances required")?;
    let config = params
        .configuration
        .as_ref()
        .map(|s| s.as_str())
        .unwrap_or("series");

    if resistances.is_empty() {
        return Err("At least one resistance value required".to_string());
    }

    let total_r = match config {
        "series" => {
            // Series: R_total = R1 + R2 + R3 + ...
            resistances.iter().sum()
        }
        "parallel" => {
            // Parallel: 1/R_total = 1/R1 + 1/R2 + 1/R3 + ...
            let sum_reciprocals: f64 = resistances.iter().map(|r| 1.0 / r).sum();
            1.0 / sum_reciprocals
        }
        _ => return Err("Configuration must be 'series' or 'parallel'".to_string()),
    };

    // Calculate heat rate if temperatures provided
    let heat_rate = if let (Some(t_hot), Some(t_cold)) = (params.temp_hot, params.temp_cold) {
        (t_hot - t_cold) / total_r
    } else {
        0.0
    };

    Ok(ThermodynamicsResult {
        value: heat_rate,
        unit: "W".to_string(),
        formula_used: format!("Thermal resistance ({} network)", config),
        interpretation: format!("Total thermal resistance: {:.4} K/W", total_r),
        additional_info: Some(serde_json::json!({
            "thermal_resistance": total_r,
            "configuration": config
        })),
    })
}

/// Entropy calculations
/// Supports: Clausius definition (ΔS = Q/T), Boltzmann entropy (S = k·ln(Ω)), and thermal entropy change
fn calculate_entropy(params: &ThermodynamicsParams) -> Result<ThermodynamicsResult, String> {
    // Method 1: Clausius definition ΔS = Q/T
    if let (Some(q), Some(t)) = (params.heat_transfer, params.temperature) {
        if t <= 0.0 {
            return Err("Temperature must be positive (use Kelvin)".to_string());
        }

        let delta_s = q / t;

        let interpretation = if delta_s > 0.0 {
            format!(
                "Entropy increases by {:.4} J/K (heat added to system)",
                delta_s
            )
        } else if delta_s < 0.0 {
            format!(
                "Entropy decreases by {:.4} J/K (heat removed from system)",
                delta_s.abs()
            )
        } else {
            "No entropy change (reversible process)".to_string()
        };

        return Ok(ThermodynamicsResult {
            value: delta_s,
            unit: "J/K".to_string(),
            formula_used: "Clausius: ΔS = Q/T".to_string(),
            interpretation,
            additional_info: Some(serde_json::json!({
                "heat_transfer": q,
                "temperature": t,
                "reversible": delta_s.abs() < 1e-10
            })),
        });
    }

    // Method 2: Boltzmann entropy S = k·ln(Ω)
    if let Some(omega) = params.num_microstates {
        if omega <= 0.0 {
            return Err("Number of microstates must be positive".to_string());
        }

        let s = BOLTZMANN_CONSTANT * omega.ln();

        return Ok(ThermodynamicsResult {
            value: s,
            unit: "J/K".to_string(),
            formula_used: "Boltzmann: S = k·ln(Ω)".to_string(),
            interpretation: format!("Statistical entropy for {} microstates", omega),
            additional_info: Some(serde_json::json!({
                "microstates": omega,
                "boltzmann_constant": BOLTZMANN_CONSTANT
            })),
        });
    }

    // Method 3: Thermal entropy change ΔS = m·c·ln(T₂/T₁)
    if let (Some(m), Some(c), Some(t1), Some(t2)) = (
        params.mass,
        params.specific_heat,
        params.initial_temp,
        params.final_temp,
    ) {
        if t1 <= 0.0 || t2 <= 0.0 {
            return Err("Temperatures must be positive (use Kelvin)".to_string());
        }

        let delta_s = m * c * (t2 / t1).ln();

        let interpretation = if t2 > t1 {
            format!(
                "Heating: entropy increases by {:.4} J/K ({:.1}K → {:.1}K)",
                delta_s, t1, t2
            )
        } else if t2 < t1 {
            format!(
                "Cooling: entropy changes by {:.4} J/K ({:.1}K → {:.1}K)",
                delta_s, t1, t2
            )
        } else {
            "No temperature change, no entropy change".to_string()
        };

        return Ok(ThermodynamicsResult {
            value: delta_s,
            unit: "J/K".to_string(),
            formula_used: "Thermal: ΔS = m·c·ln(T₂/T₁)".to_string(),
            interpretation,
            additional_info: Some(serde_json::json!({
                "mass": m,
                "specific_heat": c,
                "initial_temp": t1,
                "final_temp": t2,
                "temperature_ratio": t2 / t1
            })),
        });
    }

    Err(
        "Insufficient parameters for entropy calculation. Provide either:\n\
         1. heat_transfer + temperature (Clausius)\n\
         2. num_microstates (Boltzmann)\n\
         3. mass + specific_heat + initial_temp + final_temp (Thermal)"
            .to_string(),
    )
}


impl Default for ThermodynamicsParams {
    fn default() -> Self {
        Self {
            thermal_conductivity: None,
            area: None,
            temperature_gradient: None,
            thickness: None,
            temp_hot: None,
            temp_cold: None,
            heat_transfer_coefficient: None,
            surface_temp: None,
            fluid_temp: None,
            emissivity: None,
            surface_temp_1: None,
            surface_temp_2: None,
            resistances: None,
            configuration: None,
            heat_transfer: None,
            temperature: None,
            mass: None,
            specific_heat: None,
            initial_temp: None,
            final_temp: None,
            num_microstates: None,
        }
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/thermodynamics_tests.rs"]
mod tests;
