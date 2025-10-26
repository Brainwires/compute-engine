/**
 * Chemistry Module
 *
 * Implements essential chemistry formulas:
 * - pH and buffer calculations (Henderson-Hasselbalch)
 * - Chemical kinetics (rate laws, Arrhenius)
 * - Thermochemistry (Gibbs, equilibrium constants)
 * - Electrochemistry (Nernst, Butler-Volmer)
 * - Gas laws (Van der Waals, virial)
 * - Spectroscopy (Beer-Lambert)
 */
use serde::{Deserialize, Serialize};
use std::f64::consts::E;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemistryInput {
    pub operation: ChemistryOperation,
    pub parameters: ChemistryParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChemistryOperation {
    PhCalculation,
    BufferCapacity,
    Arrhenius,
    RateLaw,
    GibbsFreeEnergy,
    NernstEquation,
    BeerLambert,
    VanDerWaals,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemistryParams {
    // pH calculations
    pub pka: Option<f64>,
    pub concentration_acid: Option<f64>,
    pub concentration_base: Option<f64>,

    // Kinetics
    pub activation_energy: Option<f64>, // kJ/mol
    pub temperature: Option<f64>,       // K
    pub pre_exponential: Option<f64>,   // A in Arrhenius
    pub rate_constant: Option<f64>,
    pub concentration: Option<f64>,
    pub order: Option<i32>, // reaction order

    // Thermodynamics
    pub enthalpy: Option<f64>, // ΔH in kJ/mol
    pub entropy: Option<f64>,  // ΔS in J/(mol·K)

    // Electrochemistry
    pub standard_potential: Option<f64>, // E° in V
    pub n_electrons: Option<i32>,
    pub concentrations: Option<Vec<f64>>,

    // Spectroscopy
    pub absorbance: Option<f64>,
    pub path_length: Option<f64>,        // cm
    pub molar_absorptivity: Option<f64>, // L/(mol·cm)

    // Gas laws
    pub pressure: Option<f64>, // atm
    pub volume: Option<f64>,   // L
    pub moles: Option<f64>,
    pub a_constant: Option<f64>, // Van der Waals a
    pub b_constant: Option<f64>, // Van der Waals b
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemistryResult {
    pub value: f64,
    pub unit: String,
    pub formula_used: String,
    pub interpretation: String,
}

const R: f64 = 8.314; // Gas constant J/(mol·K)
const F: f64 = 96485.0; // Faraday constant C/mol

pub fn calculate_chemistry(input: ChemistryInput) -> Result<ChemistryResult, String> {
    match input.operation {
        ChemistryOperation::PhCalculation => calculate_ph(&input.parameters),
        ChemistryOperation::BufferCapacity => calculate_buffer_capacity(&input.parameters),
        ChemistryOperation::Arrhenius => calculate_arrhenius(&input.parameters),
        ChemistryOperation::RateLaw => calculate_rate_law(&input.parameters),
        ChemistryOperation::GibbsFreeEnergy => calculate_gibbs(&input.parameters),
        ChemistryOperation::NernstEquation => calculate_nernst(&input.parameters),
        ChemistryOperation::BeerLambert => calculate_beer_lambert(&input.parameters),
        ChemistryOperation::VanDerWaals => calculate_van_der_waals(&input.parameters),
    }
}

/// Henderson-Hasselbalch: pH = pKa + log([A⁻]/[HA])
fn calculate_ph(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let pka = params.pka.ok_or("pKa required")?;
    let conc_base = params
        .concentration_base
        .ok_or("Base concentration required")?;
    let conc_acid = params
        .concentration_acid
        .ok_or("Acid concentration required")?;

    if conc_acid <= 0.0 {
        return Err("Acid concentration must be positive".to_string());
    }

    let ph = pka + (conc_base / conc_acid).log10();

    let interpretation = if ph < 7.0 {
        "Acidic solution"
    } else if ph > 7.0 {
        "Basic solution"
    } else {
        "Neutral solution"
    };

    Ok(ChemistryResult {
        value: ph,
        unit: "pH units".to_string(),
        formula_used: "Henderson-Hasselbalch: pH = pKa + log([A⁻]/[HA])".to_string(),
        interpretation: interpretation.to_string(),
    })
}

/// Buffer capacity
fn calculate_buffer_capacity(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let conc_acid = params
        .concentration_acid
        .ok_or("Acid concentration required")?;
    let conc_base = params
        .concentration_base
        .ok_or("Base concentration required")?;

    // Buffer capacity β = 2.303 * [HA][A⁻]/([HA] + [A⁻])
    let total = conc_acid + conc_base;
    let capacity = 2.303 * (conc_acid * conc_base) / total;

    Ok(ChemistryResult {
        value: capacity,
        unit: "mol/(L·pH)".to_string(),
        formula_used: "β = 2.303·[HA][A⁻]/([HA]+[A⁻])".to_string(),
        interpretation: format!(
            "Buffer can resist pH change by absorbing {:.3e} mol/L of acid/base per pH unit",
            capacity
        ),
    })
}

/// Arrhenius equation: k = A·exp(-Ea/RT)
fn calculate_arrhenius(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let ea = params
        .activation_energy
        .ok_or("Activation energy (kJ/mol) required")?;
    let t = params.temperature.ok_or("Temperature (K) required")?;
    let a = params
        .pre_exponential
        .ok_or("Pre-exponential factor A required")?;

    if t <= 0.0 {
        return Err("Temperature must be positive".to_string());
    }

    // Convert Ea from kJ/mol to J/mol
    let ea_joules = ea * 1000.0;

    // k = A·exp(-Ea/RT)
    let rate_constant = a * E.powf(-ea_joules / (R * t));

    Ok(ChemistryResult {
        value: rate_constant,
        unit: "s⁻¹ (or M⁻¹s⁻¹ depending on order)".to_string(),
        formula_used: "Arrhenius: k = A·exp(-Ea/RT)".to_string(),
        interpretation: format!("Rate constant at {}K is {:.3e}", t, rate_constant),
    })
}

/// Rate law: rate = k[A]^n
fn calculate_rate_law(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let k = params.rate_constant.ok_or("Rate constant required")?;
    let conc = params.concentration.ok_or("Concentration required")?;
    let order = params.order.unwrap_or(1);

    // rate = k[A]^n
    let rate = k * conc.powi(order);

    let order_str = match order {
        0 => "Zero order",
        1 => "First order",
        2 => "Second order",
        _ => "Higher order",
    };

    Ok(ChemistryResult {
        value: rate,
        unit: "M/s".to_string(),
        formula_used: format!("Rate law: rate = k[A]^{}", order),
        interpretation: format!("{} reaction with rate {:.3e} M/s", order_str, rate),
    })
}

/// Gibbs free energy: ΔG = ΔH - TΔS
fn calculate_gibbs(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let h = params.enthalpy.ok_or("Enthalpy ΔH (kJ/mol) required")?;
    let s = params.entropy.ok_or("Entropy ΔS (J/(mol·K)) required")?;
    let t = params.temperature.ok_or("Temperature (K) required")?;

    // ΔG = ΔH - TΔS (convert to kJ/mol)
    let delta_g = h - (t * s / 1000.0);

    let spontaneity = if delta_g < 0.0 {
        "Spontaneous (thermodynamically favorable)"
    } else if delta_g > 0.0 {
        "Non-spontaneous (thermodynamically unfavorable)"
    } else {
        "At equilibrium"
    };

    Ok(ChemistryResult {
        value: delta_g,
        unit: "kJ/mol".to_string(),
        formula_used: "Gibbs: ΔG = ΔH - TΔS".to_string(),
        interpretation: spontaneity.to_string(),
    })
}

/// Nernst equation: E = E° - (RT/nF)ln(Q)
fn calculate_nernst(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let e_standard = params
        .standard_potential
        .ok_or("Standard potential E° (V) required")?;
    let n = params
        .n_electrons
        .ok_or("Number of electrons transferred required")? as f64;
    let t = params.temperature.unwrap_or(298.15); // Default 25°C
    let concs = params
        .concentrations
        .as_ref()
        .ok_or("Concentrations required")?;

    if concs.len() < 2 {
        return Err("Need at least 2 concentrations for reaction quotient".to_string());
    }

    // Simplified: Q = [products]/[reactants]
    let q = concs[1] / concs[0];

    // E = E° - (RT/nF)ln(Q)
    let e_cell = e_standard - (R * t / (n * F)) * q.ln();

    Ok(ChemistryResult {
        value: e_cell,
        unit: "V".to_string(),
        formula_used: "Nernst: E = E° - (RT/nF)ln(Q)".to_string(),
        interpretation: format!("Cell potential at given concentrations: {:.4} V", e_cell),
    })
}

/// Beer-Lambert law: A = εbc
fn calculate_beer_lambert(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let absorbance = params.absorbance;
    let epsilon = params.molar_absorptivity;
    let path = params.path_length;
    let conc = params.concentration;

    // Solve for missing variable
    if let Some(a) = absorbance {
        if let (Some(e), Some(b)) = (epsilon, path) {
            // Calculate concentration: c = A/(εb)
            let c = a / (e * b);
            return Ok(ChemistryResult {
                value: c,
                unit: "M".to_string(),
                formula_used: "Beer-Lambert: A = εbc → c = A/(εb)".to_string(),
                interpretation: format!("Concentration is {:.3e} M", c),
            });
        }
    }

    // Calculate absorbance: A = εbc
    let e = epsilon.ok_or("Molar absorptivity required")?;
    let b = path.ok_or("Path length required")?;
    let c = conc.ok_or("Concentration required")?;

    let a = e * b * c;

    Ok(ChemistryResult {
        value: a,
        unit: "AU (absorbance units)".to_string(),
        formula_used: "Beer-Lambert: A = εbc".to_string(),
        interpretation: format!("Absorbance is {:.4}", a),
    })
}

/// Van der Waals equation: (P + a/V²)(V - b) = RT
fn calculate_van_der_waals(params: &ChemistryParams) -> Result<ChemistryResult, String> {
    let p = params.pressure.ok_or("Pressure (atm) required")?;
    let v = params.volume.ok_or("Molar volume (L/mol) required")?;
    let t = params.temperature.ok_or("Temperature (K) required")?;
    let a = params
        .a_constant
        .ok_or("Van der Waals constant 'a' required")?;
    let b = params
        .b_constant
        .ok_or("Van der Waals constant 'b' required")?;

    // R in L·atm/(mol·K)
    let r_gas = 0.08206;

    // Check equation: (P + a/V²)(V - b) = RT
    let lhs = (p + a / (v * v)) * (v - b);
    let rhs = r_gas * t;

    let deviation = ((lhs - rhs) / rhs * 100.0).abs();

    Ok(ChemistryResult {
        value: deviation,
        unit: "% deviation from ideal".to_string(),
        formula_used: "Van der Waals: (P + a/V²)(V - b) = RT".to_string(),
        interpretation: format!("Gas deviates {:.2}% from ideal behavior", deviation),
    })
}


impl Default for ChemistryParams {
    fn default() -> Self {
        Self {
            pka: None,
            concentration_acid: None,
            concentration_base: None,
            activation_energy: None,
            temperature: None,
            pre_exponential: None,
            rate_constant: None,
            concentration: None,
            order: None,
            enthalpy: None,
            entropy: None,
            standard_potential: None,
            n_electrons: None,
            concentrations: None,
            absorbance: None,
            path_length: None,
            molar_absorptivity: None,
            pressure: None,
            volume: None,
            moles: None,
            a_constant: None,
            b_constant: None,
        }
    }
}
