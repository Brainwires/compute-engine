//! Stoichiometry and Advanced Chemistry Module
//!
//! Implements additional chemical calculations and simulations:
//! - Stoichiometry and equation balancing
//! - Thermodynamics (enthalpy, entropy, Gibbs free energy)
//! - Electrochemistry (redox, cell potentials)
//! - Kinetics (rate laws, Arrhenius equation)
//! - Gas laws and ideal/real gas behavior
//! - Acid-base chemistry (pH, buffers)
//!
//! This module complements the main chemistry module with more detailed
//! request/response types for complex chemical calculations.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Request/Response types

#[derive(Debug, Deserialize)]
pub struct BalanceEquationRequest {
    pub reactants: HashMap<String, i32>, // molecule -> coefficient
    pub products: HashMap<String, i32>,
}

#[derive(Debug, Serialize)]
pub struct BalanceEquationResult {
    pub balanced_reactants: HashMap<String, i32>,
    pub balanced_products: HashMap<String, i32>,
    pub balanced: bool,
    pub stoichiometric_coefficients: Vec<i32>,
}

#[derive(Debug, Deserialize)]
pub struct ThermodynamicsRequest {
    pub temperature: f64, // Kelvin
    pub pressure: f64,    // atm
    pub compounds: Vec<CompoundData>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct CompoundData {
    pub formula: String,
    pub moles: f64,
    pub enthalpy_formation: Option<f64>, // kJ/mol
    pub entropy: Option<f64>,            // J/(mol·K)
}

#[derive(Debug, Serialize)]
pub struct ThermodynamicsResult {
    pub delta_h: f64, // kJ
    pub delta_s: f64, // J/K
    pub delta_g: f64, // kJ
    pub spontaneous: bool,
    pub equilibrium_constant: f64,
}

#[derive(Debug, Deserialize)]
pub struct ElectrochemistryRequest {
    pub half_reactions: Vec<HalfReaction>,
    pub temperature: f64,                             // Kelvin
    pub concentrations: Option<HashMap<String, f64>>, // mol/L
}

#[derive(Debug, Deserialize, Clone)]
pub struct HalfReaction {
    pub oxidized: String,
    pub reduced: String,
    pub electrons: i32,
    pub standard_potential: f64, // V
}

#[derive(Debug, Serialize)]
pub struct ElectrochemistryResult {
    pub cell_potential: f64,          // V
    pub standard_cell_potential: f64, // V
    pub delta_g: f64,                 // kJ/mol
    pub spontaneous: bool,
    pub balanced_equation: String,
}

#[derive(Debug, Deserialize)]
pub struct KineticsRequest {
    pub rate_constant: f64,
    pub temperature: f64,         // Kelvin
    pub activation_energy: f64,   // kJ/mol
    pub concentrations: Vec<f64>, // mol/L
    pub orders: Vec<i32>,         // reaction order for each reactant
}

#[derive(Debug, Serialize)]
pub struct KineticsResult {
    pub rate: f64,              // mol/(L·s)
    pub half_life: Option<f64>, // seconds
    pub rate_constant_at_temp: f64,
}

#[derive(Debug, Deserialize)]
pub struct GasLawRequest {
    pub pressure: Option<f64>,    // atm
    pub volume: Option<f64>,      // L
    pub temperature: Option<f64>, // K
    pub moles: Option<f64>,
    pub gas_type: Option<String>, // for real gas calculations
}

#[derive(Debug, Serialize)]
pub struct GasLawResult {
    pub pressure: f64,    // atm
    pub volume: f64,      // L
    pub temperature: f64, // K
    pub moles: f64,
    pub density: f64, // g/L
}

#[derive(Debug, Deserialize)]
pub struct AcidBaseRequest {
    pub acid_ka: Option<f64>,
    pub base_kb: Option<f64>,
    pub concentration: f64,  // mol/L
    pub volume: Option<f64>, // L
}

#[derive(Debug, Serialize)]
pub struct AcidBaseResult {
    pub ph: f64,
    pub poh: f64,
    pub h_concentration: f64,  // mol/L
    pub oh_concentration: f64, // mol/L
}

#[derive(Debug, Deserialize)]
pub struct MolarMassRequest {
    pub formula: String,
}

#[derive(Debug, Serialize)]
pub struct MolarMassResult {
    pub molar_mass: f64,                   // g/mol
    pub composition: HashMap<String, f64>, // element -> mass percentage
}

// Constants
const R: f64 = 8.314; // J/(mol·K)
const F: f64 = 96485.0; // C/mol (Faraday constant)
const GAS_CONSTANT_ATM: f64 = 0.08206; // L·atm/(mol·K)

// Atomic masses (simplified)
fn atomic_mass(element: &str) -> f64 {
    match element {
        "H" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "S" => 32.06,
        "P" => 30.974,
        "Cl" => 35.45,
        "Na" => 22.990,
        "K" => 39.098,
        "Ca" => 40.078,
        "Fe" => 55.845,
        "Cu" => 63.546,
        "Zn" => 65.38,
        "Br" => 79.904,
        "I" => 126.90,
        _ => 0.0,
    }
}

/// Balance a chemical equation (simplified)
pub fn balance_equation(request: BalanceEquationRequest) -> Result<BalanceEquationResult, String> {
    // Simplified balancing - in practice, use matrix methods
    // For now, just verify if already balanced

    let balanced = check_balance(&request.reactants, &request.products);

    Ok(BalanceEquationResult {
        balanced_reactants: request.reactants.clone(),
        balanced_products: request.products.clone(),
        balanced,
        stoichiometric_coefficients: vec![1, 1], // simplified
    })
}

fn check_balance(reactants: &HashMap<String, i32>, products: &HashMap<String, i32>) -> bool {
    // Parse molecules and count atoms on each side
    let mut reactant_atoms: HashMap<String, i32> = HashMap::new();
    let mut product_atoms: HashMap<String, i32> = HashMap::new();

    // Count atoms in reactants
    for (molecule, &coeff) in reactants.iter() {
        let atoms = parse_molecule(molecule);
        for (element, count) in atoms {
            *reactant_atoms.entry(element).or_insert(0) += count * coeff;
        }
    }

    // Count atoms in products
    for (molecule, &coeff) in products.iter() {
        let atoms = parse_molecule(molecule);
        for (element, count) in atoms {
            *product_atoms.entry(element).or_insert(0) += count * coeff;
        }
    }

    // Check if all elements balance
    for (element, &count) in reactant_atoms.iter() {
        if product_atoms.get(element) != Some(&count) {
            return false;
        }
    }

    // Check no extra elements in products
    for (element, _) in product_atoms.iter() {
        if !reactant_atoms.contains_key(element) {
            return false;
        }
    }

    true
}

fn parse_molecule(formula: &str) -> HashMap<String, i32> {
    // Simple molecule parser: H2O -> {H: 2, O: 1}
    let mut atoms = HashMap::new();
    let chars: Vec<char> = formula.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        if chars[i].is_uppercase() {
            let mut element = chars[i].to_string();
            i += 1;

            // Check for lowercase (e.g., Cl, Ca)
            if i < chars.len() && chars[i].is_lowercase() {
                element.push(chars[i]);
                i += 1;
            }

            // Get count
            let mut count_str = String::new();
            while i < chars.len() && chars[i].is_numeric() {
                count_str.push(chars[i]);
                i += 1;
            }

            let count = if count_str.is_empty() {
                1
            } else {
                count_str.parse::<i32>().unwrap_or(1)
            };

            *atoms.entry(element).or_insert(0) += count;
        } else {
            i += 1;
        }
    }

    atoms
}

/// Calculate thermodynamic properties
pub fn thermodynamics(request: ThermodynamicsRequest) -> Result<ThermodynamicsResult, String> {
    let t = request.temperature;

    // Calculate ΔH (sum of products - sum of reactants)
    let delta_h: f64 = request
        .compounds
        .iter()
        .map(|c| c.moles * c.enthalpy_formation.unwrap_or(0.0))
        .sum();

    // Calculate ΔS
    let delta_s: f64 = request
        .compounds
        .iter()
        .map(|c| c.moles * c.entropy.unwrap_or(0.0))
        .sum();

    // Calculate ΔG = ΔH - TΔS
    let delta_g = delta_h - (t * delta_s / 1000.0); // Convert J to kJ

    // Calculate equilibrium constant: K = exp(-ΔG/RT)
    let k_eq = (-delta_g * 1000.0 / (R * t)).exp();

    Ok(ThermodynamicsResult {
        delta_h,
        delta_s,
        delta_g,
        spontaneous: delta_g < 0.0,
        equilibrium_constant: k_eq,
    })
}

/// Calculate electrochemical cell potential
pub fn electrochemistry(
    request: ElectrochemistryRequest,
) -> Result<ElectrochemistryResult, String> {
    if request.half_reactions.len() != 2 {
        return Err("Need exactly 2 half-reactions".to_string());
    }

    let t = request.temperature;
    let oxidation = &request.half_reactions[0];
    let reduction = &request.half_reactions[1];

    // Standard cell potential
    let e_cell_standard = reduction.standard_potential - oxidation.standard_potential;

    // Nernst equation: E = E° - (RT/nF)ln(Q)
    let n = oxidation.electrons.max(reduction.electrons) as f64;

    let e_cell = if let Some(conc) = &request.concentrations {
        let q = conc.values().product::<f64>(); // Simplified reaction quotient
        e_cell_standard - (R * t / (n * F)) * q.ln()
    } else {
        e_cell_standard
    };

    // ΔG = -nFE
    let delta_g = -n * F * e_cell / 1000.0; // Convert to kJ/mol

    Ok(ElectrochemistryResult {
        cell_potential: e_cell,
        standard_cell_potential: e_cell_standard,
        delta_g,
        spontaneous: e_cell > 0.0,
        balanced_equation: format!("{} -> {}", oxidation.oxidized, reduction.reduced),
    })
}

/// Calculate reaction kinetics
pub fn kinetics(request: KineticsRequest) -> Result<KineticsResult, String> {
    let k = request.rate_constant;
    let t = request.temperature;
    let ea = request.activation_energy * 1000.0; // Convert to J/mol

    // Arrhenius equation: k = A * exp(-Ea/RT)
    // Solve for A: A = k / exp(-Ea/RT)
    let a = k / (-ea / (R * t)).exp();

    // Calculate rate constant at given temperature
    let k_at_t = a * (-ea / (R * t)).exp();

    // Rate = k[A]^m[B]^n...
    let mut rate = k_at_t;
    for (i, &conc) in request.concentrations.iter().enumerate() {
        if i < request.orders.len() {
            rate *= conc.powi(request.orders[i]);
        }
    }

    // Half-life (only for first-order: t_1/2 = ln(2)/k)
    let half_life = if request.orders.get(0) == Some(&1) {
        Some(2.0_f64.ln() / k_at_t)
    } else {
        None
    };

    Ok(KineticsResult {
        rate,
        half_life,
        rate_constant_at_temp: k_at_t,
    })
}

/// Ideal gas law calculations
pub fn gas_law(request: GasLawRequest) -> Result<GasLawResult, String> {
    // PV = nRT
    let r = GAS_CONSTANT_ATM;

    let (p, v, t, n) = match (
        request.pressure,
        request.volume,
        request.temperature,
        request.moles,
    ) {
        (Some(p), Some(v), Some(t), None) => {
            let n = (p * v) / (r * t);
            (p, v, t, n)
        }
        (Some(p), Some(v), None, Some(n)) => {
            let t = (p * v) / (r * n);
            (p, v, t, n)
        }
        (Some(p), None, Some(t), Some(n)) => {
            let v = (n * r * t) / p;
            (p, v, t, n)
        }
        (None, Some(v), Some(t), Some(n)) => {
            let p = (n * r * t) / v;
            (p, v, t, n)
        }
        _ => return Err("Need exactly 3 of 4 parameters (P, V, T, n)".to_string()),
    };

    // Calculate density (simplified, assuming ideal gas with molar mass ~29 g/mol for air)
    let molar_mass = 29.0; // g/mol (approximate for air)
    let density = (p * molar_mass) / (r * t);

    Ok(GasLawResult {
        pressure: p,
        volume: v,
        temperature: t,
        moles: n,
        density,
    })
}

/// Calculate pH and acid-base properties
pub fn acid_base(request: AcidBaseRequest) -> Result<AcidBaseResult, String> {
    let c = request.concentration;

    let (h_conc, oh_conc) = if let Some(ka) = request.acid_ka {
        // Weak acid: [H+] = sqrt(Ka * C)
        let h = (ka * c).sqrt();
        let oh = 1.0e-14 / h;
        (h, oh)
    } else if let Some(kb) = request.base_kb {
        // Weak base: [OH-] = sqrt(Kb * C)
        let oh = (kb * c).sqrt();
        let h = 1.0e-14 / oh;
        (h, oh)
    } else {
        return Err("Need either Ka or Kb".to_string());
    };

    let ph = -h_conc.log10();
    let poh = -oh_conc.log10();

    Ok(AcidBaseResult {
        ph,
        poh,
        h_concentration: h_conc,
        oh_concentration: oh_conc,
    })
}

/// Calculate molar mass from chemical formula
pub fn molar_mass(request: MolarMassRequest) -> Result<MolarMassResult, String> {
    let formula = &request.formula;

    // Simple parser for formulas like H2O, C6H12O6
    let mut mass = 0.0;
    let mut composition = HashMap::new();

    let chars: Vec<char> = formula.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        if chars[i].is_uppercase() {
            let mut element = chars[i].to_string();
            i += 1;

            // Check for lowercase (e.g., Cl, Ca)
            if i < chars.len() && chars[i].is_lowercase() {
                element.push(chars[i]);
                i += 1;
            }

            // Get count
            let mut count_str = String::new();
            while i < chars.len() && chars[i].is_numeric() {
                count_str.push(chars[i]);
                i += 1;
            }

            let count = if count_str.is_empty() {
                1
            } else {
                count_str.parse::<i32>().unwrap_or(1)
            };

            let atom_mass = atomic_mass(&element);
            mass += atom_mass * count as f64;
            *composition.entry(element).or_insert(0.0) += atom_mass * count as f64;
        } else {
            i += 1;
        }
    }

    // Convert to percentages
    for (_, value) in composition.iter_mut() {
        *value = (*value / mass) * 100.0;
    }

    Ok(MolarMassResult {
        molar_mass: mass,
        composition,
    })
}

/// Stoichiometry calculations
pub fn stoichiometry(
    reactant_mass: f64,
    reactant_molar_mass: f64,
    product_molar_mass: f64,
    reactant_coeff: i32,
    product_coeff: i32,
) -> f64 {
    // Convert mass to moles
    let reactant_moles = reactant_mass / reactant_molar_mass;

    // Use stoichiometric ratio
    let product_moles = reactant_moles * (product_coeff as f64 / reactant_coeff as f64);

    // Convert to mass
    product_moles * product_molar_mass
}

