/**
 * Biology Module
 *
 * Implements essential biology formulas:
 * - Enzyme kinetics (Michaelis-Menten, Lineweaver-Burk)
 * - Pharmacokinetics (compartment models, clearance)
 * - Population genetics (Hardy-Weinberg)
 * - Membrane potentials (Goldman-Hodgkin-Katz)
 * - Allometric scaling
 */
use serde::{Deserialize, Serialize};
use std::f64::consts::E;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiologyInput {
    pub operation: BiologyOperation,
    pub parameters: BiologyParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BiologyOperation {
    MichaelisMenten,
    LineweaverBurk,
    Pharmacokinetics,
    HardyWeinberg,
    GoldmanEquation,
    AllometricScaling,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiologyParams {
    // Enzyme kinetics
    pub vmax: Option<f64>,
    pub km: Option<f64>,
    pub substrate_concentration: Option<f64>,
    pub inhibitor_concentration: Option<f64>,
    pub ki: Option<f64>, // Inhibition constant

    // Pharmacokinetics
    pub dose: Option<f64>,
    pub clearance: Option<f64>,
    pub volume_distribution: Option<f64>,
    pub bioavailability: Option<f64>,
    pub elimination_rate: Option<f64>,
    pub time: Option<f64>,

    // Population genetics
    pub allele_frequency_p: Option<f64>,

    // Membrane potentials
    pub ion_concentrations_inside: Option<Vec<f64>>,
    pub ion_concentrations_outside: Option<Vec<f64>>,
    pub permeabilities: Option<Vec<f64>>,
    pub temperature: Option<f64>,
    pub valence: Option<Vec<i32>>,

    // Allometric scaling
    pub body_mass: Option<f64>,       // kg
    pub scaling_type: Option<String>, // "metabolic", "surface_area", "lifespan"
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiologyResult {
    pub value: f64,
    pub unit: String,
    pub formula_used: String,
    pub interpretation: String,
    pub additional_data: Option<std::collections::HashMap<String, f64>>,
}

const R: f64 = 8.314; // Gas constant J/(mol·K)
const F: f64 = 96485.0; // Faraday constant C/mol

pub fn calculate_biology(input: BiologyInput) -> Result<BiologyResult, String> {
    match input.operation {
        BiologyOperation::MichaelisMenten => calculate_michaelis_menten(&input.parameters),
        BiologyOperation::LineweaverBurk => calculate_lineweaver_burk(&input.parameters),
        BiologyOperation::Pharmacokinetics => calculate_pharmacokinetics(&input.parameters),
        BiologyOperation::HardyWeinberg => calculate_hardy_weinberg(&input.parameters),
        BiologyOperation::GoldmanEquation => calculate_goldman(&input.parameters),
        BiologyOperation::AllometricScaling => calculate_allometric(&input.parameters),
    }
}

/// Michaelis-Menten: v = (Vmax·[S])/(Km + [S])
fn calculate_michaelis_menten(params: &BiologyParams) -> Result<BiologyResult, String> {
    let vmax = params.vmax.ok_or("Vmax required")?;
    let km = params.km.ok_or("Km required")?;
    let s = params
        .substrate_concentration
        .ok_or("Substrate concentration required")?;

    if km <= 0.0 || vmax <= 0.0 {
        return Err("Vmax and Km must be positive".to_string());
    }

    // v = (Vmax·[S])/(Km + [S])
    let velocity = (vmax * s) / (km + s);

    // Calculate efficiency
    let efficiency = (velocity / vmax) * 100.0;

    let mut additional = std::collections::HashMap::new();
    additional.insert("efficiency_percent".to_string(), efficiency);
    additional.insert("half_vmax_at_km".to_string(), km);

    let interpretation = if s < km {
        "Below Km - enzyme not saturated, first-order kinetics"
    } else if s > 10.0 * km {
        "Well above Km - enzyme saturated, zero-order kinetics"
    } else {
        "Near Km - mixed-order kinetics"
    };

    Ok(BiologyResult {
        value: velocity,
        unit: "μmol/(min·mg enzyme)".to_string(),
        formula_used: "Michaelis-Menten: v = (Vmax·[S])/(Km + [S])".to_string(),
        interpretation: interpretation.to_string(),
        additional_data: Some(additional),
    })
}

/// Lineweaver-Burk plot: 1/v = (Km/Vmax)·(1/[S]) + 1/Vmax
fn calculate_lineweaver_burk(params: &BiologyParams) -> Result<BiologyResult, String> {
    let vmax = params.vmax.ok_or("Vmax required")?;
    let km = params.km.ok_or("Km required")?;
    let s = params
        .substrate_concentration
        .ok_or("Substrate concentration required")?;

    if s <= 0.0 {
        return Err("Substrate concentration must be positive".to_string());
    }

    // Calculate velocity first
    let v = (vmax * s) / (km + s);

    // Lineweaver-Burk coordinates
    let x = 1.0 / s; // 1/[S]
    let y = 1.0 / v; // 1/v

    let mut additional = std::collections::HashMap::new();
    additional.insert("x_intercept_1_over_s".to_string(), x);
    additional.insert("y_intercept_1_over_v".to_string(), y);
    additional.insert("slope_km_over_vmax".to_string(), km / vmax);
    additional.insert("y_intercept_1_over_vmax".to_string(), 1.0 / vmax);

    Ok(BiologyResult {
        value: y,
        unit: "1/v".to_string(),
        formula_used: "Lineweaver-Burk: 1/v = (Km/Vmax)·(1/[S]) + 1/Vmax".to_string(),
        interpretation: format!(
            "Double reciprocal plot point: (1/[S]={:.3}, 1/v={:.3})",
            x, y
        ),
        additional_data: Some(additional),
    })
}

/// One-compartment pharmacokinetics: C(t) = (F·Dose/V)·exp(-k·t)
fn calculate_pharmacokinetics(params: &BiologyParams) -> Result<BiologyResult, String> {
    let dose = params.dose.ok_or("Dose (mg) required")?;
    let v = params
        .volume_distribution
        .ok_or("Volume of distribution (L) required")?;
    let f = params.bioavailability.unwrap_or(1.0); // Default 100%
    let k = params
        .elimination_rate
        .ok_or("Elimination rate constant (1/h) required")?;
    let t = params.time.ok_or("Time (h) required")?;

    // C(t) = (F·Dose/V)·exp(-k·t)
    let concentration = (f * dose / v) * E.powf(-k * t);

    // Calculate half-life
    let t_half = 0.693 / k;

    let mut additional = std::collections::HashMap::new();
    additional.insert("half_life_hours".to_string(), t_half);
    additional.insert("clearance_l_per_h".to_string(), k * v);
    additional.insert("auc_mg_h_per_l".to_string(), f * dose / (k * v));

    Ok(BiologyResult {
        value: concentration,
        unit: "mg/L".to_string(),
        formula_used: "One-compartment: C(t) = (F·Dose/V)·exp(-k·t)".to_string(),
        interpretation: format!(
            "Plasma concentration at t={:.1}h is {:.3} mg/L (t½={:.2}h)",
            t, concentration, t_half
        ),
        additional_data: Some(additional),
    })
}

/// Hardy-Weinberg equilibrium: p² + 2pq + q² = 1
fn calculate_hardy_weinberg(params: &BiologyParams) -> Result<BiologyResult, String> {
    let p = params
        .allele_frequency_p
        .ok_or("Allele frequency p required")?;

    if p < 0.0 || p > 1.0 {
        return Err("Allele frequency must be between 0 and 1".to_string());
    }

    let q = 1.0 - p;

    // Genotype frequencies
    let aa = p * p; // Homozygous dominant
    let aa_het = 2.0 * p * q; // Heterozygous
    let aa_rec = q * q; // Homozygous recessive

    let mut additional = std::collections::HashMap::new();
    additional.insert("freq_AA_homozygous_dominant".to_string(), aa);
    additional.insert("freq_Aa_heterozygous".to_string(), aa_het);
    additional.insert("freq_aa_homozygous_recessive".to_string(), aa_rec);
    additional.insert("allele_freq_q".to_string(), q);

    Ok(BiologyResult {
        value: aa,
        unit: "frequency".to_string(),
        formula_used: "Hardy-Weinberg: p² + 2pq + q² = 1".to_string(),
        interpretation: format!(
            "AA={:.3}, Aa={:.3}, aa={:.3} (p={:.3}, q={:.3})",
            aa, aa_het, aa_rec, p, q
        ),
        additional_data: Some(additional),
    })
}

/// Goldman-Hodgkin-Katz equation for membrane potential
fn calculate_goldman(params: &BiologyParams) -> Result<BiologyResult, String> {
    let inside = params
        .ion_concentrations_inside
        .as_ref()
        .ok_or("Inside ion concentrations required")?;
    let outside = params
        .ion_concentrations_outside
        .as_ref()
        .ok_or("Outside ion concentrations required")?;
    let perms = params
        .permeabilities
        .as_ref()
        .ok_or("Permeabilities required")?;
    let t = params.temperature.unwrap_or(310.0); // Default 37°C

    if inside.len() != outside.len() || inside.len() != perms.len() {
        return Err("Ion concentrations and permeabilities must have same length".to_string());
    }

    // Simplified for Na+, K+, Cl- (most common)
    // Em = (RT/F)·ln((PK[K+]out + PNa[Na+]out + PCl[Cl-]in)/(PK[K+]in + PNa[Na+]in + PCl[Cl-]out))

    let mut numerator = 0.0;
    let mut denominator = 0.0;

    // For cations (positive ions)
    for i in 0..inside.len().min(2) {
        numerator += perms[i] * outside[i];
        denominator += perms[i] * inside[i];
    }

    // For anions (if present) - reverse inside/outside
    if inside.len() > 2 {
        numerator += perms[2] * inside[2];
        denominator += perms[2] * outside[2];
    }

    if denominator == 0.0 {
        return Err("Invalid ion concentrations".to_string());
    }

    // Em = (RT/F)·ln(numerator/denominator)
    // Convert to mV
    let potential = (R * t / F) * (numerator / denominator).ln() * 1000.0;

    Ok(BiologyResult {
        value: potential,
        unit: "mV".to_string(),
        formula_used: "Goldman-Hodgkin-Katz equation".to_string(),
        interpretation: format!("Membrane potential: {:.1} mV", potential),
        additional_data: None,
    })
}

/// Allometric scaling: Y = a·M^b (Kleiber's law for metabolism: b ≈ 0.75)
fn calculate_allometric(params: &BiologyParams) -> Result<BiologyResult, String> {
    let mass = params.body_mass.ok_or("Body mass (kg) required")?;
    let scaling_type = params
        .scaling_type
        .as_ref()
        .ok_or("Scaling type required (metabolic, surface_area, lifespan)")?;

    let (a, b, unit, description) = match scaling_type.as_str() {
        "metabolic" => {
            // Kleiber's law: BMR = 70·M^0.75 (kcal/day for mammals)
            (70.0, 0.75, "kcal/day", "Basal metabolic rate")
        }
        "surface_area" => {
            // Body surface area: BSA = 0.1·M^0.67 (m² for mammals)
            (0.1, 0.67, "m²", "Body surface area")
        }
        "lifespan" => {
            // Lifespan scales as M^0.25 (approximate)
            (10.0, 0.25, "years", "Maximum lifespan estimate")
        }
        _ => return Err("Unknown scaling type".to_string()),
    };

    let value = a * mass.powf(b);

    Ok(BiologyResult {
        value,
        unit: unit.to_string(),
        formula_used: format!("Allometric: Y = {:.1}·M^{:.2}", a, b),
        interpretation: format!(
            "{}: {:.2} {} for {:.1} kg organism",
            description, value, unit, mass
        ),
        additional_data: None,
    })
}


impl Default for BiologyParams {
    fn default() -> Self {
        Self {
            vmax: None,
            km: None,
            substrate_concentration: None,
            inhibitor_concentration: None,
            ki: None,
            dose: None,
            clearance: None,
            volume_distribution: None,
            bioavailability: None,
            elimination_rate: None,
            time: None,
            allele_frequency_p: None,
            ion_concentrations_inside: None,
            ion_concentrations_outside: None,
            permeabilities: None,
            temperature: None,
            valence: None,
            body_mass: None,
            scaling_type: None,
        }
    }
}

// Test module
#[cfg(test)]
#[path = "../../../tests/unit/compute/biology/biology_tests.rs"]
mod tests;
