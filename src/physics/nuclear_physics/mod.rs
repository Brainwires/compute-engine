//! Nuclear Physics Module
//!
//! Implements nuclear physics operations:
//! - Radioactive Decay and Decay Chains
//! - Half-Life Calculations
//! - Binding Energy and Mass Defect
//! - Nuclear Fission and Fusion Energy
//! - Nuclear Reactions

use serde::{Deserialize, Serialize};

// Physical constants
const C: f64 = 299792458.0; // Speed of light (m/s)
const AMU: f64 = 1.66053906660e-27; // Atomic mass unit (kg)
const MEV_PER_AMU: f64 = 931.494; // MeV per amu
const PROTON_MASS: f64 = 1.007276466812; // amu
const NEUTRON_MASS: f64 = 1.00866491588; // amu
const ELECTRON_MASS: f64 = 0.00054857990946; // amu

// ============================================================================
// RADIOACTIVE DECAY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct RadioactiveDecayRequest {
    pub initial_quantity: f64, // N₀ (number of atoms or activity)
    pub decay_constant: f64,   // λ (1/s)
    pub time: f64,             // t (seconds)
}

#[derive(Debug, Serialize)]
pub struct RadioactiveDecayResult {
    pub remaining_quantity: f64, // N(t) = N₀e^(-λt)
    pub decayed_quantity: f64,   // N₀ - N(t)
    pub activity: f64,           // A(t) = λN(t)
    pub decay_rate: f64,         // dN/dt = -λN(t)
    pub fraction_remaining: f64, // N(t)/N₀
}

/// Calculate radioactive decay over time
pub fn radioactive_decay(
    request: RadioactiveDecayRequest,
) -> Result<RadioactiveDecayResult, String> {
    if request.decay_constant <= 0.0 {
        return Err("Decay constant must be positive".to_string());
    }
    if request.time < 0.0 {
        return Err("Time cannot be negative".to_string());
    }

    // N(t) = N₀ * exp(-λt)
    let remaining = request.initial_quantity * (-request.decay_constant * request.time).exp();
    let decayed = request.initial_quantity - remaining;
    let activity = request.decay_constant * remaining;
    let decay_rate = -request.decay_constant * remaining;
    let fraction = remaining / request.initial_quantity;

    Ok(RadioactiveDecayResult {
        remaining_quantity: remaining,
        decayed_quantity: decayed,
        activity,
        decay_rate,
        fraction_remaining: fraction,
    })
}

// ============================================================================
// DECAY CHAIN
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct DecayChainRequest {
    pub parent_initial: f64,          // N₁(0)
    pub parent_decay_constant: f64,   // λ₁
    pub daughter_decay_constant: f64, // λ₂
    pub time: f64,                    // t
}

#[derive(Debug, Serialize)]
pub struct DecayChainResult {
    pub parent_quantity: f64,        // N₁(t)
    pub daughter_quantity: f64,      // N₂(t)
    pub parent_activity: f64,        // A₁(t)
    pub daughter_activity: f64,      // A₂(t)
    pub secular_equilibrium: bool,   // λ₁ << λ₂
    pub transient_equilibrium: bool, // λ₁ < λ₂
}

/// Model radioactive decay chains (parent → daughter)
pub fn decay_chain(request: DecayChainRequest) -> Result<DecayChainResult, String> {
    if request.parent_decay_constant <= 0.0 || request.daughter_decay_constant <= 0.0 {
        return Err("Decay constants must be positive".to_string());
    }

    let lambda1 = request.parent_decay_constant;
    let lambda2 = request.daughter_decay_constant;
    let n10 = request.parent_initial;
    let t = request.time;

    // Parent: N₁(t) = N₁(0) * exp(-λ₁t)
    let n1 = n10 * (-lambda1 * t).exp();

    // Daughter: N₂(t) = (λ₁N₁(0))/(λ₂ - λ₁) * [exp(-λ₁t) - exp(-λ₂t)]
    let n2 = if (lambda2 - lambda1).abs() > 1e-10 {
        (lambda1 * n10 / (lambda2 - lambda1)) * ((-lambda1 * t).exp() - (-lambda2 * t).exp())
    } else {
        // When λ₁ ≈ λ₂, use limiting case
        lambda1 * n10 * t * (-lambda1 * t).exp()
    };

    let a1 = lambda1 * n1;
    let a2 = lambda2 * n2;

    // Equilibrium conditions
    let secular = lambda1 < lambda2 * 0.01; // Parent half-life >> daughter half-life
    let transient = lambda1 < lambda2 && !secular;

    Ok(DecayChainResult {
        parent_quantity: n1,
        daughter_quantity: n2.max(0.0),
        parent_activity: a1,
        daughter_activity: a2,
        secular_equilibrium: secular,
        transient_equilibrium: transient,
    })
}

// ============================================================================
// HALF-LIFE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct HalfLifeRequest {
    pub decay_constant: Option<f64>, // λ (1/s)
    pub half_life: Option<f64>,      // t₁/₂ (seconds)
    pub mean_lifetime: Option<f64>,  // τ (seconds)
}

#[derive(Debug, Serialize)]
pub struct HalfLifeResult {
    pub half_life: f64,         // t₁/₂ = ln(2)/λ
    pub decay_constant: f64,    // λ = ln(2)/t₁/₂
    pub mean_lifetime: f64,     // τ = 1/λ
    pub decay_probability: f64, // Probability of decay per unit time
}

/// Calculate half-life, decay constant, and related quantities
pub fn half_life(request: HalfLifeRequest) -> Result<HalfLifeResult, String> {
    let lambda = if let Some(l) = request.decay_constant {
        if l <= 0.0 {
            return Err("Decay constant must be positive".to_string());
        }
        l
    } else if let Some(t_half) = request.half_life {
        if t_half <= 0.0 {
            return Err("Half-life must be positive".to_string());
        }
        2.0_f64.ln() / t_half
    } else if let Some(tau) = request.mean_lifetime {
        if tau <= 0.0 {
            return Err("Mean lifetime must be positive".to_string());
        }
        1.0 / tau
    } else {
        return Err("Must provide decay_constant, half_life, or mean_lifetime".to_string());
    };

    let t_half = 2.0_f64.ln() / lambda;
    let tau = 1.0 / lambda;
    let prob = 1.0 - (-lambda).exp(); // Probability per unit time

    Ok(HalfLifeResult {
        half_life: t_half,
        decay_constant: lambda,
        mean_lifetime: tau,
        decay_probability: prob,
    })
}

// ============================================================================
// BINDING ENERGY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct BindingEnergyRequest {
    pub atomic_number: u32,       // Z (number of protons)
    pub mass_number: u32,         // A (total nucleons)
    pub atomic_mass: Option<f64>, // Measured atomic mass in amu
}

#[derive(Debug, Serialize)]
pub struct BindingEnergyResult {
    pub binding_energy_mev: f64,         // Total BE
    pub binding_energy_per_nucleon: f64, // BE/A
    pub mass_defect_amu: f64,            // Δm in amu
    pub mass_defect_kg: f64,             // Δm in kg
    pub separation_energy: f64,          // Approximate
}

/// Calculate nuclear binding energy
pub fn binding_energy(request: BindingEnergyRequest) -> Result<BindingEnergyResult, String> {
    let z = request.atomic_number as f64;
    let a = request.mass_number as f64;
    let n = a - z; // Number of neutrons

    if z <= 0.0 || a <= 0.0 || n < 0.0 {
        return Err("Invalid atomic or mass number".to_string());
    }

    // Mass defect: Δm = Zm_p + Nm_n - m_nucleus
    let expected_mass = z * PROTON_MASS + n * NEUTRON_MASS;
    let actual_mass = request.atomic_mass.unwrap_or_else(|| {
        // Use semi-empirical mass formula (SEMF) if not provided
        // This is a simplified approximation
        let a_v = 15.75; // Volume term
        let a_s = 17.8; // Surface term
        let a_c = 0.711; // Coulomb term
        let a_a = 23.7; // Asymmetry term
        let a_p = 11.18; // Pairing term

        let binding = a_v * a
            - a_s * a.powf(2.0 / 3.0)
            - a_c * z * (z - 1.0) / a.powf(1.0 / 3.0)
            - a_a * (a - 2.0 * z).powi(2) / a
            + if a as u32 % 2 == 0 {
                a_p / a.sqrt()
            } else {
                0.0
            };

        expected_mass - binding / MEV_PER_AMU
    });

    let mass_defect_amu = expected_mass - actual_mass;
    let mass_defect_kg = mass_defect_amu * AMU;

    // Binding energy: E = Δm * c²
    let be_mev = mass_defect_amu * MEV_PER_AMU;
    let be_per_nucleon = be_mev / a;

    // Approximate separation energy (last nucleon)
    let sep_energy = be_per_nucleon * 0.8; // Simplified

    Ok(BindingEnergyResult {
        binding_energy_mev: be_mev,
        binding_energy_per_nucleon: be_per_nucleon,
        mass_defect_amu,
        mass_defect_kg,
        separation_energy: sep_energy,
    })
}

// ============================================================================
// MASS DEFECT
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct MassDefectRequest {
    pub protons: u32,
    pub neutrons: u32,
    pub nuclear_mass: f64, // amu
}

#[derive(Debug, Serialize)]
pub struct MassDefectResult {
    pub mass_defect_amu: f64,
    pub mass_defect_kg: f64,
    pub expected_mass: f64,
    pub actual_mass: f64,
    pub binding_energy_mev: f64,
    pub binding_energy_joules: f64,
}

/// Calculate mass defect
pub fn mass_defect(request: MassDefectRequest) -> Result<MassDefectResult, String> {
    let expected = request.protons as f64 * PROTON_MASS + request.neutrons as f64 * NEUTRON_MASS;

    let defect_amu = expected - request.nuclear_mass;
    let defect_kg = defect_amu * AMU;

    // E = mc²
    let be_mev = defect_amu * MEV_PER_AMU;
    let be_joules = defect_kg * C * C;

    Ok(MassDefectResult {
        mass_defect_amu: defect_amu,
        mass_defect_kg: defect_kg,
        expected_mass: expected,
        actual_mass: request.nuclear_mass,
        binding_energy_mev: be_mev,
        binding_energy_joules: be_joules,
    })
}

// ============================================================================
// FISSION ENERGY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct FissionEnergyRequest {
    pub parent_mass: f64,       // amu
    pub fragment1_mass: f64,    // amu
    pub fragment2_mass: f64,    // amu
    pub neutrons_released: u32, // Number of neutrons
}

#[derive(Debug, Serialize)]
pub struct FissionEnergyResult {
    pub q_value_mev: f64,               // Total energy released
    pub kinetic_energy_mev: f64,        // KE of fragments
    pub energy_per_neutron_mev: f64,    // Average energy per neutron
    pub mass_to_energy_efficiency: f64, // Fraction of mass converted to energy
}

/// Calculate energy released in nuclear fission
pub fn fission_energy(request: FissionEnergyRequest) -> Result<FissionEnergyResult, String> {
    let neutron_mass_total = request.neutrons_released as f64 * NEUTRON_MASS;

    // Q-value: Q = (m_parent - m_products) * c²
    let mass_diff = request.parent_mass
        - (request.fragment1_mass + request.fragment2_mass + neutron_mass_total);

    let q_value_mev = mass_diff * MEV_PER_AMU;

    if q_value_mev < 0.0 {
        return Err("Endothermic reaction (Q < 0)".to_string());
    }

    // Most energy goes to kinetic energy of fragments (~167 MeV for U-235)
    let kinetic_energy = q_value_mev * 0.85; // ~85% as kinetic energy

    let energy_per_neutron = if request.neutrons_released > 0 {
        q_value_mev / request.neutrons_released as f64
    } else {
        0.0
    };

    let efficiency = mass_diff / request.parent_mass;

    Ok(FissionEnergyResult {
        q_value_mev,
        kinetic_energy_mev: kinetic_energy,
        energy_per_neutron_mev: energy_per_neutron,
        mass_to_energy_efficiency: efficiency,
    })
}

// ============================================================================
// FUSION ENERGY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct FusionEnergyRequest {
    pub reactant1_mass: f64,        // amu
    pub reactant2_mass: f64,        // amu
    pub product1_mass: f64,         // amu
    pub product2_mass: Option<f64>, // amu (if two products)
}

#[derive(Debug, Serialize)]
pub struct FusionEnergyResult {
    pub q_value_mev: f64,
    pub energy_per_nucleon_mev: f64,
    pub cross_section_barns: f64,      // Approximate
    pub ignition_temperature_kev: f64, // Approximate
}

/// Calculate energy released in nuclear fusion
pub fn fusion_energy(request: FusionEnergyRequest) -> Result<FusionEnergyResult, String> {
    let total_reactant = request.reactant1_mass + request.reactant2_mass;
    let total_product = request.product1_mass + request.product2_mass.unwrap_or(0.0);

    // Q-value: Q = (m_reactants - m_products) * c²
    let mass_diff = total_reactant - total_product;
    let q_value_mev = mass_diff * MEV_PER_AMU;

    if q_value_mev < 0.0 {
        return Err("Endothermic reaction (Q < 0)".to_string());
    }

    // Approximate number of nucleons involved
    let nucleons = (total_reactant).round();
    let energy_per_nucleon = q_value_mev / nucleons;

    // Approximate cross-section (varies greatly with energy)
    let cross_section = 0.1; // barns (very rough estimate for D-T fusion)

    // Approximate ignition temperature (D-T fusion ~15 keV)
    let ignition_temp = 15.0; // keV

    Ok(FusionEnergyResult {
        q_value_mev,
        energy_per_nucleon_mev: energy_per_nucleon,
        cross_section_barns: cross_section,
        ignition_temperature_kev: ignition_temp,
    })
}

// ============================================================================
// NUCLEAR REACTION
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct NuclearReactionRequest {
    pub reactants: Vec<f64>,                // Masses in amu
    pub products: Vec<f64>,                 // Masses in amu
    pub projectile_energy_mev: Option<f64>, // Kinetic energy of projectile
}

#[derive(Debug, Serialize)]
pub struct NuclearReactionResult {
    pub q_value_mev: f64,
    pub threshold_energy_mev: f64, // Minimum projectile energy for reaction
    pub reaction_type: String,     // "exothermic" or "endothermic"
    pub available_energy_mev: f64, // Energy available to products
    pub cross_section_barns: f64,  // Approximate
}

/// Calculate Q-value and properties of nuclear reactions
pub fn nuclear_reaction(request: NuclearReactionRequest) -> Result<NuclearReactionResult, String> {
    if request.reactants.is_empty() || request.products.is_empty() {
        return Err("Must have at least one reactant and one product".to_string());
    }

    let total_reactant: f64 = request.reactants.iter().sum();
    let total_product: f64 = request.products.iter().sum();

    // Q-value = (m_reactants - m_products) * c²
    let mass_diff = total_reactant - total_product;
    let q_value_mev = mass_diff * MEV_PER_AMU;

    let reaction_type = if q_value_mev > 0.0 {
        "exothermic"
    } else {
        "endothermic"
    };

    // Threshold energy for endothermic reactions
    // E_threshold ≈ -Q * (1 + m_products/(2*m_target))
    let threshold = if q_value_mev < 0.0 {
        let m_target = request.reactants.get(0).copied().unwrap_or(1.0);
        -q_value_mev * (1.0 + total_product / (2.0 * m_target * MEV_PER_AMU))
    } else {
        0.0
    };

    // Available energy = Q + projectile KE
    let projectile_ke = request.projectile_energy_mev.unwrap_or(0.0);
    let available = q_value_mev + projectile_ke;

    // Approximate cross-section (highly simplified)
    let cross_section = if projectile_ke > threshold {
        0.5 // barns (very rough)
    } else {
        0.0
    };

    Ok(NuclearReactionResult {
        q_value_mev,
        threshold_energy_mev: threshold,
        reaction_type: reaction_type.to_string(),
        available_energy_mev: available,
        cross_section_barns: cross_section,
    })
}

