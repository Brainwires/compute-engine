//! Nuclear physics computation operations
//!
//! This module handles nuclear physics calculations including:
//! - Radioactive decay and decay chains
//! - Binding energy and mass defect
//! - Nuclear reactions (fission and fusion)
//! - Nuclear stability and energy release

use crate::engine::*;

/// Compute nuclear physics operations
///
/// Handles nuclear physics calculations:
///
/// **Radioactive Decay:**
/// - Radioactive decay: Exponential decay law N(t) = N₀ exp(-λt)
/// - Decay chains: Parent → Daughter → Granddaughter sequences
/// - Half-life: Relationship between decay constant and half-life
///
/// **Nuclear Binding:**
/// - Binding energy: Energy required to disassemble nucleus
/// - Mass defect: Difference between nuclear mass and sum of constituents
/// - Semi-empirical mass formula (SEMF): Predict nuclear masses
///
/// **Nuclear Reactions:**
/// - Fission energy: Energy released in nuclear fission
/// - Fusion energy: Energy released in nuclear fusion
/// - Q-value: Net energy released in nuclear reaction
/// - Projectile energy effects on reaction dynamics
///
/// **Applications:**
/// - Nuclear power calculations
/// - Stellar nucleosynthesis
/// - Medical isotope production
/// - Nuclear dating (carbon-14, uranium-lead)
pub fn compute_nuclear_physics(op: &NuclearOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::physics::nuclear_physics::*;

    let result_json = match op {
        NuclearOp::RadioactiveDecay => {
            let initial_quantity = input
                .parameters
                .get("initial_quantity")
                .and_then(|v| v.as_f64())
                .ok_or("initial_quantity parameter required")?;
            let decay_constant = input
                .parameters
                .get("decay_constant")
                .and_then(|v| v.as_f64())
                .ok_or("decay_constant parameter required")?;
            let time = input
                .parameters
                .get("time")
                .and_then(|v| v.as_f64())
                .ok_or("time parameter required")?;

            let result = radioactive_decay(RadioactiveDecayRequest {
                initial_quantity,
                decay_constant,
                time,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::DecayChain => {
            let parent_initial = input
                .parameters
                .get("parent_initial")
                .and_then(|v| v.as_f64())
                .ok_or("parent_initial parameter required")?;
            let parent_decay_constant = input
                .parameters
                .get("parent_decay_constant")
                .and_then(|v| v.as_f64())
                .ok_or("parent_decay_constant parameter required")?;
            let daughter_decay_constant = input
                .parameters
                .get("daughter_decay_constant")
                .and_then(|v| v.as_f64())
                .ok_or("daughter_decay_constant parameter required")?;
            let time = input
                .parameters
                .get("time")
                .and_then(|v| v.as_f64())
                .ok_or("time parameter required")?;

            let result = decay_chain(DecayChainRequest {
                parent_initial,
                parent_decay_constant,
                daughter_decay_constant,
                time,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::HalfLife => {
            let decay_constant_val = input
                .parameters
                .get("decay_constant")
                .and_then(|v| v.as_f64());
            let half_life_val = input.parameters.get("half_life").and_then(|v| v.as_f64());
            let mean_lifetime_val = input
                .parameters
                .get("mean_lifetime")
                .and_then(|v| v.as_f64());

            let result = half_life(HalfLifeRequest {
                decay_constant: decay_constant_val,
                half_life: half_life_val,
                mean_lifetime: mean_lifetime_val,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::BindingEnergy => {
            let atomic_number = input
                .parameters
                .get("atomic_number")
                .and_then(|v| v.as_u64())
                .ok_or("atomic_number parameter required")? as u32;
            let mass_number = input
                .parameters
                .get("mass_number")
                .and_then(|v| v.as_u64())
                .ok_or("mass_number parameter required")? as u32;
            let atomic_mass = input.parameters.get("atomic_mass").and_then(|v| v.as_f64());

            let result = binding_energy(BindingEnergyRequest {
                atomic_number,
                mass_number,
                atomic_mass,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::MassDefect => {
            let protons = input
                .parameters
                .get("protons")
                .and_then(|v| v.as_u64())
                .ok_or("protons parameter required")? as u32;
            let neutrons = input
                .parameters
                .get("neutrons")
                .and_then(|v| v.as_u64())
                .ok_or("neutrons parameter required")? as u32;
            let nuclear_mass = input
                .parameters
                .get("nuclear_mass")
                .and_then(|v| v.as_f64())
                .ok_or("nuclear_mass parameter required")?;

            let result = mass_defect(MassDefectRequest {
                protons,
                neutrons,
                nuclear_mass,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::FissionEnergy => {
            let parent_mass = input
                .parameters
                .get("parent_mass")
                .and_then(|v| v.as_f64())
                .ok_or("parent_mass parameter required")?;
            let fragment1_mass = input
                .parameters
                .get("fragment1_mass")
                .and_then(|v| v.as_f64())
                .ok_or("fragment1_mass parameter required")?;
            let fragment2_mass = input
                .parameters
                .get("fragment2_mass")
                .and_then(|v| v.as_f64())
                .ok_or("fragment2_mass parameter required")?;
            let neutrons_released = input
                .parameters
                .get("neutrons_released")
                .and_then(|v| v.as_u64())
                .ok_or("neutrons_released parameter required")?
                as u32;

            let result = fission_energy(FissionEnergyRequest {
                parent_mass,
                fragment1_mass,
                fragment2_mass,
                neutrons_released,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::FusionEnergy => {
            let reactant1_mass = input
                .parameters
                .get("reactant1_mass")
                .and_then(|v| v.as_f64())
                .ok_or("reactant1_mass parameter required")?;
            let reactant2_mass = input
                .parameters
                .get("reactant2_mass")
                .and_then(|v| v.as_f64())
                .ok_or("reactant2_mass parameter required")?;
            let product1_mass = input
                .parameters
                .get("product1_mass")
                .and_then(|v| v.as_f64())
                .ok_or("product1_mass parameter required")?;
            let product2_mass = input
                .parameters
                .get("product2_mass")
                .and_then(|v| v.as_f64());

            let result = fusion_energy(FusionEnergyRequest {
                reactant1_mass,
                reactant2_mass,
                product1_mass,
                product2_mass,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        NuclearOp::NuclearReaction => {
            let reactants: Vec<f64> = input
                .parameters
                .get("reactants")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("reactants parameter required as array of masses in amu")?;
            let products: Vec<f64> = input
                .parameters
                .get("products")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("products parameter required as array of masses in amu")?;
            let projectile_energy_mev = input
                .parameters
                .get("projectile_energy_mev")
                .and_then(|v| v.as_f64());

            let result = nuclear_reaction(NuclearReactionRequest {
                reactants,
                products,
                projectile_energy_mev,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }
    };

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}
