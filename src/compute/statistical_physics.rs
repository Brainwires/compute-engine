//! Statistical physics computation operations
//!
//! This module handles statistical mechanics and thermodynamics calculations including:
//! - Partition functions for various ensembles
//! - Distribution functions (Maxwell-Boltzmann, Fermi-Dirac, Bose-Einstein)
//! - Phase transitions and critical phenomena

use crate::engine::*;

/// Compute statistical physics operations
///
/// Handles statistical mechanics calculations:
///
/// **Partition Functions:**
/// - General partition function for canonical, grand canonical, and microcanonical ensembles
/// - Canonical partition function for fixed particle number
/// - Grand canonical partition function for variable particle number
///
/// **Distribution Functions:**
/// - Maxwell-Boltzmann: Classical particle velocity/energy distribution
/// - Fermi-Dirac: Fermion occupation probability (electrons, protons)
/// - Bose-Einstein: Boson occupation probability (photons, phonons)
///
/// **Thermodynamic Properties:**
/// - Chemical potential: Gibbs free energy per particle
/// - Fugacity coefficient: Deviation from ideal gas behavior
///
/// **Phase Transitions:**
/// - Phase transition analysis using various models
/// - Critical phenomena near phase transition points
pub fn compute_statistical_physics(
    op: &StatPhysicsOp,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use super::physics::statistical_physics::*;

    let result_json = match op {
        StatPhysicsOp::PartitionFunction => {
            let ensemble = input
                .parameters
                .get("ensemble")
                .and_then(|v| v.as_str())
                .ok_or("ensemble parameter required")?;
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;

            let energy_levels: Option<Vec<f64>> = input
                .parameters
                .get("energy_levels")
                .and_then(|v| serde_json::from_value(v.clone()).ok());
            let degeneracies: Option<Vec<usize>> = input
                .parameters
                .get("degeneracies")
                .and_then(|v| serde_json::from_value(v.clone()).ok());
            let volume = input.parameters.get("volume").and_then(|v| v.as_f64());
            let num_particles = input
                .parameters
                .get("num_particles")
                .and_then(|v| v.as_u64())
                .map(|v| v as usize);
            let chemical_potential = input
                .parameters
                .get("chemical_potential")
                .and_then(|v| v.as_f64());

            let result = partition_function(PartitionFunctionRequest {
                ensemble: ensemble.to_string(),
                temperature,
                volume,
                num_particles,
                energy_levels,
                degeneracies,
                chemical_potential,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::PartitionFunctionCanonical => {
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let energy_levels: Vec<f64> = input
                .parameters
                .get("energy_levels")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("energy_levels required")?;
            let degeneracies: Vec<usize> = input
                .parameters
                .get("degeneracies")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("degeneracies required")?;

            let result = canonical_partition(CanonicalPartitionRequest {
                temperature,
                energy_levels,
                degeneracies,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::PartitionFunctionGrandCanonical => {
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let volume = input
                .parameters
                .get("volume")
                .and_then(|v| v.as_f64())
                .ok_or("volume parameter required")?;
            let chemical_potential = input
                .parameters
                .get("chemical_potential")
                .and_then(|v| v.as_f64())
                .ok_or("chemical_potential parameter required")?;
            let particle_type = input
                .parameters
                .get("particle_type")
                .and_then(|v| v.as_str())
                .unwrap_or("classical");

            let result = grand_canonical_partition(GrandCanonicalRequest {
                temperature,
                volume,
                chemical_potential,
                particle_type: particle_type.to_string(),
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::MaxwellBoltzmann => {
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;
            let velocity = input.parameters.get("velocity").and_then(|v| v.as_f64());

            let result = maxwell_boltzmann(MaxwellBoltzmannRequest {
                temperature,
                mass,
                velocity,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::FermiDirac => {
            let energy = input
                .parameters
                .get("energy")
                .and_then(|v| v.as_f64())
                .ok_or("energy parameter required")?;
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let chemical_potential = input
                .parameters
                .get("chemical_potential")
                .and_then(|v| v.as_f64())
                .ok_or("chemical_potential parameter required")?;

            let result = fermi_dirac(FermiDiracRequest {
                energy,
                temperature,
                chemical_potential,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::BoseEinstein => {
            let energy = input
                .parameters
                .get("energy")
                .and_then(|v| v.as_f64())
                .ok_or("energy parameter required")?;
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let chemical_potential = input
                .parameters
                .get("chemical_potential")
                .and_then(|v| v.as_f64())
                .ok_or("chemical_potential parameter required")?;

            let result = bose_einstein(BoseEinsteinRequest {
                energy,
                temperature,
                chemical_potential,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::ChemicalPotential => {
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let pressure = input
                .parameters
                .get("pressure")
                .and_then(|v| v.as_f64())
                .ok_or("pressure parameter required")?;
            let particle_density = input
                .parameters
                .get("particle_density")
                .and_then(|v| v.as_f64())
                .ok_or("particle_density parameter required")?;
            let particle_type = input
                .parameters
                .get("particle_type")
                .and_then(|v| v.as_str())
                .unwrap_or("classical");
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;

            let result = chemical_potential(ChemicalPotentialRequest {
                temperature,
                pressure,
                particle_density,
                particle_type: particle_type.to_string(),
                mass,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::FugacityCoefficient => {
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let pressure = input
                .parameters
                .get("pressure")
                .and_then(|v| v.as_f64())
                .ok_or("pressure parameter required")?;
            let particle_density = input
                .parameters
                .get("particle_density")
                .and_then(|v| v.as_f64())
                .ok_or("particle_density parameter required")?;

            let result = fugacity_coefficient(FugacityRequest {
                temperature,
                pressure,
                particle_density,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::PhaseTransition => {
            let model = input
                .parameters
                .get("model")
                .and_then(|v| v.as_str())
                .ok_or("model parameter required")?;
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let critical_temperature = input
                .parameters
                .get("critical_temperature")
                .and_then(|v| v.as_f64())
                .ok_or("critical_temperature parameter required")?;

            let parameters = std::collections::HashMap::new();

            let result = phase_transition(PhaseTransitionRequest {
                model: model.to_string(),
                temperature,
                critical_temperature,
                parameters,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        StatPhysicsOp::CriticalPhenomena => {
            let temperature = input
                .parameters
                .get("temperature")
                .and_then(|v| v.as_f64())
                .ok_or("temperature parameter required")?;
            let critical_temperature = input
                .parameters
                .get("critical_temperature")
                .and_then(|v| v.as_f64())
                .ok_or("critical_temperature parameter required")?;
            let model = input
                .parameters
                .get("model")
                .and_then(|v| v.as_str())
                .unwrap_or("mean_field");

            let result = critical_phenomena(CriticalPhenomenaRequest {
                temperature,
                critical_temperature,
                model: model.to_string(),
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
