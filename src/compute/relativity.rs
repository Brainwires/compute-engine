//! Relativity computation operations
//!
//! This module handles special and general relativity calculations including:
//! - Special relativity: Lorentz transformations, time dilation, length contraction
//! - General relativity: Schwarzschild metric, gravitational effects, black holes

use crate::engine::*;

/// Compute relativity operations
///
/// Handles both special and general relativity calculations:
///
/// **Special Relativity:**
/// - Lorentz transformation: Coordinate transformations between inertial frames
/// - Time dilation: Time passes slower for moving objects
/// - Length contraction: Objects contract along direction of motion
/// - Relativistic energy: Total energy including rest mass and kinetic energy
/// - Velocity addition: Relativistic velocity composition
///
/// **General Relativity:**
/// - Schwarzschild metric: Spacetime geometry around spherical masses
/// - Gravitational time dilation: Time passes slower in gravitational fields
/// - Orbital precession: Perihelion shift of orbits (e.g., Mercury)
/// - Gravitational lensing: Light bending around massive objects
/// - Black hole properties: Event horizon, Schwarzschild radius, etc.
pub fn compute_relativity(op: &RelativityOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use super::physics::relativity::*;

    let result_json = match op {
        RelativityOp::LorentzTransform => {
            let velocity = input
                .parameters
                .get("velocity")
                .and_then(|v| v.as_f64())
                .ok_or("velocity parameter required")?;
            let position: Vec<f64> = input
                .parameters
                .get("position")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("position [x, y, z] required")?;
            let time = input
                .parameters
                .get("time")
                .and_then(|v| v.as_f64())
                .ok_or("time parameter required")?;

            let result = lorentz_transform(LorentzTransformRequest {
                velocity,
                position,
                time,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::TimeDilation => {
            let proper_time = input
                .parameters
                .get("proper_time")
                .and_then(|v| v.as_f64())
                .ok_or("proper_time parameter required")?;
            let velocity = input
                .parameters
                .get("velocity")
                .and_then(|v| v.as_f64())
                .ok_or("velocity parameter required")?;

            let result = time_dilation(TimeDilationRequest {
                proper_time,
                velocity,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::LengthContraction => {
            let proper_length = input
                .parameters
                .get("proper_length")
                .and_then(|v| v.as_f64())
                .ok_or("proper_length parameter required")?;
            let velocity = input
                .parameters
                .get("velocity")
                .and_then(|v| v.as_f64())
                .ok_or("velocity parameter required")?;

            let result = length_contraction(LengthContractionRequest {
                proper_length,
                velocity,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::RelativisticEnergy => {
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;
            let velocity = input
                .parameters
                .get("velocity")
                .and_then(|v| v.as_f64())
                .ok_or("velocity parameter required")?;

            let result = relativistic_energy(RelativisticEnergyRequest { mass, velocity })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::VelocityAddition => {
            let velocity1 = input
                .parameters
                .get("velocity1")
                .and_then(|v| v.as_f64())
                .ok_or("velocity1 parameter required")?;
            let velocity2 = input
                .parameters
                .get("velocity2")
                .and_then(|v| v.as_f64())
                .ok_or("velocity2 parameter required")?;

            let result = velocity_addition(VelocityAdditionRequest {
                velocity1,
                velocity2,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::SchwarzschildMetric => {
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;
            let radius = input
                .parameters
                .get("radius")
                .and_then(|v| v.as_f64())
                .ok_or("radius parameter required")?;

            let result = schwarzschild_metric(SchwarzschildRequest { mass, radius })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::GravitationalTimeDilation => {
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;
            let radius = input
                .parameters
                .get("radius")
                .and_then(|v| v.as_f64())
                .ok_or("radius parameter required")?;
            let proper_time = input
                .parameters
                .get("proper_time")
                .and_then(|v| v.as_f64())
                .ok_or("proper_time parameter required")?;

            let result = gravitational_time_dilation(GravitationalTimeDilationRequest {
                mass,
                radius,
                proper_time,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::OrbitalPrecession => {
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;
            let semi_major_axis = input
                .parameters
                .get("semi_major_axis")
                .and_then(|v| v.as_f64())
                .ok_or("semi_major_axis parameter required")?;
            let eccentricity = input
                .parameters
                .get("eccentricity")
                .and_then(|v| v.as_f64())
                .ok_or("eccentricity parameter required")?;

            let result = orbital_precession(OrbitalPrecessionRequest {
                mass,
                semi_major_axis,
                eccentricity,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::GravitationalLensing => {
            let lens_mass = input
                .parameters
                .get("lens_mass")
                .and_then(|v| v.as_f64())
                .ok_or("lens_mass parameter required")?;
            let impact_parameter = input
                .parameters
                .get("impact_parameter")
                .and_then(|v| v.as_f64())
                .ok_or("impact_parameter parameter required")?;

            let result = gravitational_lensing(GravitationalLensingRequest {
                lens_mass,
                impact_parameter,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        RelativityOp::BlackHoleProperties => {
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;

            let result = black_hole_properties(BlackHoleRequest { mass })?;

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
