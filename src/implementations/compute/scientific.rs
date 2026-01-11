//! Scientific Computing Operations
//!
//! This module contains implementations for scientific domain-specific calculations
//! including chemistry, biology, thermodynamics, optics, geophysics, and engineering.
//!
//! Each function maps the unified API types to specialized scientific modules
//! and handles parameter conversion and result formatting.

use crate::engine::*;

/// Compute chemistry operations
///
/// Handles chemistry calculations including:
/// - pH calculations
/// - Buffer capacity
/// - Arrhenius equation
/// - Rate laws
/// - Gibbs free energy
/// - Nernst equation
/// - Beer-Lambert law
/// - Van der Waals equation
pub fn compute_chemistry(op: &ChemistryOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::chemistry::*;

    // Map operation to chemistry module operation
    let chem_op = match op {
        ChemistryOp::PhCalculation => ChemistryOperation::PhCalculation,
        ChemistryOp::BufferCapacity => ChemistryOperation::BufferCapacity,
        ChemistryOp::Arrhenius => ChemistryOperation::Arrhenius,
        ChemistryOp::RateLaw => ChemistryOperation::RateLaw,
        ChemistryOp::GibbsFreeEnergy => ChemistryOperation::GibbsFreeEnergy,
        ChemistryOp::NernstEquation => ChemistryOperation::NernstEquation,
        ChemistryOp::BeerLambert => ChemistryOperation::BeerLambert,
        ChemistryOp::VanDerWaals => ChemistryOperation::VanDerWaals,
        _ => {
            return Err(format!(
                "Legacy chemistry operation {:?} not implemented in new module",
                op
            ));
        }
    };

    // Convert parameters to ChemistryParams
    let params: ChemistryParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse chemistry parameters: {}", e))?;

    let chem_input = ChemistryInput {
        operation: chem_op,
        parameters: params,
    };

    let result = calculate_chemistry(chem_input)
        .map_err(|e| format!("Chemistry calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "interpretation": result.interpretation
        }),
        additional: None,
        metadata: None,
    })
}

/// Compute biology operations
///
/// Handles biological calculations including:
/// - Michaelis-Menten kinetics
/// - Pharmacokinetics
/// - Hardy-Weinberg equilibrium
/// - Goldman equation
/// - Allometric scaling
/// - Growth models
pub fn compute_biology(op: &BiologyOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::biology::*;

    let bio_op = match op {
        BiologyOp::MichaelisMenten => BiologyOperation::MichaelisMenten,
        BiologyOp::Pharmacokinetics => BiologyOperation::Pharmacokinetics,
        BiologyOp::HardyWeinberg => BiologyOperation::HardyWeinberg,
        BiologyOp::GoldmanEquation => BiologyOperation::GoldmanEquation,
        BiologyOp::AllometricScaling => BiologyOperation::AllometricScaling,
        BiologyOp::GrowthModel => BiologyOperation::LineweaverBurk,
    };

    let params: BiologyParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse biology parameters: {}", e))?;

    let bio_input = BiologyInput {
        operation: bio_op,
        parameters: params,
    };

    let result =
        calculate_biology(bio_input).map_err(|e| format!("Biology calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "interpretation": result.interpretation,
            "additional_data": result.additional_data
        }),
        additional: None,
        metadata: None,
    })
}

/// Compute thermodynamics operations
///
/// Handles thermodynamic calculations including:
/// - Heat conduction
/// - Convection
/// - Thermal radiation
/// - Thermal resistance
/// - Entropy calculations
pub fn compute_thermodynamics(
    op: &ThermodynamicsOp,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::thermodynamics::*;

    let operation = match op {
        ThermodynamicsOp::Conduction => ThermodynamicsOperation::Conduction,
        ThermodynamicsOp::Convection => ThermodynamicsOperation::Convection,
        ThermodynamicsOp::Radiation => ThermodynamicsOperation::Radiation,
        ThermodynamicsOp::ThermalResistance => ThermodynamicsOperation::ThermalResistance,
        ThermodynamicsOp::Entropy => ThermodynamicsOperation::Entropy,
    };

    let params: ThermodynamicsParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse thermodynamics parameters: {}", e))?;

    let thermo_input = ThermodynamicsInput {
        operation,
        parameters: params,
    };

    let result = calculate_thermodynamics(thermo_input)
        .map_err(|e| format!("Thermodynamics calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "interpretation": result.interpretation,
            "additional_info": result.additional_info
        }),
        additional: None,
        metadata: None,
    })
}

/// Compute optics operations
///
/// Handles optical physics calculations including:
/// - Thin lens equation
/// - Snell's law
/// - Diffraction grating
/// - Fresnel equations
pub fn compute_optics(op: &OpticsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::optics::*;

    let optics_op = match op {
        OpticsOp::ThinLens => OpticsOperation::ThinLens,
        OpticsOp::SnellsLaw => OpticsOperation::SnellsLaw,
        OpticsOp::DiffractionGrating => OpticsOperation::DiffractionGrating,
        OpticsOp::FresnelEquations => OpticsOperation::FresnelEquations,
    };

    let params: OpticsParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse optics parameters: {}", e))?;

    let optics_input = OpticsInput {
        operation: optics_op,
        parameters: params,
    };

    let result =
        calculate_optics(optics_input).map_err(|e| format!("Optics calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "primary_value": result.primary_value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "secondary_values": result.secondary_values,
            "interpretation": result.interpretation
        }),
        additional: None,
        metadata: None,
    })
}

/// Compute geophysics operations
///
/// Handles geophysics calculations including:
/// - Seismology
/// - Atmospheric science
/// - Radiometric dating
/// - Planetary science
pub fn compute_geophysics(op: &GeophysicsOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::geophysics::*;

    let category = match op {
        GeophysicsOp::Seismology => GeophysicsCategory::Seismology,
        GeophysicsOp::Atmosphere => GeophysicsCategory::Atmosphere,
        GeophysicsOp::RadiometricDating => GeophysicsCategory::RadiometricDating,
        GeophysicsOp::PlanetaryScience => GeophysicsCategory::PlanetaryScience,
    };

    let params: GeophysicsParams = serde_json::from_value(
        serde_json::to_value(&input.parameters)
            .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
    )
    .map_err(|e| format!("Failed to parse geophysics parameters: {}", e))?;

    let geo_input = GeophysicsInput {
        category,
        parameters: params,
    };

    let result = calculate_geophysics(geo_input)
        .map_err(|e| format!("Geophysics calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "uncertainty": result.uncertainty,
            "interpretation": result.interpretation,
            "additional_data": result.additional_data
        }),
        additional: None,
        metadata: None,
    })
}

/// Compute engineering operations
///
/// Handles engineering calculations organized by discipline:
/// - Acoustics: Sound pressure level, Doppler effect, Reverberation time
/// - Materials: Stress, Strain, Fracture mechanics
/// - Fluid Mechanics: Bernoulli, Poiseuille, Drag
/// - Control Theory: PID controller, First-order response
pub fn compute_engineering(op: &EngineeringOp, input: &ComputeInput) -> ToolResult<ComputeOutput> {
    use crate::engineering::*;

    // Engineering operations are organized by discipline, infer from operation
    let (discipline, params) = match op {
        EngineeringOp::SoundPressureLevel
        | EngineeringOp::DopplerEffect
        | EngineeringOp::ReverberationTime => (
            EngineeringDiscipline::Acoustics,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
        EngineeringOp::Stress | EngineeringOp::Strain | EngineeringOp::FractureMechanics => (
            EngineeringDiscipline::Materials,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
        EngineeringOp::Bernoulli | EngineeringOp::Poiseuille | EngineeringOp::Drag => (
            EngineeringDiscipline::FluidMechanics,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
        EngineeringOp::PidController | EngineeringOp::FirstOrderResponse => (
            EngineeringDiscipline::ControlTheory,
            serde_json::to_value(&input.parameters)
                .map_err(|e| format!("Failed to serialize parameters: {}", e))?,
        ),
    };

    let eng_params: EngineeringParams = serde_json::from_value(params)
        .map_err(|e| format!("Failed to parse engineering parameters: {}", e))?;

    let eng_input = EngineeringInput {
        discipline,
        parameters: eng_params,
    };

    let result = calculate_engineering(eng_input)
        .map_err(|e| format!("Engineering calculation error: {}", e))?;

    Ok(ComputeOutput {
        result: serde_json::json!({
            "value": result.value,
            "unit": result.unit,
            "formula_used": result.formula_used,
            "classification": result.classification,
            "interpretation": result.interpretation,
            "additional": result.additional
        }),
        additional: None,
        metadata: None,
    })
}
