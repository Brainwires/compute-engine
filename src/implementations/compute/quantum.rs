//! Quantum mechanics computation operations
//!
//! This module handles quantum mechanical calculations including:
//! - Fundamental quantum systems (harmonic oscillator, hydrogen atom)
//! - Angular momentum and spin operators
//! - Quantum states and density matrices
//! - Quantum information (entanglement, entropy, coherence)
//! - Quantum computing (gates, circuits, tomography)

use crate::engine::*;

/// Compute quantum mechanics operations
///
/// Handles quantum mechanical calculations:
///
/// **Fundamental Systems:**
/// - Schrodinger equation: Solve for wavefunctions and energies
/// - Harmonic oscillator: Energy levels and wavefunctions
/// - Hydrogen atom: Atomic orbitals and energy levels
///
/// **Angular Momentum & Spin:**
/// - Angular momentum: Eigenvalues and ladder operators
/// - Spin operators: Pauli matrices and spin states
///
/// **Perturbation Theory:**
/// - First and second-order energy corrections
/// - Matrix element calculations
///
/// **Quantum Tunneling:**
/// - Barrier penetration probability
/// - WKB approximation
///
/// **Quantum States:**
/// - Density matrices: Pure and mixed states
/// - State tomography: Reconstruct quantum states
///
/// **Quantum Information:**
/// - Entanglement measures: Concurrence, negativity
/// - Quantum entropy: Von Neumann and Renyi entropy
/// - Quantum coherence: Basis-dependent coherence measures
/// - Bell inequalities: Test for quantum correlations
///
/// **Quantum Computing:**
/// - Quantum gates: Hadamard, Pauli, phase, CNOT
/// - Quantum circuits: Multi-qubit gate sequences
pub fn compute_quantum_mechanics(
    op: &QuantumMechOp,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::physics::quantum_mechanics::*;

    let result_json = match op {
        QuantumMechOp::SchrodingerEquation => {
            let potential = input
                .parameters
                .get("potential")
                .and_then(|v| v.as_str())
                .ok_or("potential parameter required")?;
            let energy = input
                .parameters
                .get("energy")
                .and_then(|v| v.as_f64())
                .ok_or("energy parameter required")?;
            let position = input
                .parameters
                .get("position")
                .and_then(|v| v.as_f64())
                .ok_or("position parameter required")?;

            let mut parameters = std::collections::HashMap::new();
            if let Some(val) = input.parameters.get("width").and_then(|v| v.as_f64()) {
                parameters.insert("width".to_string(), val);
            }
            if let Some(val) = input
                .parameters
                .get("quantum_number")
                .and_then(|v| v.as_f64())
            {
                parameters.insert("quantum_number".to_string(), val);
            }
            if let Some(val) = input.parameters.get("omega").and_then(|v| v.as_f64()) {
                parameters.insert("omega".to_string(), val);
            }
            if let Some(val) = input.parameters.get("mass").and_then(|v| v.as_f64()) {
                parameters.insert("mass".to_string(), val);
            }
            if let Some(val) = input
                .parameters
                .get("barrier_height")
                .and_then(|v| v.as_f64())
            {
                parameters.insert("barrier_height".to_string(), val);
            }

            let result = schrodinger_equation(SchrodingerRequest {
                potential: potential.to_string(),
                energy,
                position,
                parameters,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::HarmonicOscillator => {
            let quantum_number = input
                .parameters
                .get("quantum_number")
                .and_then(|v| v.as_u64())
                .ok_or("quantum_number parameter required")?
                as usize;
            let omega = input
                .parameters
                .get("omega")
                .and_then(|v| v.as_f64())
                .ok_or("omega parameter required")?;
            let mass = input
                .parameters
                .get("mass")
                .and_then(|v| v.as_f64())
                .ok_or("mass parameter required")?;
            let position = input.parameters.get("position").and_then(|v| v.as_f64());

            let result = harmonic_oscillator(HarmonicOscillatorRequest {
                quantum_number,
                omega,
                mass,
                position,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::HydrogenAtom => {
            let n = input
                .parameters
                .get("n")
                .and_then(|v| v.as_u64())
                .ok_or("n (principal quantum number) required")? as usize;
            let l = input
                .parameters
                .get("l")
                .and_then(|v| v.as_u64())
                .ok_or("l (orbital quantum number) required")? as usize;
            let m = input
                .parameters
                .get("m")
                .and_then(|v| v.as_i64())
                .ok_or("m (magnetic quantum number) required")? as i32;
            let r = input.parameters.get("r").and_then(|v| v.as_f64());

            let result = hydrogen_atom(HydrogenAtomRequest { n, l, m, r })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::AngularMomentum => {
            let l = input
                .parameters
                .get("l")
                .and_then(|v| v.as_u64())
                .ok_or("l parameter required")? as usize;
            let m = input
                .parameters
                .get("m")
                .and_then(|v| v.as_i64())
                .ok_or("m parameter required")? as i32;
            let operation = input
                .parameters
                .get("operation")
                .and_then(|v| v.as_str())
                .unwrap_or("eigenvalue");

            let result = angular_momentum(AngularMomentumRequest {
                l,
                m,
                operation: operation.to_string(),
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::SpinOperators => {
            let spin = input
                .parameters
                .get("spin")
                .and_then(|v| v.as_f64())
                .ok_or("spin parameter required")?;
            let component = input
                .parameters
                .get("component")
                .and_then(|v| v.as_str())
                .ok_or("component parameter required (x, y, z, plus, minus)")?;
            let state: Option<Vec<f64>> = input
                .parameters
                .get("state")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let result = spin_operators(SpinRequest {
                spin,
                component: component.to_string(),
                state,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::PerturbationTheory => {
            let order = input
                .parameters
                .get("order")
                .and_then(|v| v.as_u64())
                .ok_or("order parameter required (1 or 2)")? as usize;
            let unperturbed_energy = input
                .parameters
                .get("unperturbed_energy")
                .and_then(|v| v.as_f64())
                .ok_or("unperturbed_energy required")?;
            let perturbation_matrix_element = input
                .parameters
                .get("perturbation_matrix_element")
                .and_then(|v| v.as_f64())
                .ok_or("perturbation_matrix_element required")?;
            let energy_differences: Option<Vec<f64>> = input
                .parameters
                .get("energy_differences")
                .and_then(|v| serde_json::from_value(v.clone()).ok());
            let coupling_matrix_elements: Option<Vec<f64>> = input
                .parameters
                .get("coupling_matrix_elements")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let result = perturbation_theory(PerturbationRequest {
                order,
                unperturbed_energy,
                perturbation_matrix_element,
                energy_differences,
                coupling_matrix_elements,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::TunnelingProbability => {
            let barrier_height = input
                .parameters
                .get("barrier_height")
                .and_then(|v| v.as_f64())
                .ok_or("barrier_height required")?;
            let barrier_width = input
                .parameters
                .get("barrier_width")
                .and_then(|v| v.as_f64())
                .ok_or("barrier_width required")?;
            let particle_energy = input
                .parameters
                .get("particle_energy")
                .and_then(|v| v.as_f64())
                .ok_or("particle_energy required")?;
            let particle_mass = input
                .parameters
                .get("particle_mass")
                .and_then(|v| v.as_f64())
                .ok_or("particle_mass required")?;

            let result = tunneling_probability(TunnelingRequest {
                barrier_height,
                barrier_width,
                particle_energy,
                particle_mass,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::DensityMatrix => {
            let state_type = input
                .parameters
                .get("state_type")
                .and_then(|v| v.as_str())
                .ok_or("state_type required (pure or mixed)")?;
            let state_vector: Option<Vec<f64>> = input
                .parameters
                .get("state_vector")
                .and_then(|v| serde_json::from_value(v.clone()).ok());
            let density_mat: Option<Vec<Vec<f64>>> = input
                .parameters
                .get("density_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok());

            let result = density_matrix(DensityMatrixRequest {
                state_type: state_type.to_string(),
                state_vector,
                density_matrix: density_mat,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::EntanglementMeasure => {
            let state_vector: Vec<f64> = input
                .parameters
                .get("state_vector")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("state_vector required (4 components for two-qubit state)")?;
            let measure = input
                .parameters
                .get("measure")
                .and_then(|v| v.as_str())
                .unwrap_or("concurrence");

            let result = entanglement_measure(EntanglementRequest {
                state_vector,
                measure: measure.to_string(),
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::QuantumEntropy => {
            let density_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("density_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("density_matrix required")?;
            let entropy_type = input
                .parameters
                .get("entropy_type")
                .and_then(|v| v.as_str())
                .unwrap_or("von_neumann");
            let alpha = input.parameters.get("alpha").and_then(|v| v.as_f64());

            let result = quantum_entropy(QuantumEntropyRequest {
                density_matrix,
                entropy_type: entropy_type.to_string(),
                alpha,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::QuantumCoherence => {
            let density_matrix: Vec<Vec<f64>> = input
                .parameters
                .get("density_matrix")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("density_matrix required")?;
            let basis = input
                .parameters
                .get("basis")
                .and_then(|v| v.as_str())
                .unwrap_or("computational");

            let result = quantum_coherence(CoherenceRequest {
                density_matrix,
                basis: basis.to_string(),
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::BellInequality => {
            let state_vector: Vec<f64> = input
                .parameters
                .get("state_vector")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("state_vector required (two-qubit state)")?;
            let measurement_angles: Vec<f64> = input
                .parameters
                .get("measurement_angles")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("measurement_angles required (4 angles for CHSH)")?;

            let result = bell_inequality(BellInequalityRequest {
                state_vector,
                measurement_angles,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::QuantumTomography => {
            let measurements: Vec<Vec<f64>> = input
                .parameters
                .get("measurements")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("measurements required")?;
            let measurement_bases: Vec<String> = input
                .parameters
                .get("measurement_bases")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("measurement_bases required (e.g., [X, Y, Z])")?;

            let result = quantum_tomography(QuantumTomographyRequest {
                measurements,
                measurement_bases,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::QuantumGate => {
            let gate_type = input
                .parameters
                .get("gate_type")
                .and_then(|v| v.as_str())
                .ok_or("gate_type required (hadamard, pauli_x, pauli_z, phase, cnot)")?;
            let input_state: Vec<f64> = input
                .parameters
                .get("input_state")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("input_state required")?;
            let phase_angle = input.parameters.get("phase_angle").and_then(|v| v.as_f64());

            let result = quantum_gate(QuantumGateRequest {
                gate_type: gate_type.to_string(),
                input_state,
                phase_angle,
            })?;

            serde_json::to_value(result)
                .map_err(|e| format!("Failed to serialize result: {}", e))?
        }

        QuantumMechOp::QuantumCircuit => {
            let num_qubits = input
                .parameters
                .get("num_qubits")
                .and_then(|v| v.as_u64())
                .ok_or("num_qubits required")? as usize;
            let gates: Vec<serde_json::Value> = input
                .parameters
                .get("gates")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("gates required (list of gate operations)")?;
            let initial_state: Vec<f64> = input
                .parameters
                .get("initial_state")
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .ok_or("initial_state required")?;

            let result = quantum_circuit(QuantumCircuitRequest {
                num_qubits,
                gates,
                initial_state,
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
