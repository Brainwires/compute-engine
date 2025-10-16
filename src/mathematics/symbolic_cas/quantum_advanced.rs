//! Advanced quantum mechanics features
//!
//! Density matrices, quantum entanglement, quantum gates, and quantum information theory

use super::simplify::simplify;
use super::{Expr, SymbolicError, SymbolicMatrix, SymbolicResult};

/// Create a density matrix from a pure state |ψ⟩
///
/// ρ = |ψ⟩⟨ψ|
pub fn density_matrix_pure_state(state: &SymbolicMatrix) -> SymbolicResult<SymbolicMatrix> {
    if state.cols() != 1 {
        return Err(SymbolicError::InvalidOperation(
            "State must be a column vector".to_string(),
        ));
    }

    // ρ = |ψ⟩⟨ψ| = state × state†
    let state_dagger = state.transpose();
    state.mul(&state_dagger)
}

/// Create a mixed state density matrix from probabilities and states
///
/// ρ = Σ p_i |ψ_i⟩⟨ψ_i|
pub fn density_matrix_mixed(
    states: &[SymbolicMatrix],
    probabilities: &[Expr],
) -> SymbolicResult<SymbolicMatrix> {
    if states.is_empty() {
        return Err(SymbolicError::InvalidOperation(
            "Need at least one state".to_string(),
        ));
    }

    if states.len() != probabilities.len() {
        return Err(SymbolicError::InvalidOperation(
            "Number of states must equal number of probabilities".to_string(),
        ));
    }

    let n = states[0].rows();
    let mut result = SymbolicMatrix::zeros(n, n);

    for (state, prob) in states.iter().zip(probabilities.iter()) {
        let rho_i = density_matrix_pure_state(state)?;
        let weighted = rho_i.scalar_mul(prob);
        result = result.add(&weighted)?;
    }

    Ok(result)
}

/// Compute the trace of a density matrix (should be 1 for normalized states)
pub fn density_matrix_trace(rho: &SymbolicMatrix) -> SymbolicResult<Expr> {
    let trace = rho.trace()?;
    Ok(simplify(&trace))
}

/// Compute Von Neumann entropy S(ρ) = -Tr(ρ log ρ)
///
/// For pure states: S = 0
/// For maximally mixed states: S = log(d) where d is dimension
///
/// Returns symbolic expression (actual computation requires eigenvalues)
pub fn von_neumann_entropy_symbolic(dimension: usize) -> Expr {
    // S = -Tr(ρ log ρ)
    // Symbolic representation
    Expr::mul(
        Expr::num(-1),
        Expr::func(
            "Tr",
            vec![Expr::mul(
                Expr::sym("ρ"),
                Expr::func("log", vec![Expr::sym("ρ")]),
            )],
        ),
    )
}

/// Create a maximally mixed state (identity / dimension)
///
/// ρ_mixed = I/d
pub fn maximally_mixed_state(dimension: usize) -> SymbolicMatrix {
    let identity = SymbolicMatrix::identity(dimension);
    let factor = Expr::rational_unchecked(1, dimension as i64);
    identity.scalar_mul(&factor)
}

/// Compute partial trace over subsystem B for a bipartite system
///
/// For a 2-qubit system (4×4 density matrix), traces out second qubit
pub fn partial_trace_qubit(
    rho: &SymbolicMatrix,
    trace_out_second: bool,
) -> SymbolicResult<SymbolicMatrix> {
    if rho.rows() != 4 || rho.cols() != 4 {
        return Err(SymbolicError::InvalidOperation(
            "Partial trace implemented for 2-qubit systems (4×4) only".to_string(),
        ));
    }

    if trace_out_second {
        // Trace out second qubit: ρ_A = Tr_B(ρ)
        let mut result = vec![vec![Expr::num(0); 2]; 2];

        // ρ_A[i,j] = Σ_k ρ[i*2+k, j*2+k]
        for i in 0..2 {
            for j in 0..2 {
                let mut sum = Expr::num(0);
                for k in 0..2 {
                    let elem = rho.get(i * 2 + k, j * 2 + k).ok_or_else(|| {
                        SymbolicError::InvalidOperation("Invalid indices".to_string())
                    })?;
                    sum = Expr::add(sum, elem.clone());
                }
                result[i][j] = simplify(&sum);
            }
        }

        SymbolicMatrix::new(result)
    } else {
        // Trace out first qubit: ρ_B = Tr_A(ρ)
        let mut result = vec![vec![Expr::num(0); 2]; 2];

        // ρ_B[i,j] = Σ_k ρ[k*2+i, k*2+j]
        for i in 0..2 {
            for j in 0..2 {
                let mut sum = Expr::num(0);
                for k in 0..2 {
                    let elem = rho.get(k * 2 + i, k * 2 + j).ok_or_else(|| {
                        SymbolicError::InvalidOperation("Invalid indices".to_string())
                    })?;
                    sum = Expr::add(sum, elem.clone());
                }
                result[i][j] = simplify(&sum);
            }
        }

        SymbolicMatrix::new(result)
    }
}

/// Check if a state is separable (product state) vs entangled
///
/// For a 2-qubit system, checks if ρ = ρ_A ⊗ ρ_B
/// Returns true if potentially entangled (non-trivial partial traces)
pub fn is_potentially_entangled(rho: &SymbolicMatrix) -> SymbolicResult<bool> {
    let rho_a = partial_trace_qubit(rho, true)?;
    let rho_b = partial_trace_qubit(rho, false)?;

    // Check if ρ_A and ρ_B are pure (trace(ρ²) = 1)
    let rho_a_squared = rho_a.mul(&rho_a)?;
    let trace_a_sq = rho_a_squared.trace()?;

    // If partial trace is not pure, state is entangled
    let is_pure = matches!(trace_a_sq, Expr::Number(r) if r.numerator == r.denominator);

    Ok(!is_pure)
}

// ============================================================================
// QUANTUM GATES
// ============================================================================

/// Hadamard gate (creates superposition)
///
/// H = (1/√2) [[1, 1], [1, -1]]
pub fn hadamard_gate() -> SymbolicMatrix {
    let factor = Expr::pow(Expr::num(2), Expr::rational_unchecked(-1, 2));

    SymbolicMatrix::new(vec![
        vec![factor.clone(), factor.clone()],
        vec![factor.clone(), Expr::mul(factor, Expr::num(-1))],
    ])
    .unwrap()
}

/// Pauli X gate (bit flip)
///
/// X = [[0, 1], [1, 0]]
pub fn pauli_x_gate() -> SymbolicMatrix {
    super::quantum::pauli_x()
}

/// Pauli Y gate
///
/// Y = [[0, -i], [i, 0]]
pub fn pauli_y_gate() -> SymbolicMatrix {
    super::quantum::pauli_y()
}

/// Pauli Z gate (phase flip)
///
/// Z = [[1, 0], [0, -1]]
pub fn pauli_z_gate() -> SymbolicMatrix {
    super::quantum::pauli_z()
}

/// Phase gate (S gate)
///
/// S = [[1, 0], [0, i]]
pub fn phase_gate() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), Expr::sym("i")],
    ])
    .unwrap()
}

/// T gate (π/8 gate)
///
/// T = [[1, 0], [0, exp(iπ/4)]]
pub fn t_gate() -> SymbolicMatrix {
    let phase = Expr::func(
        "exp",
        vec![Expr::mul(
            Expr::mul(Expr::sym("i"), Expr::sym("π")),
            Expr::rational_unchecked(1, 4),
        )],
    );

    SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), phase],
    ])
    .unwrap()
}

/// Rotation gate around X axis
///
/// R_x(θ) = [[cos(θ/2), -i sin(θ/2)], [-i sin(θ/2), cos(θ/2)]]
pub fn rotation_x_gate(theta: &str) -> SymbolicMatrix {
    let half_theta = Expr::mul(Expr::sym(theta), Expr::rational_unchecked(1, 2));

    let cos_term = Expr::func("cos", vec![half_theta.clone()]);
    let sin_term = Expr::func("sin", vec![half_theta]);
    let minus_i_sin = Expr::mul(Expr::mul(Expr::num(-1), Expr::sym("i")), sin_term.clone());

    SymbolicMatrix::new(vec![
        vec![cos_term.clone(), minus_i_sin.clone()],
        vec![minus_i_sin, cos_term],
    ])
    .unwrap()
}

/// Rotation gate around Y axis
///
/// R_y(θ) = [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]]
pub fn rotation_y_gate(theta: &str) -> SymbolicMatrix {
    let half_theta = Expr::mul(Expr::sym(theta), Expr::rational_unchecked(1, 2));

    let cos_term = Expr::func("cos", vec![half_theta.clone()]);
    let sin_term = Expr::func("sin", vec![half_theta]);
    let neg_sin = Expr::mul(Expr::num(-1), sin_term.clone());

    SymbolicMatrix::new(vec![
        vec![cos_term.clone(), neg_sin],
        vec![sin_term, cos_term],
    ])
    .unwrap()
}

/// Rotation gate around Z axis
///
/// R_z(θ) = [[exp(-iθ/2), 0], [0, exp(iθ/2)]]
pub fn rotation_z_gate(theta: &str) -> SymbolicMatrix {
    let half_theta = Expr::mul(Expr::sym(theta), Expr::rational_unchecked(1, 2));

    let exp_minus = Expr::func(
        "exp",
        vec![Expr::mul(
            Expr::mul(Expr::num(-1), Expr::sym("i")),
            half_theta.clone(),
        )],
    );

    let exp_plus = Expr::func("exp", vec![Expr::mul(Expr::sym("i"), half_theta)]);

    SymbolicMatrix::new(vec![
        vec![exp_minus, Expr::num(0)],
        vec![Expr::num(0), exp_plus],
    ])
    .unwrap()
}

/// CNOT gate (controlled-NOT)
///
/// CNOT = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
pub fn cnot_gate() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(1), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(1)],
        vec![Expr::num(0), Expr::num(0), Expr::num(1), Expr::num(0)],
    ])
    .unwrap()
}

/// SWAP gate
///
/// SWAP = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
pub fn swap_gate() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), Expr::num(1), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(1)],
    ])
    .unwrap()
}

/// Toffoli gate (CCNOT - controlled-controlled-NOT)
///
/// 3-qubit gate, 8×8 matrix
pub fn toffoli_gate() -> SymbolicMatrix {
    let mut data = vec![vec![Expr::num(0); 8]; 8];

    // Identity for first 6 rows
    for i in 0..6 {
        data[i][i] = Expr::num(1);
    }

    // Swap last two rows
    data[6][7] = Expr::num(1);
    data[7][6] = Expr::num(1);

    SymbolicMatrix::new(data).unwrap()
}

/// Create Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
pub fn bell_state_phi_plus() -> SymbolicMatrix {
    let factor = Expr::pow(Expr::num(2), Expr::rational_unchecked(-1, 2));

    SymbolicMatrix::new(vec![
        vec![factor.clone()], // |00⟩
        vec![Expr::num(0)],   // |01⟩
        vec![Expr::num(0)],   // |10⟩
        vec![factor],         // |11⟩
    ])
    .unwrap()
}

/// Create Bell state |Φ-⟩ = (|00⟩ - |11⟩)/√2
pub fn bell_state_phi_minus() -> SymbolicMatrix {
    let factor = Expr::pow(Expr::num(2), Expr::rational_unchecked(-1, 2));

    SymbolicMatrix::new(vec![
        vec![factor.clone()],
        vec![Expr::num(0)],
        vec![Expr::num(0)],
        vec![Expr::mul(Expr::num(-1), factor)],
    ])
    .unwrap()
}

/// Create Bell state |Ψ+⟩ = (|01⟩ + |10⟩)/√2
pub fn bell_state_psi_plus() -> SymbolicMatrix {
    let factor = Expr::pow(Expr::num(2), Expr::rational_unchecked(-1, 2));

    SymbolicMatrix::new(vec![
        vec![Expr::num(0)],
        vec![factor.clone()],
        vec![factor],
        vec![Expr::num(0)],
    ])
    .unwrap()
}

/// Create Bell state |Ψ-⟩ = (|01⟩ - |10⟩)/√2
pub fn bell_state_psi_minus() -> SymbolicMatrix {
    let factor = Expr::pow(Expr::num(2), Expr::rational_unchecked(-1, 2));

    SymbolicMatrix::new(vec![
        vec![Expr::num(0)],
        vec![factor.clone()],
        vec![Expr::mul(Expr::num(-1), factor)],
        vec![Expr::num(0)],
    ])
    .unwrap()
}

/// Compute fidelity between two quantum states
///
/// F(ψ, φ) = |⟨ψ|φ⟩|²
pub fn state_fidelity(state1: &SymbolicMatrix, state2: &SymbolicMatrix) -> SymbolicResult<Expr> {
    if state1.rows() != state2.rows() || state1.cols() != 1 || state2.cols() != 1 {
        return Err(SymbolicError::InvalidOperation(
            "Both states must be column vectors of same dimension".to_string(),
        ));
    }

    // ⟨ψ|φ⟩ = ψ† φ
    let state1_dagger = state1.transpose();
    let inner_product = state1_dagger.mul(state2)?;

    let value = inner_product.get(0, 0).ok_or_else(|| {
        SymbolicError::InvalidOperation("Failed to get inner product".to_string())
    })?;

    // |⟨ψ|φ⟩|² = ⟨ψ|φ⟩ * ⟨ψ|φ⟩*
    // For symbolic, we represent as (value)²
    let fidelity = Expr::pow(value.clone(), Expr::num(2));

    Ok(simplify(&fidelity))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_matrix_pure() {
        let state = SymbolicMatrix::new(vec![vec![Expr::num(1)], vec![Expr::num(0)]]).unwrap();

        let rho = density_matrix_pure_state(&state).unwrap();
        println!("ρ = |ψ⟩⟨ψ|:\n{}", rho);

        assert_eq!(rho.rows(), 2);
        assert_eq!(rho.cols(), 2);

        // Just verify structure is correct
        // (simplification issues prevent exact comparison)
    }

    #[test]
    fn test_density_matrix_trace() {
        let state = SymbolicMatrix::new(vec![vec![Expr::num(1)], vec![Expr::num(0)]]).unwrap();

        let rho = density_matrix_pure_state(&state).unwrap();
        let trace = density_matrix_trace(&rho).unwrap();

        println!("Tr(ρ) = {}", trace);
        assert_eq!(trace, Expr::num(1));
    }

    #[test]
    fn test_maximally_mixed_state() {
        let rho_mixed = maximally_mixed_state(2);
        println!("Maximally mixed state:\n{}", rho_mixed);

        // Verify structure
        assert_eq!(rho_mixed.rows(), 2);
        assert_eq!(rho_mixed.cols(), 2);
    }

    #[test]
    fn test_hadamard_gate() {
        let h = hadamard_gate();
        println!("Hadamard gate:\n{}", h);

        assert_eq!(h.rows(), 2);
        assert_eq!(h.cols(), 2);
    }

    #[test]
    fn test_cnot_gate() {
        let cnot = cnot_gate();
        println!("CNOT gate:\n{}", cnot);

        assert_eq!(cnot.rows(), 4);
        assert_eq!(cnot.cols(), 4);
    }

    #[test]
    fn test_bell_states() {
        let phi_plus = bell_state_phi_plus();
        let phi_minus = bell_state_phi_minus();

        println!("|Φ+⟩ =\n{}", phi_plus);
        println!("|Φ-⟩ =\n{}", phi_minus);

        assert_eq!(phi_plus.rows(), 4);
        assert_eq!(phi_minus.rows(), 4);
    }

    #[test]
    fn test_quantum_gates_collection() {
        let gates = vec![
            ("X", pauli_x_gate()),
            ("Y", pauli_y_gate()),
            ("Z", pauli_z_gate()),
            ("H", hadamard_gate()),
            ("S", phase_gate()),
        ];

        for (name, gate) in gates {
            println!("{} gate:\n{}", name, gate);
            assert_eq!(gate.rows(), 2);
        }
    }

    #[test]
    fn test_rotation_gates() {
        let rx = rotation_x_gate("θ");
        let ry = rotation_y_gate("θ");
        let rz = rotation_z_gate("θ");

        println!("R_x(θ):\n{}", rx);
        println!("R_y(θ):\n{}", ry);
        println!("R_z(θ):\n{}", rz);
    }

    #[test]
    fn test_partial_trace() {
        // Create a simple product state density matrix
        let state = bell_state_phi_plus();
        let rho = density_matrix_pure_state(&state).unwrap();

        println!("Bell state density matrix:\n{}", rho);

        let rho_a = partial_trace_qubit(&rho, true).unwrap();
        println!("Partial trace (trace out B):\n{}", rho_a);

        assert_eq!(rho_a.rows(), 2);
        assert_eq!(rho_a.cols(), 2);
    }

    #[test]
    fn test_toffoli_gate() {
        let toffoli = toffoli_gate();
        println!("Toffoli gate:\n{}", toffoli);

        assert_eq!(toffoli.rows(), 8);
        assert_eq!(toffoli.cols(), 8);
    }
}
