//! Quantum mechanics operations
//!
//! Provides symbolic computation for quantum mechanics including commutators,
//! quantum matrices (Pauli, Dirac), operators, and quantum state manipulations.

use super::simplify::simplify;
use super::{Expr, SymbolicError, SymbolicMatrix, SymbolicResult};
use std::collections::HashMap;

/// Compute the commutator [A, B] = AB - BA
///
/// # Arguments
/// * `a` - First matrix operator
/// * `b` - Second matrix operator
///
/// # Returns
/// The commutator [A, B]
pub fn commutator(a: &SymbolicMatrix, b: &SymbolicMatrix) -> SymbolicResult<SymbolicMatrix> {
    // AB
    let ab = a.mul(b)?;

    // BA
    let ba = b.mul(a)?;

    // AB - BA
    let result = ab.add(&ba.scalar_mul(&Expr::num(-1)))?;

    // Simplify each component
    let mut simplified_data = vec![vec![Expr::num(0); result.cols()]; result.rows()];
    for i in 0..result.rows() {
        for j in 0..result.cols() {
            let elem = result.get(i, j).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid matrix indices".to_string())
            })?;
            simplified_data[i][j] = simplify(elem);
        }
    }

    SymbolicMatrix::new(simplified_data)
}

/// Compute the anticommutator {A, B} = AB + BA
///
/// # Arguments
/// * `a` - First matrix operator
/// * `b` - Second matrix operator
///
/// # Returns
/// The anticommutator {A, B}
pub fn anticommutator(a: &SymbolicMatrix, b: &SymbolicMatrix) -> SymbolicResult<SymbolicMatrix> {
    // AB
    let ab = a.mul(b)?;

    // BA
    let ba = b.mul(a)?;

    // AB + BA
    let result = ab.add(&ba)?;

    // Simplify each component
    let mut simplified_data = vec![vec![Expr::num(0); result.cols()]; result.rows()];
    for i in 0..result.rows() {
        for j in 0..result.cols() {
            let elem = result.get(i, j).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid matrix indices".to_string())
            })?;
            simplified_data[i][j] = simplify(elem);
        }
    }

    SymbolicMatrix::new(simplified_data)
}

/// Create the Pauli X matrix (σ_x)
///
/// σ_x = [[0, 1], [1, 0]]
pub fn pauli_x() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(0), Expr::num(1)],
        vec![Expr::num(1), Expr::num(0)],
    ])
    .unwrap()
}

/// Create the Pauli Y matrix (σ_y)
///
/// σ_y = [[0, -i], [i, 0]]
pub fn pauli_y() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(0), Expr::mul(Expr::num(-1), Expr::sym("i"))],
        vec![Expr::sym("i"), Expr::num(0)],
    ])
    .unwrap()
}

/// Create the Pauli Z matrix (σ_z)
///
/// σ_z = [[1, 0], [0, -1]]
pub fn pauli_z() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), Expr::num(-1)],
    ])
    .unwrap()
}

/// Get all three Pauli matrices
pub fn pauli_matrices() -> [SymbolicMatrix; 3] {
    [pauli_x(), pauli_y(), pauli_z()]
}

/// Create a Dirac gamma matrix (γ^0) for 4D Dirac equation
///
/// γ^0 = [[I_2, 0], [0, -I_2]]
pub fn dirac_gamma_0() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(1), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(-1), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(-1)],
    ])
    .unwrap()
}

/// Create a Dirac gamma matrix (γ^1)
///
/// γ^1 = [[0, σ_x], [-σ_x, 0]]
pub fn dirac_gamma_1() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(1)],
        vec![Expr::num(0), Expr::num(0), Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), Expr::num(-1), Expr::num(0), Expr::num(0)],
        vec![Expr::num(-1), Expr::num(0), Expr::num(0), Expr::num(0)],
    ])
    .unwrap()
}

/// Create a Dirac gamma matrix (γ^2)
///
/// γ^2 = [[0, σ_y], [-σ_y, 0]]
pub fn dirac_gamma_2() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![
            Expr::num(0),
            Expr::num(0),
            Expr::num(0),
            Expr::mul(Expr::num(-1), Expr::sym("i")),
        ],
        vec![Expr::num(0), Expr::num(0), Expr::sym("i"), Expr::num(0)],
        vec![
            Expr::num(0),
            Expr::mul(Expr::num(-1), Expr::sym("i")),
            Expr::num(0),
            Expr::num(0),
        ],
        vec![Expr::sym("i"), Expr::num(0), Expr::num(0), Expr::num(0)],
    ])
    .unwrap()
}

/// Create a Dirac gamma matrix (γ^3)
///
/// γ^3 = [[0, σ_z], [-σ_z, 0]]
pub fn dirac_gamma_3() -> SymbolicMatrix {
    SymbolicMatrix::new(vec![
        vec![Expr::num(0), Expr::num(0), Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(-1)],
        vec![Expr::num(-1), Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(1), Expr::num(0), Expr::num(0)],
    ])
    .unwrap()
}

/// Get all four Dirac gamma matrices
pub fn dirac_gamma_matrices() -> [SymbolicMatrix; 4] {
    [
        dirac_gamma_0(),
        dirac_gamma_1(),
        dirac_gamma_2(),
        dirac_gamma_3(),
    ]
}

/// Create the angular momentum operator L_x (for spin)
///
/// L_x = (ℏ/2) σ_x
pub fn angular_momentum_x(hbar: Option<Expr>) -> SymbolicMatrix {
    let h = hbar.unwrap_or_else(|| Expr::sym("ℏ"));
    let factor = Expr::mul(h, Expr::rational_unchecked(1, 2));
    pauli_x().scalar_mul(&factor)
}

/// Create the angular momentum operator L_y (for spin)
///
/// L_y = (ℏ/2) σ_y
pub fn angular_momentum_y(hbar: Option<Expr>) -> SymbolicMatrix {
    let h = hbar.unwrap_or_else(|| Expr::sym("ℏ"));
    let factor = Expr::mul(h, Expr::rational_unchecked(1, 2));
    pauli_y().scalar_mul(&factor)
}

/// Create the angular momentum operator L_z (for spin)
///
/// L_z = (ℏ/2) σ_z
pub fn angular_momentum_z(hbar: Option<Expr>) -> SymbolicMatrix {
    let h = hbar.unwrap_or_else(|| Expr::sym("ℏ"));
    let factor = Expr::mul(h, Expr::rational_unchecked(1, 2));
    pauli_z().scalar_mul(&factor)
}

/// Create the quantum harmonic oscillator ladder operator (creation operator a†)
///
/// For position operator x and momentum operator p:
/// a† = (1/√(2mωℏ)) (mωx - ip)
///
/// Returns the symbolic representation
pub fn creation_operator_symbolic() -> Expr {
    // (1/√(2mωℏ)) (mωx - ip)
    let prefactor = Expr::pow(
        Expr::mul(
            Expr::mul(Expr::num(2), Expr::sym("m")),
            Expr::mul(Expr::sym("ω"), Expr::sym("ℏ")),
        ),
        Expr::rational_unchecked(-1, 2),
    );

    let term = Expr::add(
        Expr::mul(Expr::mul(Expr::sym("m"), Expr::sym("ω")), Expr::sym("x")),
        Expr::mul(Expr::mul(Expr::num(-1), Expr::sym("i")), Expr::sym("p")),
    );

    Expr::mul(prefactor, term)
}

/// Create the quantum harmonic oscillator ladder operator (annihilation operator a)
///
/// a = (1/√(2mωℏ)) (mωx + ip)
pub fn annihilation_operator_symbolic() -> Expr {
    // (1/√(2mωℏ)) (mωx + ip)
    let prefactor = Expr::pow(
        Expr::mul(
            Expr::mul(Expr::num(2), Expr::sym("m")),
            Expr::mul(Expr::sym("ω"), Expr::sym("ℏ")),
        ),
        Expr::rational_unchecked(-1, 2),
    );

    let term = Expr::add(
        Expr::mul(Expr::mul(Expr::sym("m"), Expr::sym("ω")), Expr::sym("x")),
        Expr::mul(Expr::sym("i"), Expr::sym("p")),
    );

    Expr::mul(prefactor, term)
}

/// Compute the time evolution operator U(t) = exp(-iHt/ℏ)
///
/// For a Hamiltonian H and time t, returns the symbolic expression
pub fn time_evolution_operator(hamiltonian: &str, time: &str) -> Expr {
    // exp(-iHt/ℏ)
    let exponent = Expr::mul(
        Expr::mul(
            Expr::mul(Expr::num(-1), Expr::sym("i")),
            Expr::sym(hamiltonian),
        ),
        Expr::mul(Expr::sym(time), Expr::pow(Expr::sym("ℏ"), Expr::num(-1))),
    );

    Expr::func("exp", vec![exponent])
}

/// Compute expectation value <ψ|A|ψ>
///
/// # Arguments
/// * `state` - Quantum state vector (column vector)
/// * `operator` - Operator matrix
///
/// # Returns
/// The expectation value as a symbolic expression
pub fn expectation_value(
    state: &SymbolicMatrix,
    operator: &SymbolicMatrix,
) -> SymbolicResult<Expr> {
    // Check state is a column vector
    if state.cols() != 1 {
        return Err(SymbolicError::InvalidOperation(
            "State must be a column vector".to_string(),
        ));
    }

    // A|ψ>
    let a_psi = operator.mul(state)?;

    // <ψ| = |ψ>†
    let psi_dagger = state.transpose();

    // <ψ|A|ψ>
    let result = psi_dagger.mul(&a_psi)?;

    // Should be a 1x1 matrix
    if result.rows() != 1 || result.cols() != 1 {
        return Err(SymbolicError::InvalidOperation(
            "Expectation value should be a scalar".to_string(),
        ));
    }

    let value = result.get(0, 0).ok_or_else(|| {
        SymbolicError::InvalidOperation("Failed to get expectation value".to_string())
    })?;

    Ok(simplify(value))
}

/// Compute the uncertainty relation ΔA·ΔB ≥ (1/2)|<[A,B]>|
///
/// Returns the commutator [A, B] for uncertainty analysis
pub fn uncertainty_commutator(
    a: &SymbolicMatrix,
    b: &SymbolicMatrix,
) -> SymbolicResult<SymbolicMatrix> {
    commutator(a, b)
}

/// Verify Pauli matrix properties
///
/// Returns a map of property names to boolean verification results
pub fn verify_pauli_properties() -> HashMap<String, bool> {
    let mut results = HashMap::new();

    let sigma_x = pauli_x();
    let sigma_y = pauli_y();
    let sigma_z = pauli_z();
    let identity = SymbolicMatrix::identity(2);

    // σ_x² = I
    let sigma_x_squared = sigma_x.mul(&sigma_x).unwrap();
    results.insert(
        "σ_x² = I".to_string(),
        matrices_equal(&sigma_x_squared, &identity),
    );

    // σ_y² = I
    let sigma_y_squared = sigma_y.mul(&sigma_y).unwrap();
    results.insert(
        "σ_y² = I".to_string(),
        matrices_equal(&sigma_y_squared, &identity),
    );

    // σ_z² = I
    let sigma_z_squared = sigma_z.mul(&sigma_z).unwrap();
    results.insert(
        "σ_z² = I".to_string(),
        matrices_equal(&sigma_z_squared, &identity),
    );

    // [σ_x, σ_y] = 2i σ_z
    let comm_xy = commutator(&sigma_x, &sigma_y).unwrap();
    let expected_xy = sigma_z.scalar_mul(&Expr::mul(Expr::num(2), Expr::sym("i")));
    results.insert(
        "[σ_x, σ_y] = 2iσ_z".to_string(),
        matrices_equal(&comm_xy, &expected_xy),
    );

    results
}

/// Helper function to check if two matrices are symbolically equal
fn matrices_equal(a: &SymbolicMatrix, b: &SymbolicMatrix) -> bool {
    if a.rows() != b.rows() || a.cols() != b.cols() {
        return false;
    }

    for i in 0..a.rows() {
        for j in 0..a.cols() {
            let elem_a = simplify(a.get(i, j).unwrap());
            let elem_b = simplify(b.get(i, j).unwrap());
            if elem_a != elem_b {
                return false;
            }
        }
    }

    true
}

