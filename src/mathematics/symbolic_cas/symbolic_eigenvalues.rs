//! Symbolic eigenvalue and eigenvector computation
//!
//! Computes characteristic polynomials and eigenvalues symbolically

use super::{Expr, SymbolicError, SymbolicMatrix, SymbolicResult};

/// Compute the characteristic polynomial of a matrix
/// det(A - λI) where λ is a symbol
pub fn characteristic_polynomial(matrix: &SymbolicMatrix) -> SymbolicResult<Expr> {
    if matrix.rows() != matrix.cols() {
        return Err(SymbolicError::InvalidOperation(
            "Characteristic polynomial only defined for square matrices".to_string(),
        ));
    }

    let n = matrix.rows();

    // Create (A - λI)
    let mut a_minus_lambda_i = matrix.clone();
    for i in 0..n {
        let diagonal_elem = a_minus_lambda_i.get(i, i).unwrap().clone();
        let new_elem = Expr::add(diagonal_elem, Expr::mul(Expr::num(-1), Expr::sym("λ")));
        a_minus_lambda_i.set(i, i, new_elem)?;
    }

    // Compute det(A - λI)
    a_minus_lambda_i.determinant()
}

/// Compute eigenvalues for 2x2 matrices symbolically
/// Returns solutions to the characteristic equation
pub fn eigenvalues_2x2(matrix: &SymbolicMatrix) -> SymbolicResult<Vec<Expr>> {
    if matrix.rows() != 2 || matrix.cols() != 2 {
        return Err(SymbolicError::InvalidOperation(
            "This function only works for 2x2 matrices".to_string(),
        ));
    }

    // For a 2x2 matrix [[a, b], [c, d]]
    // Characteristic equation: λ² - (a+d)λ + (ad-bc) = 0
    // Eigenvalues: λ = (a+d ± √((a+d)² - 4(ad-bc))) / 2

    let a = matrix.get(0, 0).unwrap();
    let b = matrix.get(0, 1).unwrap();
    let c = matrix.get(1, 0).unwrap();
    let d = matrix.get(1, 1).unwrap();

    // Trace: a + d
    let trace = Expr::add(a.clone(), d.clone());

    // Determinant: ad - bc
    let det = Expr::add(
        Expr::mul(a.clone(), d.clone()),
        Expr::mul(Expr::num(-1), Expr::mul(b.clone(), c.clone())),
    );

    // Discriminant: (a+d)² - 4(ad-bc) = trace² - 4*det
    let discriminant = Expr::add(
        Expr::pow(trace.clone(), Expr::num(2)),
        Expr::mul(Expr::num(-4), det),
    );

    // λ₁ = (trace + √discriminant) / 2
    let lambda1 = Expr::mul(
        Expr::rational_unchecked(1, 2),
        Expr::add(
            trace.clone(),
            Expr::func("sqrt", vec![discriminant.clone()]),
        ),
    );

    // λ₂ = (trace - √discriminant) / 2
    let lambda2 = Expr::mul(
        Expr::rational_unchecked(1, 2),
        Expr::add(
            trace,
            Expr::mul(Expr::num(-1), Expr::func("sqrt", vec![discriminant])),
        ),
    );

    Ok(vec![lambda1, lambda2])
}

/// Compute the adjugate (classical adjoint) matrix
/// Used for matrix inversion: A⁻¹ = adj(A) / det(A)
pub fn adjugate(matrix: &SymbolicMatrix) -> SymbolicResult<SymbolicMatrix> {
    if matrix.rows() != matrix.cols() {
        return Err(SymbolicError::InvalidOperation(
            "Adjugate only defined for square matrices".to_string(),
        ));
    }

    let n = matrix.rows();
    let mut adj_data = vec![vec![Expr::num(0); n]; n];

    for i in 0..n {
        for j in 0..n {
            // Cofactor C_ij = (-1)^(i+j) * M_ij
            let minor = matrix.minor(i, j)?;
            let minor_det = minor.determinant()?;

            let sign = if (i + j) % 2 == 0 { 1 } else { -1 };
            let cofactor = Expr::mul(Expr::num(sign), minor_det);

            // Adjugate is transpose of cofactor matrix
            adj_data[j][i] = cofactor;
        }
    }

    SymbolicMatrix::new(adj_data)
}

/// Compute the inverse of a matrix symbolically
/// A⁻¹ = adj(A) / det(A)
pub fn matrix_inverse(matrix: &SymbolicMatrix) -> SymbolicResult<SymbolicMatrix> {
    if matrix.rows() != matrix.cols() {
        return Err(SymbolicError::InvalidOperation(
            "Inverse only defined for square matrices".to_string(),
        ));
    }

    let det = matrix.determinant()?;

    // Check if determinant is symbolically zero
    if det.is_zero() {
        return Err(SymbolicError::InvalidOperation(
            "Matrix is singular (determinant is zero)".to_string(),
        ));
    }

    let adj = adjugate(matrix)?;

    // Divide each element by determinant
    let inv_det = Expr::pow(det, Expr::num(-1));
    Ok(adj.scalar_mul(&inv_det))
}

impl SymbolicMatrix {
    /// Helper for eigenvalue computation - compute minor
    pub(crate) fn minor(&self, row: usize, col: usize) -> SymbolicResult<Self> {
        if self.rows() != self.cols() {
            return Err(SymbolicError::InvalidOperation(
                "Minor only defined for square matrices".to_string(),
            ));
        }

        let n = self.rows() - 1;
        let mut data = vec![vec![Expr::num(0); n]; n];

        let mut mi = 0;
        for i in 0..self.rows() {
            if i == row {
                continue;
            }
            let mut mj = 0;
            for j in 0..self.cols() {
                if j == col {
                    continue;
                }
                data[mi][mj] = self.get(i, j).unwrap().clone();
                mj += 1;
            }
            mi += 1;
        }

        Ok(SymbolicMatrix::new(data)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_characteristic_polynomial() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let char_poly = characteristic_polynomial(&mat).unwrap();
        println!("Characteristic polynomial: {}", char_poly);
    }

    #[test]
    fn test_eigenvalues_2x2_numeric() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(3), Expr::num(1)],
            vec![Expr::num(1), Expr::num(3)],
        ])
        .unwrap();

        let eigenvalues = eigenvalues_2x2(&mat).unwrap();
        assert_eq!(eigenvalues.len(), 2);
        println!("Eigenvalue 1: {}", eigenvalues[0]);
        println!("Eigenvalue 2: {}", eigenvalues[1]);
    }

    #[test]
    fn test_eigenvalues_2x2_symbolic() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("b"), Expr::sym("a")],
        ])
        .unwrap();

        let eigenvalues = eigenvalues_2x2(&mat).unwrap();
        println!("Symbolic eigenvalue 1: {}", eigenvalues[0]);
        println!("Symbolic eigenvalue 2: {}", eigenvalues[1]);
    }

    #[test]
    fn test_matrix_inverse() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ])
        .unwrap();

        let inv = matrix_inverse(&mat).unwrap();
        println!("Inverse matrix:\n{}", inv);
    }

    #[test]
    fn test_symbolic_inverse() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("c"), Expr::sym("d")],
        ])
        .unwrap();

        let inv = matrix_inverse(&mat).unwrap();
        println!("Symbolic inverse:\n{}", inv);
    }
}
