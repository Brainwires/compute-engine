//! Symbolic matrix operations
//!
//! Matrices where elements can be symbolic expressions

use super::{Expr, SymbolicError, SymbolicResult};
use std::fmt;

/// A matrix with symbolic expression elements
#[derive(Debug, Clone, PartialEq)]
pub struct SymbolicMatrix {
    rows: usize,
    cols: usize,
    data: Vec<Vec<Expr>>,
}

impl SymbolicMatrix {
    /// Create a new symbolic matrix from a 2D vector of expressions
    pub fn new(data: Vec<Vec<Expr>>) -> SymbolicResult<Self> {
        if data.is_empty() {
            return Err(SymbolicError::InvalidOperation("Empty matrix".to_string()));
        }

        let rows = data.len();
        let cols = data[0].len();

        if cols == 0 {
            return Err(SymbolicError::InvalidOperation("Matrix has no columns".to_string()));
        }

        // Verify all rows have the same number of columns
        for (i, row) in data.iter().enumerate() {
            if row.len() != cols {
                return Err(SymbolicError::InvalidOperation(
                    format!("Row {} has {} columns, expected {}", i, row.len(), cols)
                ));
            }
        }

        Ok(SymbolicMatrix { rows, cols, data })
    }

    /// Create an identity matrix of size n x n
    pub fn identity(n: usize) -> Self {
        let mut data = vec![vec![Expr::num(0); n]; n];
        for i in 0..n {
            data[i][i] = Expr::num(1);
        }
        SymbolicMatrix { rows: n, cols: n, data }
    }

    /// Create a zero matrix of size rows x cols
    pub fn zeros(rows: usize, cols: usize) -> Self {
        let data = vec![vec![Expr::num(0); cols]; rows];
        SymbolicMatrix { rows, cols, data }
    }

    /// Get the number of rows
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Get the number of columns
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Get element at (row, col)
    pub fn get(&self, row: usize, col: usize) -> Option<&Expr> {
        self.data.get(row)?.get(col)
    }

    /// Set element at (row, col)
    pub fn set(&mut self, row: usize, col: usize, expr: Expr) -> SymbolicResult<()> {
        if row >= self.rows || col >= self.cols {
            return Err(SymbolicError::InvalidOperation("Index out of bounds".to_string()));
        }
        self.data[row][col] = expr;
        Ok(())
    }

    /// Transpose the matrix
    pub fn transpose(&self) -> Self {
        let mut data = vec![vec![Expr::num(0); self.rows]; self.cols];
        for i in 0..self.rows {
            for j in 0..self.cols {
                data[j][i] = self.data[i][j].clone();
            }
        }
        SymbolicMatrix {
            rows: self.cols,
            cols: self.rows,
            data,
        }
    }

    /// Matrix addition
    pub fn add(&self, other: &SymbolicMatrix) -> SymbolicResult<Self> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(SymbolicError::InvalidOperation(
                format!("Matrix dimensions don't match: {}x{} vs {}x{}",
                    self.rows, self.cols, other.rows, other.cols)
            ));
        }

        let mut result = vec![vec![Expr::num(0); self.cols]; self.rows];
        for i in 0..self.rows {
            for j in 0..self.cols {
                result[i][j] = Expr::add(self.data[i][j].clone(), other.data[i][j].clone());
            }
        }

        Ok(SymbolicMatrix {
            rows: self.rows,
            cols: self.cols,
            data: result,
        })
    }

    /// Matrix multiplication
    pub fn mul(&self, other: &SymbolicMatrix) -> SymbolicResult<Self> {
        if self.cols != other.rows {
            return Err(SymbolicError::InvalidOperation(
                format!("Matrix dimensions incompatible for multiplication: {}x{} × {}x{}",
                    self.rows, self.cols, other.rows, other.cols)
            ));
        }

        let mut result = vec![vec![Expr::num(0); other.cols]; self.rows];
        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut sum = Expr::num(0);
                for k in 0..self.cols {
                    let term = Expr::mul(
                        self.data[i][k].clone(),
                        other.data[k][j].clone()
                    );
                    sum = Expr::add(sum, term);
                }
                result[i][j] = sum;
            }
        }

        Ok(SymbolicMatrix {
            rows: self.rows,
            cols: other.cols,
            data: result,
        })
    }

    /// Scalar multiplication
    pub fn scalar_mul(&self, scalar: &Expr) -> Self {
        let mut result = vec![vec![Expr::num(0); self.cols]; self.rows];
        for i in 0..self.rows {
            for j in 0..self.cols {
                result[i][j] = Expr::mul(scalar.clone(), self.data[i][j].clone());
            }
        }
        SymbolicMatrix {
            rows: self.rows,
            cols: self.cols,
            data: result,
        }
    }

    /// Compute determinant (for square matrices)
    pub fn determinant(&self) -> SymbolicResult<Expr> {
        if self.rows != self.cols {
            return Err(SymbolicError::InvalidOperation(
                "Determinant only defined for square matrices".to_string()
            ));
        }

        match self.rows {
            0 => Ok(Expr::num(1)),
            1 => Ok(self.data[0][0].clone()),
            2 => {
                // det = a*d - b*c
                let ad = Expr::mul(self.data[0][0].clone(), self.data[1][1].clone());
                let bc = Expr::mul(self.data[0][1].clone(), self.data[1][0].clone());
                Ok(Expr::add(ad, Expr::mul(Expr::num(-1), bc)))
            }
            3 => {
                // Sarrus rule for 3x3
                let mut sum = Expr::num(0);

                // Positive terms
                for i in 0..3 {
                    let mut prod = self.data[0][i].clone();
                    prod = Expr::mul(prod, self.data[1][(i + 1) % 3].clone());
                    prod = Expr::mul(prod, self.data[2][(i + 2) % 3].clone());
                    sum = Expr::add(sum, prod);
                }

                // Negative terms
                for i in 0..3 {
                    let mut prod = self.data[0][i].clone();
                    prod = Expr::mul(prod, self.data[1][(i + 2) % 3].clone());
                    prod = Expr::mul(prod, self.data[2][(i + 1) % 3].clone());
                    sum = Expr::add(sum, Expr::mul(Expr::num(-1), prod));
                }

                Ok(sum)
            }
            _ => {
                // Laplace expansion for larger matrices
                let mut det = Expr::num(0);
                for j in 0..self.cols {
                    let minor = self.minor(0, j)?;
                    let cofactor = minor.determinant()?;

                    let sign = if j % 2 == 0 { 1 } else { -1 };
                    let term = Expr::mul(
                        Expr::num(sign),
                        Expr::mul(self.data[0][j].clone(), cofactor)
                    );
                    det = Expr::add(det, term);
                }
                Ok(det)
            }
        }
    }

    /// Get trace (sum of diagonal elements)
    pub fn trace(&self) -> SymbolicResult<Expr> {
        if self.rows != self.cols {
            return Err(SymbolicError::InvalidOperation(
                "Trace only defined for square matrices".to_string()
            ));
        }

        let mut sum = Expr::num(0);
        for i in 0..self.rows {
            sum = Expr::add(sum, self.data[i][i].clone());
        }
        Ok(sum)
    }

    // Note: minor() method is defined in symbolic_eigenvalues.rs as pub(crate)
}

impl fmt::Display for SymbolicMatrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "[")?;
        for row in &self.data {
            write!(f, "  [")?;
            for (j, elem) in row.iter().enumerate() {
                if j > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", elem)?;
            }
            writeln!(f, "]")?;
        }
        write!(f, "]")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mathematics::symbolic_cas::Expr;

    #[test]
    fn test_identity_matrix() {
        let identity = SymbolicMatrix::identity(3);
        assert_eq!(identity.rows(), 3);
        assert_eq!(identity.cols(), 3);
        assert_eq!(identity.get(0, 0), Some(&Expr::num(1)));
        assert_eq!(identity.get(0, 1), Some(&Expr::num(0)));
    }

    #[test]
    fn test_transpose() {
        let data = vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ];
        let mat = SymbolicMatrix::new(data).unwrap();
        let transposed = mat.transpose();

        assert_eq!(transposed.rows(), 2);
        assert_eq!(transposed.cols(), 2);
        assert_eq!(transposed.get(0, 1), Some(&Expr::num(3)));
    }

    #[test]
    fn test_matrix_mul() {
        let a = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ]).unwrap();

        let b = SymbolicMatrix::new(vec![
            vec![Expr::num(5), Expr::num(6)],
            vec![Expr::num(7), Expr::num(8)],
        ]).unwrap();

        let c = a.mul(&b).unwrap();
        assert_eq!(c.rows(), 2);
        assert_eq!(c.cols(), 2);
    }

    #[test]
    fn test_determinant_2x2() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::num(1), Expr::num(2)],
            vec![Expr::num(3), Expr::num(4)],
        ]).unwrap();

        let det = mat.determinant().unwrap();
        // det = 1*4 - 2*3 = -2
        println!("Determinant: {}", det);
    }

    #[test]
    fn test_symbolic_determinant() {
        let mat = SymbolicMatrix::new(vec![
            vec![Expr::sym("a"), Expr::sym("b")],
            vec![Expr::sym("c"), Expr::sym("d")],
        ]).unwrap();

        let det = mat.determinant().unwrap();
        // det = a*d - b*c
        println!("Symbolic determinant: {}", det);
    }
}
