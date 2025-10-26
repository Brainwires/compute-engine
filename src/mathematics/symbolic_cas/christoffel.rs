//! Christoffel symbols and tensor calculus operations
//!
//! Provides symbolic computation of connection coefficients, curvature tensors,
//! and other differential geometry operations for general relativity.

use super::differentiate::differentiate as diff;
use super::simplify::simplify;
use super::symbolic_eigenvalues::matrix_inverse;
use super::{Expr, IndexType, SymbolicError, SymbolicMatrix, SymbolicResult, SymbolicTensor};

/// Compute Christoffel symbols of the second kind from a metric tensor
///
/// Γ^μ_νλ = (1/2) g^μσ (∂_ν g_σλ + ∂_λ g_σν - ∂_σ g_νλ)
///
/// # Arguments
/// * `metric` - The metric tensor g_μν (covariant)
/// * `coords` - Coordinate variable names (e.g., ["t", "r", "θ", "φ"])
///
/// # Returns
/// A rank-3 tensor with index structure (contravariant, covariant, covariant)
pub fn christoffel_symbols(
    metric: &SymbolicMatrix,
    coords: &[&str],
) -> SymbolicResult<SymbolicTensor> {
    let n = metric.rows();

    if coords.len() != n {
        return Err(SymbolicError::InvalidOperation(format!(
            "Number of coordinates ({}) must match metric dimension ({})",
            coords.len(),
            n
        )));
    }

    // Compute inverse metric g^μν
    let inverse_metric = matrix_inverse(metric)?;

    // Compute partial derivatives of metric components
    let mut metric_derivatives = vec![vec![vec![Expr::num(0); n]; n]; n];

    for i in 0..n {
        for j in 0..n {
            let g_ij = metric.get(i, j).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid metric indices".to_string())
            })?;

            for k in 0..n {
                // ∂_k g_ij
                metric_derivatives[k][i][j] = diff(g_ij, coords[k]);
            }
        }
    }

    // Compute Christoffel symbols: Γ^μ_νλ = (1/2) g^μσ (∂_ν g_σλ + ∂_λ g_σν - ∂_σ g_νλ)
    let mut christoffel_data = Vec::new();

    for mu in 0..n {
        for nu in 0..n {
            for lambda in 0..n {
                let mut sum = Expr::num(0);

                for sigma in 0..n {
                    let g_inv = inverse_metric.get(mu, sigma).ok_or_else(|| {
                        SymbolicError::InvalidOperation(
                            "Invalid inverse metric indices".to_string(),
                        )
                    })?;

                    // ∂_ν g_σλ
                    let term1 = metric_derivatives[nu][sigma][lambda].clone();

                    // ∂_λ g_σν
                    let term2 = metric_derivatives[lambda][sigma][nu].clone();

                    // ∂_σ g_νλ
                    let term3 = metric_derivatives[sigma][nu][lambda].clone();

                    // (∂_ν g_σλ + ∂_λ g_σν - ∂_σ g_νλ)
                    let bracket =
                        Expr::add(Expr::add(term1, term2), Expr::mul(Expr::num(-1), term3));

                    // g^μσ * bracket
                    let contribution = Expr::mul(g_inv.clone(), bracket);
                    sum = Expr::add(sum, contribution);
                }

                // Multiply by 1/2 and simplify
                let gamma = Expr::mul(Expr::rational_unchecked(1, 2), sum);
                let gamma = simplify(&gamma);
                christoffel_data.push(gamma);
            }
        }
    }

    // Create rank-3 tensor: Γ^μ_νλ (contravariant, covariant, covariant)
    SymbolicTensor::new(
        vec![n, n, n],
        vec![
            IndexType::Contravariant,
            IndexType::Covariant,
            IndexType::Covariant,
        ],
        christoffel_data,
    )
}

/// Compute Christoffel symbols of the first kind
///
/// Γ_μνλ = (1/2) (∂_ν g_μλ + ∂_λ g_μν - ∂_μ g_νλ)
pub fn christoffel_first_kind(
    metric: &SymbolicMatrix,
    coords: &[&str],
) -> SymbolicResult<SymbolicTensor> {
    let n = metric.rows();

    if coords.len() != n {
        return Err(SymbolicError::InvalidOperation(format!(
            "Number of coordinates ({}) must match metric dimension ({})",
            coords.len(),
            n
        )));
    }

    // Compute partial derivatives of metric components
    let mut metric_derivatives = vec![vec![vec![Expr::num(0); n]; n]; n];

    for i in 0..n {
        for j in 0..n {
            let g_ij = metric.get(i, j).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid metric indices".to_string())
            })?;

            for k in 0..n {
                metric_derivatives[k][i][j] = diff(g_ij, coords[k]);
            }
        }
    }

    // Compute Γ_μνλ = (1/2) (∂_ν g_μλ + ∂_λ g_μν - ∂_μ g_νλ)
    let mut christoffel_data = Vec::new();

    for mu in 0..n {
        for nu in 0..n {
            for lambda in 0..n {
                let term1 = metric_derivatives[nu][mu][lambda].clone();
                let term2 = metric_derivatives[lambda][mu][nu].clone();
                let term3 = metric_derivatives[mu][nu][lambda].clone();

                let gamma = Expr::mul(
                    Expr::rational_unchecked(1, 2),
                    Expr::add(Expr::add(term1, term2), Expr::mul(Expr::num(-1), term3)),
                );
                let gamma = simplify(&gamma);

                christoffel_data.push(gamma);
            }
        }
    }

    // All indices are covariant
    SymbolicTensor::new(
        vec![n, n, n],
        vec![
            IndexType::Covariant,
            IndexType::Covariant,
            IndexType::Covariant,
        ],
        christoffel_data,
    )
}

/// Compute the geodesic equation coefficients from a metric
///
/// The geodesic equation is: d²x^μ/dτ² + Γ^μ_νλ (dx^ν/dτ)(dx^λ/dτ) = 0
///
/// Returns the Christoffel symbols which are the coefficients in this equation
pub fn geodesic_coefficients(
    metric: &SymbolicMatrix,
    coords: &[&str],
) -> SymbolicResult<SymbolicTensor> {
    christoffel_symbols(metric, coords)
}

/// Compute the Riemann curvature tensor
///
/// R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
///
/// # Arguments
/// * `metric` - The metric tensor g_μν
/// * `coords` - Coordinate variable names
///
/// # Returns
/// A rank-4 tensor with index structure (contravariant, covariant, covariant, covariant)
pub fn riemann_tensor(metric: &SymbolicMatrix, coords: &[&str]) -> SymbolicResult<SymbolicTensor> {
    let n = metric.rows();

    // First compute Christoffel symbols
    let christoffel = christoffel_symbols(metric, coords)?;

    // Compute partial derivatives of Christoffel symbols
    // ∂_μ Γ^ρ_νσ
    let mut christoffel_derivatives = vec![vec![vec![vec![Expr::num(0); n]; n]; n]; n];

    for rho in 0..n {
        for nu in 0..n {
            for sigma in 0..n {
                let gamma = christoffel.get(&[rho, nu, sigma])?;

                for mu in 0..n {
                    christoffel_derivatives[mu][rho][nu][sigma] = diff(gamma, coords[mu]);
                }
            }
        }
    }

    // Compute Riemann tensor components
    let mut riemann_data = Vec::new();

    for rho in 0..n {
        for sigma in 0..n {
            for mu in 0..n {
                for nu in 0..n {
                    // ∂_μ Γ^ρ_νσ
                    let term1 = christoffel_derivatives[mu][rho][nu][sigma].clone();

                    // ∂_ν Γ^ρ_μσ
                    let term2 = christoffel_derivatives[nu][rho][mu][sigma].clone();

                    // Γ^ρ_μλ Γ^λ_νσ (sum over λ)
                    let mut term3 = Expr::num(0);
                    for lambda in 0..n {
                        let gamma1 = christoffel.get(&[rho, mu, lambda])?;
                        let gamma2 = christoffel.get(&[lambda, nu, sigma])?;
                        term3 = Expr::add(term3, Expr::mul(gamma1.clone(), gamma2.clone()));
                    }

                    // Γ^ρ_νλ Γ^λ_μσ (sum over λ)
                    let mut term4 = Expr::num(0);
                    for lambda in 0..n {
                        let gamma1 = christoffel.get(&[rho, nu, lambda])?;
                        let gamma2 = christoffel.get(&[lambda, mu, sigma])?;
                        term4 = Expr::add(term4, Expr::mul(gamma1.clone(), gamma2.clone()));
                    }

                    // R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
                    let riemann_component = Expr::add(
                        Expr::add(term1, Expr::mul(Expr::num(-1), term2)),
                        Expr::add(term3, Expr::mul(Expr::num(-1), term4)),
                    );

                    let riemann_component = simplify(&riemann_component);
                    riemann_data.push(riemann_component);
                }
            }
        }
    }

    // Create rank-4 tensor: R^ρ_σμν
    SymbolicTensor::new(
        vec![n, n, n, n],
        vec![
            IndexType::Contravariant,
            IndexType::Covariant,
            IndexType::Covariant,
            IndexType::Covariant,
        ],
        riemann_data,
    )
}

/// Compute the Ricci tensor by contracting the Riemann tensor
///
/// R_μν = R^ρ_μρν (sum over ρ)
///
/// # Arguments
/// * `metric` - The metric tensor g_μν
/// * `coords` - Coordinate variable names
///
/// # Returns
/// A rank-2 covariant tensor (Ricci tensor)
pub fn ricci_tensor(metric: &SymbolicMatrix, coords: &[&str]) -> SymbolicResult<SymbolicMatrix> {
    let n = metric.rows();

    // Compute Riemann tensor
    let riemann = riemann_tensor(metric, coords)?;

    // Contract: R_μν = R^ρ_μρν (sum over ρ)
    let mut ricci_data = vec![vec![Expr::num(0); n]; n];

    for mu in 0..n {
        for nu in 0..n {
            let mut sum = Expr::num(0);

            for rho in 0..n {
                // R^ρ_μρν
                let component = riemann.get(&[rho, mu, rho, nu])?;
                sum = Expr::add(sum, component.clone());
            }

            ricci_data[mu][nu] = simplify(&sum);
        }
    }

    SymbolicMatrix::new(ricci_data)
}

/// Compute the Ricci scalar (curvature scalar)
///
/// R = g^μν R_μν (sum over μ, ν)
///
/// # Arguments
/// * `metric` - The metric tensor g_μν
/// * `coords` - Coordinate variable names
///
/// # Returns
/// The Ricci scalar as a symbolic expression
pub fn ricci_scalar(metric: &SymbolicMatrix, coords: &[&str]) -> SymbolicResult<Expr> {
    let n = metric.rows();

    // Compute Ricci tensor
    let ricci = ricci_tensor(metric, coords)?;

    // Compute inverse metric
    let inv_metric = matrix_inverse(metric)?;

    // R = g^μν R_μν
    let mut scalar = Expr::num(0);

    for mu in 0..n {
        for nu in 0..n {
            let g_inv = inv_metric.get(mu, nu).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid metric indices".to_string())
            })?;
            let r_mn = ricci.get(mu, nu).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid Ricci indices".to_string())
            })?;

            scalar = Expr::add(scalar, Expr::mul(g_inv.clone(), r_mn.clone()));
        }
    }

    Ok(simplify(&scalar))
}

/// Compute the Einstein tensor
///
/// G_μν = R_μν - (1/2) g_μν R
///
/// # Arguments
/// * `metric` - The metric tensor g_μν
/// * `coords` - Coordinate variable names
///
/// # Returns
/// The Einstein tensor as a symbolic matrix
pub fn einstein_tensor(metric: &SymbolicMatrix, coords: &[&str]) -> SymbolicResult<SymbolicMatrix> {
    let n = metric.rows();

    // Compute Ricci tensor and scalar
    let ricci = ricci_tensor(metric, coords)?;
    let scalar = ricci_scalar(metric, coords)?;

    // G_μν = R_μν - (1/2) g_μν R
    let mut einstein_data = vec![vec![Expr::num(0); n]; n];

    for mu in 0..n {
        for nu in 0..n {
            let r_mn = ricci.get(mu, nu).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid Ricci indices".to_string())
            })?;
            let g_mn = metric.get(mu, nu).ok_or_else(|| {
                SymbolicError::InvalidOperation("Invalid metric indices".to_string())
            })?;

            let term1 = r_mn.clone();
            let term2 = Expr::mul(
                Expr::mul(Expr::rational_unchecked(1, 2), g_mn.clone()),
                scalar.clone(),
            );

            einstein_data[mu][nu] = simplify(&Expr::add(term1, Expr::mul(Expr::num(-1), term2)));
        }
    }

    SymbolicMatrix::new(einstein_data)
}

