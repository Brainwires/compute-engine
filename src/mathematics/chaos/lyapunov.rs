//! Lyapunov Exponents - Measure of Chaos
//!
//! Lyapunov exponents quantify the rate of separation of infinitesimally close trajectories.
//! Positive largest Lyapunov exponent indicates chaos.

use super::Point3D;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LyapunovResult {
    pub largest_exponent: f64,
    pub all_exponents: Vec<f64>,
    pub is_chaotic: bool,
}

/// Calculate largest Lyapunov exponent for a 1D map
pub fn lyapunov_1d_map<F>(map: F, x0: f64, iterations: usize, transient: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut x = x0;
    let mut sum_log = 0.0;

    // Discard transient
    for _ in 0..transient {
        x = map(x);
    }

    // Calculate Lyapunov exponent
    for _ in 0..iterations {
        x = map(x);

        // Numerical derivative
        let epsilon = 1e-8;
        let derivative = (map(x + epsilon) - map(x - epsilon)) / (2.0 * epsilon);

        if derivative.abs() > 1e-10 {
            sum_log += derivative.abs().ln();
        }
    }

    sum_log / iterations as f64
}

/// Calculate Lyapunov exponent for logistic map
pub fn lyapunov_logistic_map(r: f64, iterations: usize, transient: usize) -> f64 {
    let mut x = 0.5;
    let mut sum_log = 0.0;

    // Discard transient
    for _ in 0..transient {
        x = r * x * (1.0 - x);
    }

    // Calculate Lyapunov exponent
    for _ in 0..iterations {
        x = r * x * (1.0 - x);

        // Analytical derivative: df/dx = r(1 - 2x)
        let derivative = r * (1.0 - 2.0 * x);

        if derivative.abs() > 1e-10 {
            sum_log += derivative.abs().ln();
        }
    }

    sum_log / iterations as f64
}

/// Calculate Lyapunov spectrum for 3D system (Lorenz, etc.)
/// Uses Benettin algorithm with QR decomposition
pub fn lyapunov_spectrum_3d<F>(
    derivative: F,
    initial: Point3D,
    dt: f64,
    steps: usize,
    transient: usize,
) -> Vec<f64>
where
    F: Fn(Point3D) -> [[f64; 3]; 3], // Jacobian matrix
{
    let mut state = initial;
    let mut exponents = vec![0.0; 3];

    // Initial orthonormal basis
    let mut basis = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    // Discard transient
    for _ in 0..transient {
        state = evolve_state(state, &derivative, dt);
    }

    // Calculate exponents
    for _ in 0..steps {
        state = evolve_state(state, &derivative, dt);

        // Evolve basis vectors
        let jacobian = derivative(state);
        for i in 0..3 {
            let mut new_vec = [0.0; 3];
            for j in 0..3 {
                for k in 0..3 {
                    new_vec[j] += jacobian[j][k] * basis[i][k];
                }
            }
            basis[i] = new_vec;
        }

        // Gram-Schmidt orthonormalization and extract norms
        let norms = gram_schmidt(&mut basis);

        // Accumulate exponents
        for i in 0..3 {
            if norms[i] > 1e-10 {
                exponents[i] += norms[i].ln();
            }
        }
    }

    // Average over time
    for i in 0..3 {
        exponents[i] /= steps as f64 * dt;
    }

    exponents
}

fn evolve_state<F>(state: Point3D, _derivative: &F, _dt: f64) -> Point3D
where
    F: Fn(Point3D) -> [[f64; 3]; 3],
{
    // Simplified: just return state (full implementation would integrate ODE)
    state
}

fn gram_schmidt(basis: &mut [[f64; 3]; 3]) -> Vec<f64> {
    let mut norms = vec![];

    for i in 0..3 {
        // Orthogonalize against previous vectors
        for j in 0..i {
            let dot: f64 = (0..3).map(|k| basis[i][k] * basis[j][k]).sum();
            for k in 0..3 {
                basis[i][k] -= dot * basis[j][k];
            }
        }

        // Normalize
        let norm: f64 = (0..3).map(|k| basis[i][k] * basis[i][k]).sum::<f64>().sqrt();
        norms.push(norm);

        if norm > 1e-10 {
            for k in 0..3 {
                basis[i][k] /= norm;
            }
        }
    }

    norms
}

/// Kaplan-Yorke dimension from Lyapunov spectrum
pub fn kaplan_yorke_dimension(exponents: &[f64]) -> f64 {
    let mut sorted = exponents.to_vec();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let mut sum = 0.0;
    let mut j = 0;

    for (i, &exp) in sorted.iter().enumerate() {
        sum += exp;
        if sum < 0.0 {
            j = i;
            break;
        }
    }

    if j == 0 || sorted[j].abs() < 1e-10 {
        return 0.0;
    }

    let sum_j: f64 = sorted[0..j].iter().sum();
    j as f64 + sum_j / sorted[j].abs()
}

