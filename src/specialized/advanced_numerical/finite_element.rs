//! Finite Element Method (FEM)
//!
//! Basic finite element method for solving boundary value problems.

use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use super::BoundaryCondition;

/// 1D Finite Element solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FiniteElementSolver1D {
    /// Number of elements
    pub num_elements: usize,
    /// Domain bounds [a, b]
    pub domain: (f64, f64),
    /// Boundary conditions at x=a and x=b
    pub boundary_conditions: (BoundaryCondition, BoundaryCondition),
}

impl FiniteElementSolver1D {
    /// Create a new 1D FEM solver
    pub fn new(
        num_elements: usize,
        domain: (f64, f64),
        boundary_conditions: (BoundaryCondition, BoundaryCondition),
    ) -> Self {
        Self {
            num_elements,
            domain,
            boundary_conditions,
        }
    }

    /// Solve 1D boundary value problem: -d²u/dx² = f(x)
    /// with given boundary conditions
    pub fn solve(&self, source: fn(f64) -> f64) -> Vec<f64> {
        let n = self.num_elements + 1; // Number of nodes
        let (a, b) = self.domain;
        let h = (b - a) / self.num_elements as f64; // Element size

        // Create nodes
        let nodes: Vec<f64> = (0..n).map(|i| a + i as f64 * h).collect();

        // Assemble stiffness matrix and load vector
        let mut K = DMatrix::<f64>::zeros(n, n);
        let mut F = DVector::<f64>::zeros(n);

        // Assemble element contributions
        for i in 0..self.num_elements {
            let x_left = nodes[i];
            let x_right = nodes[i + 1];
            let x_mid = (x_left + x_right) / 2.0;

            // Element stiffness matrix (for linear elements)
            // K_e = (1/h) * [[1, -1], [-1, 1]]
            let k_local = 1.0 / h;

            K[(i, i)] += k_local;
            K[(i, i + 1)] += -k_local;
            K[(i + 1, i)] += -k_local;
            K[(i + 1, i + 1)] += k_local;

            // Element load vector (using midpoint rule)
            // F_e = f(x_mid) * h * [0.5, 0.5]
            let f_mid = source(x_mid);
            F[i] += f_mid * h / 2.0;
            F[i + 1] += f_mid * h / 2.0;
        }

        // Apply boundary conditions
        match &self.boundary_conditions.0 {
            BoundaryCondition::Dirichlet(value) => {
                // Set first row to identity and RHS to boundary value
                for j in 0..n {
                    K[(0, j)] = 0.0;
                }
                K[(0, 0)] = 1.0;
                F[0] = *value;
            }
            _ => {} // Neumann and other BCs would be implemented here
        }

        match &self.boundary_conditions.1 {
            BoundaryCondition::Dirichlet(value) => {
                // Set last row to identity and RHS to boundary value
                for j in 0..n {
                    K[(n - 1, j)] = 0.0;
                }
                K[(n - 1, n - 1)] = 1.0;
                F[n - 1] = *value;
            }
            _ => {} // Neumann and other BCs would be implemented here
        }

        // Solve linear system K * u = F
        match K.lu().solve(&F) {
            Some(solution) => solution.as_slice().to_vec(),
            None => vec![0.0; n], // Return zeros if singular
        }
    }

    /// Get node positions
    pub fn get_nodes(&self) -> Vec<f64> {
        let (a, b) = self.domain;
        let n = self.num_elements + 1;
        let h = (b - a) / self.num_elements as f64;
        (0..n).map(|i| a + i as f64 * h).collect()
    }
}

/// 2D Finite Element solver (simplified)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FiniteElementSolver2D {
    /// Number of elements in x-direction
    pub nx: usize,
    /// Number of elements in y-direction
    pub ny: usize,
    /// Domain bounds [[x_min, x_max], [y_min, y_max]]
    pub domain: ((f64, f64), (f64, f64)),
}

impl FiniteElementSolver2D {
    /// Create a new 2D FEM solver
    pub fn new(nx: usize, ny: usize, domain: ((f64, f64), (f64, f64))) -> Self {
        Self { nx, ny, domain }
    }

    /// Solve 2D Poisson equation: -∇²u = f(x,y)
    /// with Dirichlet boundary conditions u=0 on boundary
    pub fn solve_poisson(&self, source: fn(f64, f64) -> f64) -> Vec<Vec<f64>> {
        let n_nodes_x = self.nx + 1;
        let n_nodes_y = self.ny + 1;
        let n_total = n_nodes_x * n_nodes_y;

        let ((x_min, x_max), (y_min, y_max)) = self.domain;
        let hx = (x_max - x_min) / self.nx as f64;
        let hy = (y_max - y_min) / self.ny as f64;

        // Using finite difference approximation for simplicity
        // (Full FEM with triangular elements would be more complex)
        let mut K = DMatrix::<f64>::zeros(n_total, n_total);
        let mut F = DVector::<f64>::zeros(n_total);

        // Node indexing: node(i,j) = i + j * n_nodes_x
        for j in 0..n_nodes_y {
            for i in 0..n_nodes_x {
                let idx = i + j * n_nodes_x;
                let x = x_min + i as f64 * hx;
                let y = y_min + j as f64 * hy;

                // Boundary nodes
                if i == 0 || i == n_nodes_x - 1 || j == 0 || j == n_nodes_y - 1 {
                    K[(idx, idx)] = 1.0;
                    F[idx] = 0.0; // Dirichlet BC: u = 0
                } else {
                    // Interior nodes: 5-point stencil
                    // -u_xx - u_yy = f
                    K[(idx, idx)] = 2.0 / (hx * hx) + 2.0 / (hy * hy);
                    K[(idx, idx - 1)] = -1.0 / (hx * hx); // left
                    K[(idx, idx + 1)] = -1.0 / (hx * hx); // right
                    K[(idx, idx - n_nodes_x)] = -1.0 / (hy * hy); // down
                    K[(idx, idx + n_nodes_x)] = -1.0 / (hy * hy); // up

                    F[idx] = source(x, y);
                }
            }
        }

        // Solve linear system
        let solution = match K.lu().solve(&F) {
            Some(sol) => sol,
            None => DVector::zeros(n_total),
        };

        // Convert to 2D grid
        let mut result = vec![vec![0.0; n_nodes_x]; n_nodes_y];
        for j in 0..n_nodes_y {
            for i in 0..n_nodes_x {
                let idx = i + j * n_nodes_x;
                result[j][i] = solution[idx];
            }
        }

        result
    }
}

