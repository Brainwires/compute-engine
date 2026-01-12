//! Unified Solver Module
//!
//! Routes solve requests to appropriate domain modules based on equation type.
//!
//! This module provides the SOLVE tool implementation for:
//! - Einstein field equations (general relativity)
//! - Fluid dynamics (Navier-Stokes, cavity flow)
//! - Electromagnetic equations (Maxwell, wave, Helmholtz)
//! - Chemical equations (balancing, thermodynamics, kinetics)
//! - Differential equations (ODE, PDE, boundary value problems)
//! - Linear systems
//! - Number theory (primality, factorization)
//! - Differential geometry (geodesics, parallel transport)
//! - Optimization (curve fitting, minimization)

pub mod differential;
pub mod equations;
pub mod optimization;
pub mod physics;
pub mod specialized;

use crate::engine::*;

pub struct UnifiedSolver;

impl UnifiedSolver {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UnifiedSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl Solve for UnifiedSolver {
    fn solve(&self, input: &SolveInput) -> ToolResult<SolveOutput> {
        match &input.equation_type {
            EquationType::Einstein(eq) => physics::solve_einstein(eq, input),
            EquationType::Fluid(eq) => physics::solve_fluid(eq, input),
            EquationType::Electromagnetic(em_eq) => physics::solve_electromagnetic(em_eq, input),
            EquationType::Chemical(chem_eq) => physics::solve_chemical(chem_eq, input),
            EquationType::Differential(diff_eq) => differential::solve_differential(diff_eq, input),
            EquationType::LinearSystem => equations::solve_linear_system(input),
            EquationType::RootFinding => equations::solve_root_finding(input),
            EquationType::NumberTheory(nt_prob) => specialized::solve_number_theory(nt_prob, input),
            EquationType::DifferentialGeometry(dg_prob) => {
                equations::solve_diff_geometry(dg_prob, input)
            }
            EquationType::Optimize(opt_method) => optimization::solve_optimization(opt_method, input),
        }
    }
}
