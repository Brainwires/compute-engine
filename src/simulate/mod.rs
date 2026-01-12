//! Unified Simulator Module
//!
//! Routes simulation requests to ODE/PDE solvers, stochastic process generators,
//! fluid dynamics solvers, and financial model simulators.
//!
//! This module provides the SIMULATE tool implementation for:
//! - Time evolution (Euler, RK4, adaptive step, implicit Euler)
//! - Stochastic processes (Brownian, geometric Brownian, OU, Poisson, Lévy, etc.)
//! - Fluid dynamics (Lattice Boltzmann, Quantum NS 1D/2D, NS 3D)
//! - Financial models (Black-Scholes, Heston, SABR, stochastic volatility)

pub mod finance;
pub mod fluids;
pub mod ode;
pub mod stochastic;

use crate::engine::*;

pub struct UnifiedSimulator;

impl UnifiedSimulator {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UnifiedSimulator {
    fn default() -> Self {
        Self::new()
    }
}

impl Simulate for UnifiedSimulator {
    fn simulate(&self, input: &SimulateInput) -> ToolResult<SimulateOutput> {
        match &input.model {
            SimulationModel::TimeEvolution(method) => ode::simulate_time_evolution(method, input),
            SimulationModel::Stochastic(process) => stochastic::simulate_stochastic(process, input),
            SimulationModel::FluidDynamics(fluid_sim) => fluids::simulate_fluid(fluid_sim, input),
            SimulationModel::Finance(finance_model) => finance::simulate_finance(finance_model, input),
        }
    }
}
