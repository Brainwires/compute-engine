//! Stochastic processes module
//!
//! Provides stochastic simulation capabilities:
//! - Random walks and Brownian motion
//! - Monte Carlo simulations
//! - Markov chains
//! - Stochastic differential equations

mod lib;
pub use lib::{
    BrownianMotionParams, MarkovChainParams, StochasticIntegralParams,
    compute_stochastic_integral_monte_carlo, generate_brownian_motion,
    generate_ornstein_uhlenbeck_process, simulate_markov_chain, simulate_poisson_process,
};
