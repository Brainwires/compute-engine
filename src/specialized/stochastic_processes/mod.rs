//! Stochastic processes module
//!
//! Provides stochastic simulation capabilities:
//! - Random walks and Brownian motion
//! - Monte Carlo simulations
//! - Markov chains
//! - Stochastic differential equations

mod lib;
pub use lib::{
    BrownianMotionParams, generate_brownian_motion,
    MarkovChainParams, simulate_markov_chain,
    StochasticIntegralParams, compute_stochastic_integral_monte_carlo,
    generate_ornstein_uhlenbeck_process,
    simulate_poisson_process,
};
