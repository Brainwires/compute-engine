//! Game Theory Module
//!
//! Strategic decision-making and equilibrium analysis:
//! - Normal-form games and Nash equilibrium
//! - Extensive-form games
//! - Cooperative games and Shapley value
//! - Auction theory
//! - Evolutionary game theory

pub mod normal_form;
pub mod extensive_form;
pub mod cooperative;
pub mod evolutionary;

pub use normal_form::*;
pub use extensive_form::*;
pub use cooperative::*;
pub use evolutionary::*;

use serde::{Deserialize, Serialize};

/// Player in a game
pub type Player = usize;

/// Strategy profile (one strategy choice per player)
pub type StrategyProfile = Vec<usize>;

/// Payoff vector (one payoff per player)
pub type Payoff = Vec<f64>;

/// Game solution concept
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum SolutionConcept {
    /// Nash equilibrium
    Nash,
    /// Dominant strategy equilibrium
    Dominant,
    /// Pareto optimal
    Pareto,
    /// Minimax strategy
    Minimax,
}
