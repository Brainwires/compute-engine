//! Extensive-Form Games
//!
//! Games with sequential moves and imperfect information

use serde::{Deserialize, Serialize};

/// Game tree node
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GameNode {
    /// Player to move (None if terminal)
    pub player: Option<usize>,
    /// Available actions
    pub actions: Vec<String>,
    /// Child nodes for each action
    pub children: Vec<usize>,
    /// Payoffs (if terminal node)
    pub payoffs: Option<Vec<f64>>,
}

/// Extensive-form game
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtensiveGame {
    /// Game tree nodes
    pub nodes: Vec<GameNode>,
    /// Root node index
    pub root: usize,
}

impl ExtensiveGame {
    /// Create new extensive game
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: 0,
        }
    }

    /// Solve using backward induction (for perfect information games)
    pub fn backward_induction(&self) -> Vec<f64> {
        // Placeholder: would implement recursive backward induction
        vec![0.0; 2]
    }
}

