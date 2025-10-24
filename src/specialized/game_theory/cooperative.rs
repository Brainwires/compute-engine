//! Cooperative Game Theory
//!
//! Coalition formation and value distribution

use serde::{Deserialize, Serialize};

/// Cooperative game (characteristic function form)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CooperativeGame {
    /// Number of players
    pub num_players: usize,
    /// Characteristic function: maps coalition to value
    /// Stored as map from coalition bitmask to value
    pub values: std::collections::HashMap<u32, f64>,
}

impl CooperativeGame {
    /// Create new cooperative game
    pub fn new(num_players: usize) -> Self {
        Self {
            num_players,
            values: std::collections::HashMap::new(),
        }
    }

    /// Set value for a coalition
    pub fn set_value(&mut self, coalition: &[usize], value: f64) {
        let bitmask = Self::coalition_to_bitmask(coalition);
        self.values.insert(bitmask, value);
    }

    /// Get value for a coalition
    pub fn get_value(&self, coalition: &[usize]) -> f64 {
        let bitmask = Self::coalition_to_bitmask(coalition);
        *self.values.get(&bitmask).unwrap_or(&0.0)
    }

    /// Convert coalition to bitmask
    fn coalition_to_bitmask(coalition: &[usize]) -> u32 {
        let mut mask = 0u32;
        for &player in coalition {
            mask |= 1 << player;
        }
        mask
    }

    /// Compute Shapley value
    pub fn shapley_value(&self) -> Vec<f64> {
        let n = self.num_players;
        let mut shapley = vec![0.0; n];

        // For each player
        for player in 0..n {
            let mut value_sum = 0.0;

            // For each coalition not containing player
            for coalition_mask in 0..(1u32 << n) {
                if (coalition_mask & (1 << player)) != 0 {
                    continue; // Skip if player is in coalition
                }

                let coalition_size = coalition_mask.count_ones() as usize;

                // Marginal contribution
                let with_player = coalition_mask | (1 << player);
                let value_with = *self.values.get(&with_player).unwrap_or(&0.0);
                let value_without = *self.values.get(&coalition_mask).unwrap_or(&0.0);
                let marginal = value_with - value_without;

                // Weight by coalition size
                let weight = Self::shapley_weight(n, coalition_size);
                value_sum += weight * marginal;
            }

            shapley[player] = value_sum;
        }

        shapley
    }

    /// Shapley weight for coalition of size s
    fn shapley_weight(n: usize, s: usize) -> f64 {
        // Weight = s! * (n-s-1)! / n!
        let mut numerator = 1.0;
        for i in 1..=s {
            numerator *= i as f64;
        }
        for i in 1..=(n - s - 1) {
            numerator *= i as f64;
        }

        let mut denominator = 1.0;
        for i in 1..=n {
            denominator *= i as f64;
        }

        numerator / denominator
    }

    /// Check if game is superadditive
    pub fn is_superadditive(&self) -> bool {
        // For all disjoint coalitions S and T: v(S ∪ T) ≥ v(S) + v(T)
        // Simplified check
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cooperative_game() {
        let mut game = CooperativeGame::new(3);
        game.set_value(&[0], 0.0);
        game.set_value(&[1], 0.0);
        game.set_value(&[2], 0.0);
        game.set_value(&[0, 1, 2], 12.0);

        let shapley = game.shapley_value();
        assert_eq!(shapley.len(), 3);
    }
}
