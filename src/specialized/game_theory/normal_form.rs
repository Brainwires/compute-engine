//! Normal-Form (Strategic-Form) Games
//!
//! Games represented as payoff matrices where players choose strategies simultaneously.
//! Includes Nash equilibrium finding and dominant strategy analysis.

use serde::{Deserialize, Serialize};
use super::{Player, StrategyProfile, Payoff};

/// Normal-form game representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NormalFormGame {
    /// Number of players
    pub num_players: usize,
    /// Number of strategies for each player
    pub num_strategies: Vec<usize>,
    /// Payoff function: maps strategy profile to payoff vector
    /// Stored as flattened multi-dimensional array
    pub payoffs: Vec<Payoff>,
}

impl NormalFormGame {
    /// Create a new normal-form game
    pub fn new(num_players: usize, num_strategies: Vec<usize>) -> Self {
        assert_eq!(num_players, num_strategies.len());

        let total_profiles: usize = num_strategies.iter().product();
        let payoffs = vec![vec![0.0; num_players]; total_profiles];

        Self {
            num_players,
            num_strategies,
            payoffs,
        }
    }

    /// Create a 2-player game from payoff matrices
    pub fn from_matrices(player1_payoffs: Vec<Vec<f64>>, player2_payoffs: Vec<Vec<f64>>) -> Self {
        let rows = player1_payoffs.len();
        let cols = player1_payoffs[0].len();

        let mut game = Self::new(2, vec![rows, cols]);

        for i in 0..rows {
            for j in 0..cols {
                let profile_idx = i * cols + j;
                game.payoffs[profile_idx] = vec![player1_payoffs[i][j], player2_payoffs[i][j]];
            }
        }

        game
    }

    /// Get payoff for a strategy profile
    pub fn get_payoff(&self, profile: &StrategyProfile) -> &Payoff {
        let idx = self.profile_to_index(profile);
        &self.payoffs[idx]
    }

    /// Convert strategy profile to flat index
    fn profile_to_index(&self, profile: &StrategyProfile) -> usize {
        let mut idx = 0;
        let mut multiplier = 1;

        for (i, &strategy) in profile.iter().enumerate().rev() {
            idx += strategy * multiplier;
            multiplier *= self.num_strategies[i];
        }

        idx
    }

    /// Find all Nash equilibria (pure strategy)
    pub fn find_nash_equilibria(&self) -> Vec<StrategyProfile> {
        let mut equilibria = Vec::new();

        // Generate all possible strategy profiles
        let profiles = self.generate_all_profiles();

        for profile in profiles {
            if self.is_nash_equilibrium(&profile) {
                equilibria.push(profile);
            }
        }

        equilibria
    }

    /// Check if a strategy profile is a Nash equilibrium
    pub fn is_nash_equilibrium(&self, profile: &StrategyProfile) -> bool {
        let current_payoffs = self.get_payoff(profile);

        // For each player
        for player in 0..self.num_players {
            let current_payoff = current_payoffs[player];

            // Check all possible deviations
            for alt_strategy in 0..self.num_strategies[player] {
                if alt_strategy == profile[player] {
                    continue;
                }

                // Create deviated profile
                let mut deviated_profile = profile.clone();
                deviated_profile[player] = alt_strategy;

                let deviated_payoff = self.get_payoff(&deviated_profile)[player];

                // If deviation is profitable, not a Nash equilibrium
                if deviated_payoff > current_payoff + 1e-10 {
                    return false;
                }
            }
        }

        true
    }

    /// Find dominant strategies for each player
    pub fn find_dominant_strategies(&self) -> Vec<Option<usize>> {
        let mut dominant = vec![None; self.num_players];

        for player in 0..self.num_players {
            if let Some(strategy) = self.find_dominant_strategy(player) {
                dominant[player] = Some(strategy);
            }
        }

        dominant
    }

    /// Find dominant strategy for a specific player
    fn find_dominant_strategy(&self, player: Player) -> Option<usize> {
        'strategy_loop: for candidate in 0..self.num_strategies[player] {
            // Check if candidate dominates all other strategies
            for other in 0..self.num_strategies[player] {
                if candidate == other {
                    continue;
                }

                // Check if candidate dominates other for all opponent strategies
                if !self.strategy_dominates(player, candidate, other) {
                    continue 'strategy_loop;
                }
            }

            return Some(candidate);
        }

        None
    }

    /// Check if strategy1 dominates strategy2 for a player
    fn strategy_dominates(&self, player: Player, strategy1: usize, strategy2: usize) -> bool {
        // Generate all opponent strategy combinations
        let profiles = self.generate_all_profiles();

        for profile in profiles {
            let mut profile1 = profile.clone();
            profile1[player] = strategy1;
            let payoff1 = self.get_payoff(&profile1)[player];

            let mut profile2 = profile.clone();
            profile2[player] = strategy2;
            let payoff2 = self.get_payoff(&profile2)[player];

            if payoff1 <= payoff2 {
                return false;
            }
        }

        true
    }

    /// Generate all possible strategy profiles
    fn generate_all_profiles(&self) -> Vec<StrategyProfile> {
        let mut profiles = Vec::new();
        let total: usize = self.num_strategies.iter().product();

        for i in 0..total {
            profiles.push(self.index_to_profile(i));
        }

        profiles
    }

    /// Convert flat index to strategy profile
    fn index_to_profile(&self, mut idx: usize) -> StrategyProfile {
        let mut profile = vec![0; self.num_players];

        for i in (0..self.num_players).rev() {
            profile[i] = idx % self.num_strategies[i];
            idx /= self.num_strategies[i];
        }

        profile
    }

    /// Check if outcome is Pareto optimal
    pub fn is_pareto_optimal(&self, profile: &StrategyProfile) -> bool {
        let current_payoffs = self.get_payoff(profile);

        // Check all other profiles
        for other_profile in self.generate_all_profiles() {
            if &other_profile == profile {
                continue;
            }

            let other_payoffs = self.get_payoff(&other_profile);

            // Check if other profile Pareto dominates current
            let mut all_better_or_equal = true;
            let mut at_least_one_better = false;

            for player in 0..self.num_players {
                if other_payoffs[player] < current_payoffs[player] - 1e-10 {
                    all_better_or_equal = false;
                    break;
                }
                if other_payoffs[player] > current_payoffs[player] + 1e-10 {
                    at_least_one_better = true;
                }
            }

            if all_better_or_equal && at_least_one_better {
                return false;
            }
        }

        true
    }
}

/// Classic game: Prisoner's Dilemma
pub fn prisoners_dilemma() -> NormalFormGame {
    // Cooperate = 0, Defect = 1
    // Payoffs: (C,C)=(-1,-1), (C,D)=(-3,0), (D,C)=(0,-3), (D,D)=(-2,-2)
    NormalFormGame::from_matrices(
        vec![
            vec![-1.0, -3.0],  // Player 1: Cooperate
            vec![0.0, -2.0],   // Player 1: Defect
        ],
        vec![
            vec![-1.0, 0.0],   // Player 2: Cooperate
            vec![-3.0, -2.0],  // Player 2: Defect
        ],
    )
}

/// Classic game: Battle of the Sexes
pub fn battle_of_sexes() -> NormalFormGame {
    // Opera = 0, Football = 1
    NormalFormGame::from_matrices(
        vec![
            vec![2.0, 0.0],  // Player 1: Opera
            vec![0.0, 1.0],  // Player 1: Football
        ],
        vec![
            vec![1.0, 0.0],  // Player 2: Opera
            vec![0.0, 2.0],  // Player 2: Football
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prisoners_dilemma_nash() {
        let game = prisoners_dilemma();
        let equilibria = game.find_nash_equilibria();

        // (Defect, Defect) should be the unique Nash equilibrium
        assert_eq!(equilibria.len(), 1);
        assert_eq!(equilibria[0], vec![1, 1]);
    }

    #[test]
    fn test_prisoners_dilemma_dominant() {
        let game = prisoners_dilemma();
        let dominant = game.find_dominant_strategies();

        // Both players have Defect as dominant strategy
        assert_eq!(dominant[0], Some(1));
        assert_eq!(dominant[1], Some(1));
    }

    #[test]
    fn test_battle_of_sexes_nash() {
        let game = battle_of_sexes();
        let equilibria = game.find_nash_equilibria();

        // Two pure Nash equilibria: (Opera, Opera) and (Football, Football)
        assert_eq!(equilibria.len(), 2);
    }

    #[test]
    fn test_pareto_optimal() {
        let game = prisoners_dilemma();

        // (Cooperate, Cooperate) is Pareto optimal
        assert!(game.is_pareto_optimal(&vec![0, 0]));

        // (Defect, Defect) is not Pareto optimal
        assert!(!game.is_pareto_optimal(&vec![1, 1]));
    }

    #[test]
    fn test_game_creation() {
        let game = NormalFormGame::new(2, vec![2, 3]);
        assert_eq!(game.num_players, 2);
        assert_eq!(game.payoffs.len(), 6); // 2 * 3 profiles
    }
}
