//! Evolutionary Game Theory
//!
//! Population dynamics and evolutionary stable strategies

use serde::{Deserialize, Serialize};

/// Evolutionary game with replicator dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvolutionaryGame {
    /// Payoff matrix (strategies × strategies)
    pub payoffs: Vec<Vec<f64>>,
    /// Population distribution over strategies
    pub population: Vec<f64>,
}

impl EvolutionaryGame {
    /// Create new evolutionary game
    pub fn new(payoffs: Vec<Vec<f64>>, initial_population: Vec<f64>) -> Self {
        Self {
            payoffs,
            population: initial_population,
        }
    }

    /// Simulate replicator dynamics for one step
    pub fn step(&mut self, dt: f64) {
        let n = self.population.len();
        let mut new_pop = vec![0.0; n];

        // Compute average payoff for each strategy
        let mut strategy_payoffs = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                strategy_payoffs[i] += self.payoffs[i][j] * self.population[j];
            }
        }

        // Compute average population payoff
        let avg_payoff: f64 = strategy_payoffs
            .iter()
            .zip(&self.population)
            .map(|(&pay, &pop)| pay * pop)
            .sum();

        // Update populations using replicator equation
        // dx_i/dt = x_i * (payoff_i - avg_payoff)
        for i in 0..n {
            let growth = self.population[i] * (strategy_payoffs[i] - avg_payoff);
            new_pop[i] = (self.population[i] + dt * growth).max(0.0);
        }

        // Normalize
        let total: f64 = new_pop.iter().sum();
        if total > 0.0 {
            for p in &mut new_pop {
                *p /= total;
            }
        }

        self.population = new_pop;
    }

    /// Check if current population is an ESS (Evolutionary Stable Strategy)
    pub fn is_ess(&self) -> bool {
        // Simplified check: dominant strategy
        let n = self.population.len();

        for i in 0..n {
            if self.population[i] > 0.5 {
                // Check if strategy i dominates
                for j in 0..n {
                    if i != j && self.payoffs[i][i] <= self.payoffs[j][i] {
                        return false;
                    }
                }
                return true;
            }
        }

        false
    }
}

