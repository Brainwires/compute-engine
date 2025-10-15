//! Main function approximator using genetic programming

use super::expression::ExpressionNode;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FunctionApproximatorConfig {
    pub population_size: usize,
    pub generations: usize,
    pub complexity_penalty: f64,
    pub death_probability: f64,
    pub diversity_factor: usize,
    pub max_tree_depth: usize,
}

impl Default for FunctionApproximatorConfig {
    fn default() -> Self {
        Self {
            population_size: 600,
            generations: 100,
            complexity_penalty: 2.0,
            death_probability: 0.5,
            diversity_factor: 10,
            max_tree_depth: 5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Individual {
    pub expression: ExpressionNode,
    pub fitness: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ApproximationResult {
    pub best_functions: Vec<BestFunction>,
    pub generations_run: usize,
    pub final_population_stats: PopulationStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BestFunction {
    pub expression: String,
    pub fitness: f64,
    pub complexity: usize,
    pub predictions: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PopulationStats {
    pub min_fitness: f64,
    pub max_fitness: f64,
    pub avg_fitness: f64,
    pub fitness_distribution: Vec<f64>,
}

pub struct FunctionApproximator {
    config: FunctionApproximatorConfig,
    inputs: Vec<f64>,
    outputs: Vec<f64>,
    variables: Vec<String>,
}

impl FunctionApproximator {
    pub fn new(
        inputs: Vec<f64>,
        outputs: Vec<f64>,
        config: FunctionApproximatorConfig,
    ) -> Result<Self, String> {
        if inputs.len() != outputs.len() {
            return Err("Inputs and outputs must have the same length".to_string());
        }
        if inputs.is_empty() {
            return Err("Must provide at least one input/output pair".to_string());
        }

        Ok(Self {
            config,
            inputs,
            outputs,
            variables: vec!["x".to_string()],
        })
    }

    /// Calculate fitness (lower is better)
    fn calculate_fitness(&self, expr: &ExpressionNode) -> f64 {
        let mut total_error = 0.0;
        let mut success_count = 0;

        for (input, output) in self.inputs.iter().zip(self.outputs.iter()) {
            let mut vars = HashMap::new();
            vars.insert("x".to_string(), *input);

            match expr.evaluate(&vars) {
                Ok(prediction) => {
                    let error = (prediction - output).abs();
                    total_error += error;
                    success_count += 1;
                }
                Err(_) => {
                    // Penalize failed evaluations heavily
                    total_error += 1000.0;
                }
            }
        }

        if success_count == 0 {
            return f64::MAX;
        }

        // Add complexity penalty
        let complexity = expr.complexity() as f64;
        let mse = total_error / success_count as f64;

        mse + (complexity / self.config.complexity_penalty)
    }

    /// Run the genetic algorithm
    pub fn approximate(&self) -> Result<ApproximationResult, String> {
        let mut rng = rand::thread_rng();

        // Initialize population
        let mut population: Vec<Individual> = (0..self.config.population_size)
            .map(|_| {
                let expr = ExpressionNode::random(
                    &mut rng,
                    &self.variables,
                    self.config.max_tree_depth,
                    0,
                );
                let fitness = self.calculate_fitness(&expr);
                Individual {
                    expression: expr,
                    fitness,
                }
            })
            .collect();

        // Evolution loop
        for generation in 0..self.config.generations {
            // Sort by fitness (lower is better)
            population.sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap());

            // Calculate how many to keep
            let survivors_count =
                (population.len() as f64 * (1.0 - self.config.death_probability)) as usize;
            let survivors_count = survivors_count.max(10).min(population.len());

            // Keep best individuals
            population.truncate(survivors_count);

            // Reproduce to fill population
            while population.len() < self.config.population_size {
                // Select random parent from survivors (bias towards better fitness)
                let random_val = rng.sample::<f64, _>(rand::distributions::Standard);
                let parent_idx = (random_val.powf(2.0) * survivors_count as f64) as usize;
                let parent_idx = parent_idx.min(survivors_count - 1);

                let parent = &population[parent_idx];

                // Mutate parent
                let child_expr = parent.expression.mutate(&mut rng, &self.variables);
                let child_fitness = self.calculate_fitness(&child_expr);

                population.push(Individual {
                    expression: child_expr,
                    fitness: child_fitness,
                });
            }

            // Optional: log progress
            if generation % 10 == 0 || generation == self.config.generations - 1 {
                let best_fitness = population[0].fitness;
                eprintln!("Generation {}: Best fitness = {:.6}", generation, best_fitness);
            }
        }

        // Final sort
        population.sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap());

        // Extract best functions
        let num_best = 10.min(population.len());
        let best_functions: Vec<BestFunction> = population
            .iter()
            .take(num_best)
            .map(|individual| {
                let predictions: Vec<f64> = self
                    .inputs
                    .iter()
                    .map(|input| {
                        let mut vars = HashMap::new();
                        vars.insert("x".to_string(), *input);
                        individual.expression.evaluate(&vars).unwrap_or(f64::NAN)
                    })
                    .collect();

                BestFunction {
                    expression: individual.expression.to_string(),
                    fitness: individual.fitness,
                    complexity: individual.expression.complexity(),
                    predictions,
                }
            })
            .collect();

        // Calculate population stats
        let fitnesses: Vec<f64> = population.iter().map(|i| i.fitness).collect();
        let min_fitness = fitnesses.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_fitness = fitnesses.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let avg_fitness = fitnesses.iter().sum::<f64>() / fitnesses.len() as f64;

        Ok(ApproximationResult {
            best_functions,
            generations_run: self.config.generations,
            final_population_stats: PopulationStats {
                min_fitness,
                max_fitness,
                avg_fitness,
                fitness_distribution: fitnesses.into_iter().take(20).collect(),
            },
        })
    }
}
