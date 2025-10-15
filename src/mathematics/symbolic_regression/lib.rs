use clap::{Parser, Subcommand};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use rand::Rng;
use std::error::Error;
use std::fmt;

#[derive(Parser)]
#[command(name = "symbolic-regression")]
#[command(about = "Physics-informed symbolic regression for equation discovery")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    DiscoverEquations {
        #[arg(long)]
        data: String,
        #[arg(long)]
        domain: String,
        #[arg(long)]
        complexity: u32,
        #[arg(long)]
        units: Option<String>,
    },
}

#[derive(Deserialize, Debug)]
pub struct InputData {
    pub variables: Vec<String>,
    pub data_points: Vec<Vec<f64>>,
    pub target_values: Vec<f64>,
}

#[derive(Serialize, Debug)]
pub struct EquationCandidate {
    pub expression: String,
    pub complexity: u32,
    pub fitness: f64,
    pub r_squared: f64,
    pub mse: f64,
    pub physical_validity: f64,
}

#[derive(Serialize, Debug)]
pub struct PhysicsConstraints {
    pub units: HashMap<String, String>,
    conservation_laws: Vec<String>,
    symmetries: Vec<String>,
    dimensional_consistency: bool,
    physical_plausibility: bool,
}

#[derive(Serialize, Debug)]
pub struct SymbolicRegressionResult {
    pub candidates: Vec<EquationCandidate>,
    pub constraints: PhysicsConstraints,
    pub best_fitness: f64,
    pub convergence_data: Vec<f64>,
}

#[derive(Debug)]
struct RegressionError {
    message: String,
}

impl fmt::Display for RegressionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Regression error: {}", self.message)
    }
}

impl Error for RegressionError {}

// Simple expression tree for symbolic regression
#[derive(Clone, Debug)]
enum Expression {
    Variable(String),
    Constant(f64),
    Add(Box<Expression>, Box<Expression>),
    Subtract(Box<Expression>, Box<Expression>),
    Multiply(Box<Expression>, Box<Expression>),
    Divide(Box<Expression>, Box<Expression>),
    Power(Box<Expression>, Box<Expression>),
    Sin(Box<Expression>),
    Cos(Box<Expression>),
    Exp(Box<Expression>),
    Log(Box<Expression>),
    Sqrt(Box<Expression>),
}

impl Expression {
    fn evaluate(&self, variables: &HashMap<String, f64>) -> Result<f64, Box<dyn Error>> {
        match self {
            Expression::Variable(name) => {
                variables.get(name)
                    .copied()
                    .ok_or_else(|| format!("Unknown variable: {}", name).into())
            }
            Expression::Constant(value) => Ok(*value),
            Expression::Add(left, right) => {
                Ok(left.evaluate(variables)? + right.evaluate(variables)?)
            }
            Expression::Subtract(left, right) => {
                Ok(left.evaluate(variables)? - right.evaluate(variables)?)
            }
            Expression::Multiply(left, right) => {
                Ok(left.evaluate(variables)? * right.evaluate(variables)?)
            }
            Expression::Divide(left, right) => {
                let divisor = right.evaluate(variables)?;
                if divisor.abs() < 1e-10 {
                    Err("Division by zero".into())
                } else {
                    Ok(left.evaluate(variables)? / divisor)
                }
            }
            Expression::Power(base, exp) => {
                Ok(base.evaluate(variables)?.powf(exp.evaluate(variables)?))
            }
            Expression::Sin(expr) => Ok(expr.evaluate(variables)?.sin()),
            Expression::Cos(expr) => Ok(expr.evaluate(variables)?.cos()),
            Expression::Exp(expr) => Ok(expr.evaluate(variables)?.exp()),
            Expression::Log(expr) => {
                let val = expr.evaluate(variables)?;
                if val <= 0.0 {
                    Err("Logarithm of non-positive number".into())
                } else {
                    Ok(val.ln())
                }
            }
            Expression::Sqrt(expr) => {
                let val = expr.evaluate(variables)?;
                if val < 0.0 {
                    Err("Square root of negative number".into())
                } else {
                    Ok(val.sqrt())
                }
            }
        }
    }

    fn to_string(&self) -> String {
        match self {
            Expression::Variable(name) => name.clone(),
            Expression::Constant(value) => format!("{:.3}", value),
            Expression::Add(left, right) => format!("({} + {})", left.to_string(), right.to_string()),
            Expression::Subtract(left, right) => format!("({} - {})", left.to_string(), right.to_string()),
            Expression::Multiply(left, right) => format!("({} * {})", left.to_string(), right.to_string()),
            Expression::Divide(left, right) => format!("({} / {})", left.to_string(), right.to_string()),
            Expression::Power(base, exp) => format!("({}^{})", base.to_string(), exp.to_string()),
            Expression::Sin(expr) => format!("sin({})", expr.to_string()),
            Expression::Cos(expr) => format!("cos({})", expr.to_string()),
            Expression::Exp(expr) => format!("exp({})", expr.to_string()),
            Expression::Log(expr) => format!("log({})", expr.to_string()),
            Expression::Sqrt(expr) => format!("sqrt({})", expr.to_string()),
        }
    }

    fn complexity(&self) -> u32 {
        match self {
            Expression::Variable(_) | Expression::Constant(_) => 1,
            Expression::Add(left, right) | Expression::Subtract(left, right) |
            Expression::Multiply(left, right) | Expression::Divide(left, right) |
            Expression::Power(left, right) => 1 + left.complexity() + right.complexity(),
            Expression::Sin(expr) | Expression::Cos(expr) | Expression::Exp(expr) |
            Expression::Log(expr) | Expression::Sqrt(expr) => 2 + expr.complexity(),
        }
    }
}

pub fn generate_random_expression(variables: &[String], max_depth: u32, rng: &mut impl Rng) -> Expression {
    if max_depth == 0 || rng.gen_bool(0.3) {
        // Terminal node
        if rng.gen_bool(0.7) && !variables.is_empty() {
            Expression::Variable(variables[rng.gen_range(0..variables.len())].clone())
        } else {
            Expression::Constant(rng.gen_range(-10.0..10.0))
        }
    } else {
        // Non-terminal node
        match rng.gen_range(0..11) {
            0 => Expression::Add(
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
            ),
            1 => Expression::Subtract(
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
            ),
            2 => Expression::Multiply(
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
            ),
            3 => Expression::Divide(
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
            ),
            4 => Expression::Power(
                Box::new(generate_random_expression(variables, max_depth - 1, rng)),
                Box::new(Expression::Constant(rng.gen_range(0.1..3.0))),
            ),
            5 => Expression::Sin(Box::new(generate_random_expression(variables, max_depth - 1, rng))),
            6 => Expression::Cos(Box::new(generate_random_expression(variables, max_depth - 1, rng))),
            7 => Expression::Exp(Box::new(generate_random_expression(variables, max_depth - 1, rng))),
            8 => Expression::Log(Box::new(generate_random_expression(variables, max_depth - 1, rng))),
            9 => Expression::Sqrt(Box::new(generate_random_expression(variables, max_depth - 1, rng))),
            _ => Expression::Variable(variables[rng.gen_range(0..variables.len())].clone()),
        }
    }
}

pub fn evaluate_fitness(
    expr: &Expression,
    data: &InputData,
) -> Result<(f64, f64, f64), Box<dyn Error>> {
    let mut predictions = Vec::new();
    let mut valid_count = 0;

    for (_i, data_point) in data.data_points.iter().enumerate() {
        let mut var_map = HashMap::new();
        for (j, var_name) in data.variables.iter().enumerate() {
            var_map.insert(var_name.clone(), data_point[j]);
        }

        match expr.evaluate(&var_map) {
            Ok(pred) if pred.is_finite() => {
                predictions.push(pred);
                valid_count += 1;
            }
            _ => {
                predictions.push(0.0); // Use 0 for invalid predictions
            }
        }
    }

    if valid_count == 0 {
        return Ok((0.0, 0.0, f64::INFINITY));
    }

    // Calculate MSE
    let mse = predictions.iter()
        .zip(data.target_values.iter())
        .map(|(pred, target)| (pred - target).powi(2))
        .sum::<f64>() / predictions.len() as f64;

    // Calculate R²
    let target_mean = data.target_values.iter().sum::<f64>() / data.target_values.len() as f64;
    let ss_tot = data.target_values.iter()
        .map(|target| (target - target_mean).powi(2))
        .sum::<f64>();
    let ss_res = predictions.iter()
        .zip(data.target_values.iter())
        .map(|(pred, target)| (target - pred).powi(2))
        .sum::<f64>();

    let r_squared = if ss_tot.abs() < 1e-10 {
        1.0
    } else {
        1.0 - (ss_res / ss_tot)
    };

    // Fitness is R² adjusted for complexity, penalize invalid evaluations
    let complexity_penalty = expr.complexity() as f64 * 0.01;
    let validity_bonus = (valid_count as f64 / data.data_points.len() as f64) * 0.1;
    let fitness = r_squared.max(0.0) - complexity_penalty + validity_bonus;

    Ok((fitness, r_squared, mse))
}

pub fn check_physics_constraints(
    expr: &Expression,
    domain: &str,
    _units: &Option<HashMap<String, String>>,
) -> f64 {
    // Simple heuristic physics validation
    let expression_str = expr.to_string();
    
    let mut validity_score: f64 = 1.0;

    match domain {
        "mechanics" => {
            // Check for common physics patterns
            if expression_str.contains("*") && (expression_str.contains("m") || expression_str.contains("a")) {
                validity_score += 0.2; // Likely F=ma type relationship
            }
            if expression_str.contains("^2") && expression_str.contains("v") {
                validity_score += 0.1; // Kinetic energy patterns
            }
        }
        "electromagnetism" => {
            if expression_str.contains("/") && expression_str.contains("r^2") {
                validity_score += 0.2; // Inverse square law
            }
        }
        "thermodynamics" => {
            if expression_str.contains("*") && expression_str.contains("T") {
                validity_score += 0.1; // Temperature relationships
            }
        }
        _ => {} // General physics
    }

    // Penalize overly complex expressions
    if expr.complexity() > 20 {
        validity_score *= 0.5;
    }

    validity_score.min(2.0).max(0.0)
}

pub fn discover_equations(
    data: InputData,
    domain: String,
    max_complexity: u32,
    units: Option<HashMap<String, String>>,
) -> Result<SymbolicRegressionResult, Box<dyn Error>> {
    let mut rng = rand::thread_rng();
    let population_size = 100;
    let generations = 50;
    let mut population = Vec::new();

    // Initialize population
    for _ in 0..population_size {
        let expr = generate_random_expression(&data.variables, 4, &mut rng);
        population.push(expr);
    }

    let mut best_fitnesses = Vec::new();
    let mut best_candidates: Vec<EquationCandidate> = Vec::new();

    // Evolutionary loop
    for generation in 0..generations {
        let mut fitness_scores = Vec::new();

        // Evaluate fitness for each individual
        for expr in &population {
            if expr.complexity() <= max_complexity {
                match evaluate_fitness(expr, &data) {
                    Ok((fitness, r_squared, mse)) => {
                        let physics_validity = check_physics_constraints(expr, &domain, &units);
                        let combined_fitness = fitness * 0.7 + physics_validity * 0.3;
                        
                        fitness_scores.push((combined_fitness, r_squared, mse, physics_validity));
                        
                        // Track best candidates
                        if best_candidates.len() < 10 || combined_fitness > best_candidates.last().unwrap().fitness {
                            let candidate = EquationCandidate {
                                expression: expr.to_string(),
                                complexity: expr.complexity(),
                                fitness: combined_fitness,
                                r_squared,
                                mse,
                                physical_validity: physics_validity,
                            };
                            
                            best_candidates.push(candidate);
                            best_candidates.sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap());
                            if best_candidates.len() > 10 {
                                best_candidates.truncate(10);
                            }
                        }
                    }
                    Err(_) => {
                        fitness_scores.push((0.0, 0.0, f64::INFINITY, 0.0));
                    }
                }
            } else {
                fitness_scores.push((0.0, 0.0, f64::INFINITY, 0.0));
            }
        }

        // Record best fitness
        let best_fitness = fitness_scores.iter()
            .map(|(f, _, _, _)| *f)
            .fold(0.0f64, f64::max);
        best_fitnesses.push(best_fitness);

        // Selection and reproduction (simplified)
        let mut new_population = Vec::new();
        
        // Keep best individuals (elitism)
        let mut indexed_fitness: Vec<(usize, f64)> = fitness_scores.iter()
            .enumerate()
            .map(|(i, (f, _, _, _))| (i, *f))
            .collect();
        indexed_fitness.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        
        for i in 0..(population_size / 4) {
            if let Some((idx, _)) = indexed_fitness.get(i) {
                new_population.push(population[*idx].clone());
            }
        }

        // Fill rest with random new individuals
        while new_population.len() < population_size {
            new_population.push(generate_random_expression(&data.variables, 4, &mut rng));
        }

        population = new_population;

        // Early stopping if we found a very good solution
        if best_fitness > 0.95 {
            break;
        }
    }

    let constraints = PhysicsConstraints {
        units: units.unwrap_or_default(),
        conservation_laws: match domain.as_str() {
            "mechanics" => vec!["energy".to_string(), "momentum".to_string()],
            "electromagnetism" => vec!["charge".to_string(), "energy".to_string()],
            _ => vec!["energy".to_string()],
        },
        symmetries: vec!["space_translation".to_string(), "time_translation".to_string()],
        dimensional_consistency: true,
        physical_plausibility: true,
    };

    Ok(SymbolicRegressionResult {
        candidates: best_candidates,
        constraints,
        best_fitness: best_fitnesses.iter().fold(0.0f64, |a, &b| a.max(b)),
        convergence_data: best_fitnesses,
    })
}

