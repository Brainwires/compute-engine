use serde::{Deserialize, Serialize};
use rand::prelude::*;
use rand_distr::{Normal, Uniform};
use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
struct StochasticRequest {
    operation: String,
    parameters: HashMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
struct StochasticResult {
    success: bool,
    result: Option<serde_json::Value>,
    error: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BrownianMotionParams {
    pub time_steps: usize,
    pub dt: f64,
    pub initial_value: f64,
    pub drift: f64,
    pub volatility: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MarkovChainParams {
    pub states: Vec<String>,
    pub transition_matrix: Vec<Vec<f64>>,
    pub initial_state: String,
    pub steps: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct StochasticIntegralParams {
    pub integrand_function: String,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub time_steps: usize,
    pub monte_carlo_samples: usize,
}

pub fn generate_brownian_motion(params: BrownianMotionParams) -> Result<Vec<f64>, String> {
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, (params.dt).sqrt()).map_err(|e| e.to_string())?;
    
    let mut path = Vec::with_capacity(params.time_steps + 1);
    path.push(params.initial_value);
    
    let mut current_value = params.initial_value;
    
    for _ in 0..params.time_steps {
        let dw = normal.sample(&mut rng);
        let drift_term = params.drift * params.dt;
        let stochastic_term = params.volatility * dw;
        current_value += drift_term + stochastic_term;
        path.push(current_value);
    }
    
    Ok(path)
}

pub fn simulate_markov_chain(params: MarkovChainParams) -> Result<Vec<String>, String> {
    if params.transition_matrix.len() != params.states.len() {
        return Err("Transition matrix dimensions don't match states".to_string());
    }
    
    // Validate transition matrix rows sum to 1
    for (i, row) in params.transition_matrix.iter().enumerate() {
        if row.len() != params.states.len() {
            return Err(format!("Row {} has wrong length", i));
        }
        let sum: f64 = row.iter().sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(format!("Row {} doesn't sum to 1.0", i));
        }
    }
    
    let mut current_state_idx = params.states.iter()
        .position(|s| s == &params.initial_state)
        .ok_or("Initial state not found")?;
    
    let mut path = Vec::with_capacity(params.steps + 1);
    path.push(params.states[current_state_idx].clone());
    
    let mut rng = thread_rng();
    let uniform = Uniform::new(0.0, 1.0);
    
    for _ in 0..params.steps {
        let random_value = uniform.sample(&mut rng);
        let mut cumulative_prob = 0.0;
        
        for (next_state_idx, &prob) in params.transition_matrix[current_state_idx].iter().enumerate() {
            cumulative_prob += prob;
            if random_value <= cumulative_prob {
                current_state_idx = next_state_idx;
                break;
            }
        }
        
        path.push(params.states[current_state_idx].clone());
    }
    
    Ok(path)
}

pub fn compute_stochastic_integral_monte_carlo(params: StochasticIntegralParams) -> Result<f64, String> {
    let dt = (params.upper_bound - params.lower_bound) / params.time_steps as f64;
    
    // Simplified Monte Carlo estimation for stochastic integral
    // For complex functions, would need expression parser
    let samples: Vec<f64> = (0..params.monte_carlo_samples)
        .into_par_iter()
        .map(|_| {
            let mut rng = thread_rng();
            let mut integral_sum = 0.0;
            let normal = Normal::new(0.0, dt.sqrt()).unwrap();
            
            for i in 0..params.time_steps {
                let t = params.lower_bound + i as f64 * dt;
                let dw = normal.sample(&mut rng);
                
                // Simple example: integrand = t (would need parser for general functions)
                let integrand_value = match params.integrand_function.as_str() {
                    "t" => t,
                    "t*t" => t * t,
                    "1" => 1.0,
                    _ => t, // Default fallback
                };
                
                integral_sum += integrand_value * dw;
            }
            integral_sum
        })
        .collect();
    
    Ok(samples.iter().sum::<f64>() / samples.len() as f64)
}

pub fn generate_ornstein_uhlenbeck_process(
    time_steps: usize,
    dt: f64,
    initial_value: f64,
    mean_reversion_rate: f64,
    long_term_mean: f64,
    volatility: f64,
) -> Result<Vec<f64>, String> {
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, (dt).sqrt()).map_err(|e| e.to_string())?;
    
    let mut path = Vec::with_capacity(time_steps + 1);
    path.push(initial_value);
    
    let mut current_value = initial_value;
    
    for _ in 0..time_steps {
        let dw = normal.sample(&mut rng);
        let drift_term = mean_reversion_rate * (long_term_mean - current_value) * dt;
        let stochastic_term = volatility * dw;
        current_value += drift_term + stochastic_term;
        path.push(current_value);
    }
    
    Ok(path)
}

pub fn simulate_poisson_process(rate: f64, time_horizon: f64) -> Result<Vec<f64>, String> {
    let mut rng = thread_rng();
    let exponential = rand_distr::Exp::new(rate).map_err(|e| e.to_string())?;
    
    let mut events = Vec::new();
    let mut current_time = 0.0;
    
    while current_time < time_horizon {
        let inter_arrival_time = exponential.sample(&mut rng);
        current_time += inter_arrival_time;
        if current_time < time_horizon {
            events.push(current_time);
        }
    }
    
    Ok(events)
}

fn process_request(input: &str) -> String {
    let request: StochasticRequest = match serde_json::from_str(input) {
        Ok(req) => req,
        Err(e) => return serde_json::to_string(&StochasticResult {
            success: false,
            result: None,
            error: Some(format!("Invalid JSON: {}", e)),
        }).unwrap(),
    };

    let result = match request.operation.as_str() {
        "brownian_motion" => {
            let params: BrownianMotionParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.into_iter().collect()
            )) {
                Ok(p) => p,
                Err(e) => return serde_json::to_string(&StochasticResult {
                    success: false,
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                }).unwrap(),
            };
            
            match generate_brownian_motion(params) {
                Ok(path) => StochasticResult {
                    success: true,
                    result: Some(serde_json::json!({ "path": path })),
                    error: None,
                },
                Err(e) => StochasticResult {
                    success: false,
                    result: None,
                    error: Some(e),
                },
            }
        },
        
        "markov_chain" => {
            let params: MarkovChainParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.into_iter().collect()
            )) {
                Ok(p) => p,
                Err(e) => return serde_json::to_string(&StochasticResult {
                    success: false,
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                }).unwrap(),
            };
            
            match simulate_markov_chain(params) {
                Ok(path) => StochasticResult {
                    success: true,
                    result: Some(serde_json::json!({ "path": path })),
                    error: None,
                },
                Err(e) => StochasticResult {
                    success: false,
                    result: None,
                    error: Some(e),
                },
            }
        },
        
        "stochastic_integral" => {
            let params: StochasticIntegralParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.into_iter().collect()
            )) {
                Ok(p) => p,
                Err(e) => return serde_json::to_string(&StochasticResult {
                    success: false,
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                }).unwrap(),
            };
            
            match compute_stochastic_integral_monte_carlo(params) {
                Ok(integral_value) => StochasticResult {
                    success: true,
                    result: Some(serde_json::json!({ "integral": integral_value })),
                    error: None,
                },
                Err(e) => StochasticResult {
                    success: false,
                    result: None,
                    error: Some(e),
                },
            }
        },
        
        "ornstein_uhlenbeck" => {
            let time_steps = request.parameters.get("time_steps")
                .and_then(|v| v.as_u64()).unwrap_or(1000) as usize;
            let dt = request.parameters.get("dt")
                .and_then(|v| v.as_f64()).unwrap_or(0.01);
            let initial_value = request.parameters.get("initial_value")
                .and_then(|v| v.as_f64()).unwrap_or(0.0);
            let mean_reversion_rate = request.parameters.get("mean_reversion_rate")
                .and_then(|v| v.as_f64()).unwrap_or(1.0);
            let long_term_mean = request.parameters.get("long_term_mean")
                .and_then(|v| v.as_f64()).unwrap_or(0.0);
            let volatility = request.parameters.get("volatility")
                .and_then(|v| v.as_f64()).unwrap_or(1.0);
            
            match generate_ornstein_uhlenbeck_process(
                time_steps, dt, initial_value, mean_reversion_rate, long_term_mean, volatility
            ) {
                Ok(path) => StochasticResult {
                    success: true,
                    result: Some(serde_json::json!({ "path": path })),
                    error: None,
                },
                Err(e) => StochasticResult {
                    success: false,
                    result: None,
                    error: Some(e),
                },
            }
        },
        
        "poisson_process" => {
            let rate = request.parameters.get("rate")
                .and_then(|v| v.as_f64()).unwrap_or(1.0);
            let time_horizon = request.parameters.get("time_horizon")
                .and_then(|v| v.as_f64()).unwrap_or(10.0);
            
            match simulate_poisson_process(rate, time_horizon) {
                Ok(events) => StochasticResult {
                    success: true,
                    result: Some(serde_json::json!({ "events": events, "count": events.len() })),
                    error: None,
                },
                Err(e) => StochasticResult {
                    success: false,
                    result: None,
                    error: Some(e),
                },
            }
        },
        
        _ => StochasticResult {
            success: false,
            result: None,
            error: Some(format!("Unknown operation: {}", request.operation)),
        },
    };

    serde_json::to_string(&result).unwrap()
}

