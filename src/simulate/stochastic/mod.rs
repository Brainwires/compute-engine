//! Stochastic Process Simulations
//!
//! Brownian motion, geometric Brownian, Ornstein-Uhlenbeck, Poisson, Lévy, and more
//!
//! This module provides:
//! - Random walks and Brownian motion
//! - Monte Carlo simulations
//! - Markov chains
//! - Stochastic differential equations

use crate::engine::*;

// ============================================================================
// Legacy Types and Functions (from specialized/stochastic_processes)
// ============================================================================

use rand::prelude::*;
use rand_distr::{Normal, Uniform};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Legacy request/response structs for backwards compatibility
#[derive(Debug, Serialize, Deserialize)]
#[allow(dead_code)]
struct StochasticRequest {
    operation: String,
    parameters: HashMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
#[allow(dead_code)]
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

    let mut current_state_idx = params
        .states
        .iter()
        .position(|s| s == &params.initial_state)
        .ok_or("Initial state not found")?;

    let mut path = Vec::with_capacity(params.steps + 1);
    path.push(params.states[current_state_idx].clone());

    let mut rng = thread_rng();
    let uniform = Uniform::new(0.0, 1.0);

    for _ in 0..params.steps {
        let random_value = uniform.sample(&mut rng);
        let mut cumulative_prob = 0.0;

        for (next_state_idx, &prob) in params.transition_matrix[current_state_idx]
            .iter()
            .enumerate()
        {
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

pub fn compute_stochastic_integral_monte_carlo(
    params: StochasticIntegralParams,
) -> Result<f64, String> {
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

// Legacy request processor for backwards compatibility
#[allow(dead_code)]
fn process_request(input: &str) -> String {
    let request: StochasticRequest = match serde_json::from_str(input) {
        Ok(req) => req,
        Err(e) => {
            return serde_json::to_string(&StochasticResult {
                success: false,
                result: None,
                error: Some(format!("Invalid JSON: {}", e)),
            })
            .unwrap();
        }
    };

    let result = match request.operation.as_str() {
        "brownian_motion" => {
            let params: BrownianMotionParams = match serde_json::from_value(
                serde_json::Value::Object(request.parameters.into_iter().collect()),
            ) {
                Ok(p) => p,
                Err(e) => {
                    return serde_json::to_string(&StochasticResult {
                        success: false,
                        result: None,
                        error: Some(format!("Invalid parameters: {}", e)),
                    })
                    .unwrap();
                }
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
        }

        "markov_chain" => {
            let params: MarkovChainParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.into_iter().collect(),
            )) {
                Ok(p) => p,
                Err(e) => {
                    return serde_json::to_string(&StochasticResult {
                        success: false,
                        result: None,
                        error: Some(format!("Invalid parameters: {}", e)),
                    })
                    .unwrap();
                }
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
        }

        "stochastic_integral" => {
            let params: StochasticIntegralParams = match serde_json::from_value(
                serde_json::Value::Object(request.parameters.into_iter().collect()),
            ) {
                Ok(p) => p,
                Err(e) => {
                    return serde_json::to_string(&StochasticResult {
                        success: false,
                        result: None,
                        error: Some(format!("Invalid parameters: {}", e)),
                    })
                    .unwrap();
                }
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
        }

        "ornstein_uhlenbeck" => {
            let time_steps = request
                .parameters
                .get("time_steps")
                .and_then(|v| v.as_u64())
                .unwrap_or(1000) as usize;
            let dt = request
                .parameters
                .get("dt")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.01);
            let initial_value = request
                .parameters
                .get("initial_value")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let mean_reversion_rate = request
                .parameters
                .get("mean_reversion_rate")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);
            let long_term_mean = request
                .parameters
                .get("long_term_mean")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let volatility = request
                .parameters
                .get("volatility")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);

            match generate_ornstein_uhlenbeck_process(
                time_steps,
                dt,
                initial_value,
                mean_reversion_rate,
                long_term_mean,
                volatility,
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
        }

        "poisson_process" => {
            let rate = request
                .parameters
                .get("rate")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);
            let time_horizon = request
                .parameters
                .get("time_horizon")
                .and_then(|v| v.as_f64())
                .unwrap_or(10.0);

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
        }

        _ => StochasticResult {
            success: false,
            result: None,
            error: Some(format!("Unknown operation: {}", request.operation)),
        },
    };

    serde_json::to_string(&result).unwrap()
}

// Test module for legacy functions
#[cfg(test)]
#[path = "../../../tests/unit/specialized/stochastic_processes_tests.rs"]
mod legacy_tests;

// ============================================================================
// Engine-Based Stochastic Simulation (new API)
// ============================================================================

/// Simulate stochastic processes
pub fn simulate_stochastic(
    process: &StochasticProcess,
    input: &SimulateInput,
) -> ToolResult<SimulateOutput> {
    use crate::compute::calculus::stochastic;

    let range = input.range.ok_or("range [start, end] required")?;
    let steps = input.steps.unwrap_or(1000);
    let _num_paths = input.num_paths.unwrap_or(1);

    match process {
        StochasticProcess::BrownianMotion => {
            let initial_value = input.parameters.get("initial_value").unwrap_or(&0.0);
            let drift = input.parameters.get("drift").unwrap_or(&0.0);
            let volatility = input.parameters.get("volatility").unwrap_or(&1.0);

            let path = stochastic::generate_brownian_motion(
                range[1] - range[0],
                steps,
                *initial_value,
                *drift,
                *volatility,
            );

            let times: Vec<f64> = path.iter().map(|(t, _)| *t).collect();
            let values: Vec<f64> = path.iter().map(|(_, v)| *v).collect();

            let mut results = std::collections::HashMap::new();
            results.insert("W".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "brownian_motion",
                    "drift": drift,
                    "volatility": volatility
                })),
            })
        }

        StochasticProcess::GeometricBrownian => {
            let initial_value = input.parameters.get("initial_value").unwrap_or(&100.0);
            let drift = input.parameters.get("drift").unwrap_or(&0.05);
            let volatility = input.parameters.get("volatility").unwrap_or(&0.2);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut s = *initial_value;
            let mut t = range[0];
            times.push(t);
            values.push(s);

            for _ in 0..steps {
                let dw = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                s = s * (1.0 + drift * dt + volatility * dw);
                t += dt;
                times.push(t);
                values.push(s);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("S".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "geometric_brownian",
                    "drift": drift,
                    "volatility": volatility,
                    "sde": "dS = μS dt + σS dW"
                })),
            })
        }

        StochasticProcess::OrnsteinUhlenbeck => {
            let theta = input
                .parameters
                .get("theta")
                .ok_or("theta (mean reversion rate) required")?;
            let mu = input
                .parameters
                .get("mu")
                .ok_or("mu (long-term mean) required")?;
            let sigma = input
                .parameters
                .get("sigma")
                .ok_or("sigma (volatility) required")?;
            let initial_value = input.parameters.get("initial_value").unwrap_or(&0.0);

            let path = stochastic::ornstein_uhlenbeck_process(
                *theta,
                *mu,
                *sigma,
                *initial_value,
                range[1] - range[0],
                steps,
            );

            let times: Vec<f64> = path.iter().map(|(t, _)| *t).collect();
            let values: Vec<f64> = path.iter().map(|(_, v)| *v).collect();

            let mut results = std::collections::HashMap::new();
            results.insert("X".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "ornstein_uhlenbeck",
                    "theta": theta,
                    "mu": mu,
                    "sigma": sigma
                })),
            })
        }

        StochasticProcess::Poisson => {
            let lambda = input.parameters.get("lambda").unwrap_or(&1.0);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut n = 0.0;
            let mut t = range[0];
            times.push(t);
            values.push(n);

            for _ in 0..steps {
                let prob = lambda * dt;
                if rand::random::<f64>() < prob {
                    n += 1.0;
                }
                t += dt;
                times.push(t);
                values.push(n);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("N".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "poisson",
                    "lambda": lambda,
                    "formula": "P(N(t)=k) = (λt)^k e^(-λt) / k!"
                })),
            })
        }

        StochasticProcess::Levy => {
            let alpha = input.parameters.get("alpha").unwrap_or(&1.5);
            let beta = input.parameters.get("beta").unwrap_or(&0.0);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut x = 0.0;
            let mut t = range[0];
            times.push(t);
            values.push(x);

            for _ in 0..steps {
                let dx = dt.powf(1.0 / alpha) * (rand::random::<f64>() - 0.5) * (1.0 + beta);
                x += dx;
                t += dt;
                times.push(t);
                values.push(x);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("L".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "levy",
                    "alpha": alpha,
                    "beta": beta,
                    "note": "Stable Lévy process approximation"
                })),
            })
        }

        StochasticProcess::JumpDiffusion => {
            let initial_value = input.parameters.get("initial_value").unwrap_or(&100.0);
            let drift = input.parameters.get("drift").unwrap_or(&0.05);
            let volatility = input.parameters.get("volatility").unwrap_or(&0.2);
            let jump_intensity = input.parameters.get("jump_intensity").unwrap_or(&0.1);
            let jump_mean = input.parameters.get("jump_mean").unwrap_or(&0.0);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut s = *initial_value;
            let mut t = range[0];
            times.push(t);
            values.push(s);

            for _ in 0..steps {
                let dw = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let jump = if rand::random::<f64>() < jump_intensity * dt {
                    jump_mean * (rand::random::<f64>() - 0.5)
                } else {
                    0.0
                };

                s = s * (1.0 + drift * dt + volatility * dw + jump);
                t += dt;
                times.push(t);
                values.push(s);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("S".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "jump_diffusion",
                    "jump_intensity": jump_intensity,
                    "sde": "dS = μS dt + σS dW + S dJ"
                })),
            })
        }

        StochasticProcess::FractionalBrownian => {
            let hurst = input.parameters.get("hurst").unwrap_or(&0.7);
            let initial_value = input.parameters.get("initial_value").unwrap_or(&0.0);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut b = *initial_value;
            let mut t = range[0];
            times.push(t);
            values.push(b);

            for _ in 0..steps {
                let db = dt.powf(*hurst) * (rand::random::<f64>() - 0.5) * 2.0;
                b += db;
                t += dt;
                times.push(t);
                values.push(b);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("B_H".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "fractional_brownian",
                    "hurst": hurst,
                    "note": "H=0.5 is standard Brownian, H>0.5 is persistent, H<0.5 is anti-persistent"
                })),
            })
        }

        StochasticProcess::MeanReverting => {
            let kappa = input.parameters.get("kappa").unwrap_or(&0.5);
            let theta = input.parameters.get("theta").unwrap_or(&1.0);
            let sigma = input.parameters.get("sigma").unwrap_or(&0.2);
            let initial_value = input.parameters.get("initial_value").unwrap_or(&1.0);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut x = *initial_value;
            let mut t = range[0];
            times.push(t);
            values.push(x);

            for _ in 0..steps {
                let dw = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let dx = kappa * (theta - x) * dt + sigma * x.abs().sqrt() * dw;
                x += dx;
                x = x.max(0.0);
                t += dt;
                times.push(t);
                values.push(x);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("X".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "mean_reverting",
                    "kappa": kappa,
                    "theta": theta,
                    "sde": "dx = κ(θ - x)dt + σ√x dW"
                })),
            })
        }

        StochasticProcess::VarianceGamma => {
            let nu = input.parameters.get("nu").unwrap_or(&0.2);
            let theta_vg = input.parameters.get("theta").unwrap_or(&0.0);
            let sigma = input.parameters.get("sigma").unwrap_or(&0.3);

            let dt = (range[1] - range[0]) / steps as f64;
            let mut times = Vec::with_capacity(steps + 1);
            let mut values = Vec::with_capacity(steps + 1);

            let mut x = 0.0;
            let mut t = range[0];
            times.push(t);
            values.push(x);

            for _ in 0..steps {
                let gamma_increment = nu * dt;
                let dw = gamma_increment.sqrt() * (rand::random::<f64>() - 0.5) * 2.0;
                let dx = theta_vg * gamma_increment + sigma * dw;
                x += dx;
                t += dt;
                times.push(t);
                values.push(x);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("X_VG".to_string(), values);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "process": "variance_gamma",
                    "nu": nu,
                    "theta": theta_vg,
                    "sigma": sigma,
                    "note": "Pure jump Lévy process"
                })),
            })
        }
    }
}
