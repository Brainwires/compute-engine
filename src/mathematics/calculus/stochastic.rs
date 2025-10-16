use rand::prelude::*;
use rand_distr::Normal;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
pub struct StochasticCalculusResult {
    pub success: bool,
    pub operation: String,
    pub result: Option<Value>,
    pub error: Option<String>,
}

// Itô integral approximation using Euler-Maruyama scheme
pub fn ito_integral(integrand_values: &[f64], brownian_increments: &[f64], dt: f64) -> f64 {
    if integrand_values.len() != brownian_increments.len() {
        return 0.0;
    }

    let mut integral = 0.0;
    for i in 0..integrand_values.len() {
        integral += integrand_values[i] * brownian_increments[i];
    }

    integral
}

// Stratonovich integral approximation
pub fn stratonovich_integral(
    integrand_values: &[f64],
    brownian_increments: &[f64],
    dt: f64,
) -> f64 {
    if integrand_values.len() != brownian_increments.len() || integrand_values.len() < 2 {
        return 0.0;
    }

    let mut integral = 0.0;

    for i in 0..integrand_values.len() - 1 {
        let midpoint_integrand = (integrand_values[i] + integrand_values[i + 1]) / 2.0;
        integral += midpoint_integrand * brownian_increments[i];
    }

    integral
}

// Generate Brownian motion path
pub fn generate_brownian_motion(
    t_max: f64,
    n_steps: usize,
    initial_value: f64,
    drift: f64,
    volatility: f64,
) -> Vec<(f64, f64)> {
    let dt = t_max / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut path = Vec::with_capacity(n_steps + 1);
    let mut current_value = initial_value;
    path.push((0.0, current_value));

    for i in 1..=n_steps {
        let dw = normal.sample(&mut rng) * sqrt_dt;
        current_value += drift * dt + volatility * dw;
        path.push((i as f64 * dt, current_value));
    }

    path
}

// Solve SDE using Euler-Maruyama method
pub fn solve_sde_euler_maruyama(
    drift_function: fn(f64, f64) -> f64,
    diffusion_function: fn(f64, f64) -> f64,
    initial_value: f64,
    t_max: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let dt = t_max / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut solution = Vec::with_capacity(n_steps + 1);
    let mut current_value = initial_value;
    let mut current_time = 0.0;

    solution.push((current_time, current_value));

    for _ in 1..=n_steps {
        let drift = drift_function(current_time, current_value);
        let diffusion = diffusion_function(current_time, current_value);
        let dw = normal.sample(&mut rng) * sqrt_dt;

        current_value += drift * dt + diffusion * dw;
        current_time += dt;

        solution.push((current_time, current_value));
    }

    solution
}

// Milstein method for SDE solving (higher order accuracy)
pub fn solve_sde_milstein(
    drift_function: fn(f64, f64) -> f64,
    diffusion_function: fn(f64, f64) -> f64,
    diffusion_derivative: fn(f64, f64) -> f64,
    initial_value: f64,
    t_max: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let dt = t_max / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut solution = Vec::with_capacity(n_steps + 1);
    let mut current_value = initial_value;
    let mut current_time = 0.0;

    solution.push((current_time, current_value));

    for _ in 1..=n_steps {
        let drift = drift_function(current_time, current_value);
        let diffusion = diffusion_function(current_time, current_value);
        let diffusion_deriv = diffusion_derivative(current_time, current_value);
        let dw = normal.sample(&mut rng) * sqrt_dt;

        // Milstein correction term
        let correction = 0.5 * diffusion * diffusion_deriv * (dw * dw - dt);

        current_value += drift * dt + diffusion * dw + correction;
        current_time += dt;

        solution.push((current_time, current_value));
    }

    solution
}

// Black-Scholes option pricing using Monte Carlo
pub fn black_scholes_monte_carlo(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_expiry: f64,
    n_simulations: usize,
    option_type: &str,
) -> f64 {
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let drift = (risk_free_rate - 0.5 * volatility * volatility) * time_to_expiry;
    let diffusion = volatility * time_to_expiry.sqrt();

    let mut payoff_sum = 0.0;

    for _ in 0..n_simulations {
        let z = normal.sample(&mut rng);
        let final_price = spot_price * (drift + diffusion * z).exp();

        let payoff = match option_type {
            "call" => (final_price - strike_price).max(0.0),
            "put" => (strike_price - final_price).max(0.0),
            _ => 0.0,
        };

        payoff_sum += payoff;
    }

    let option_value =
        (payoff_sum / n_simulations as f64) * (-risk_free_rate * time_to_expiry).exp();
    option_value
}

// Ornstein-Uhlenbeck process simulation
pub fn ornstein_uhlenbeck_process(
    theta: f64,
    mu: f64,
    sigma: f64,
    initial_value: f64,
    t_max: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let dt = t_max / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut path = Vec::with_capacity(n_steps + 1);
    let mut current_value = initial_value;
    path.push((0.0, current_value));

    for i in 1..=n_steps {
        let dw = normal.sample(&mut rng) * sqrt_dt;
        current_value += theta * (mu - current_value) * dt + sigma * dw;
        path.push((i as f64 * dt, current_value));
    }

    path
}

// Predefined drift and diffusion functions for common SDEs
fn linear_drift(t: f64, x: f64) -> f64 {
    -0.1 * x // Mean-reverting with rate 0.1
}

fn constant_diffusion(t: f64, x: f64) -> f64 {
    0.2 // Constant volatility
}

fn linear_diffusion_derivative(t: f64, x: f64) -> f64 {
    0.0 // Derivative of constant is zero
}

fn geometric_drift(t: f64, x: f64) -> f64 {
    0.05 * x // Geometric Brownian motion drift
}

fn geometric_diffusion(t: f64, x: f64) -> f64 {
    0.2 * x // Geometric Brownian motion volatility
}

fn geometric_diffusion_derivative(t: f64, x: f64) -> f64 {
    0.2 // Derivative of 0.2*x with respect to x
}

// Handler functions
pub fn handle_ito_integral(params: &HashMap<String, Value>) -> StochasticCalculusResult {
    let integrand_values: Vec<f64> = params
        .get("integrand_values")
        .and_then(|v| v.as_array())
        .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect())
        .unwrap_or_else(|| vec![1.0, 1.0, 1.0]);

    let dt = params.get("dt").and_then(|v| v.as_f64()).unwrap_or(0.01);

    // Generate random Brownian increments for demo
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, dt.sqrt()).unwrap();
    let brownian_increments: Vec<f64> = (0..integrand_values.len())
        .map(|_| normal.sample(&mut rng))
        .collect();

    let result = ito_integral(&integrand_values, &brownian_increments, dt);

    StochasticCalculusResult {
        success: true,
        operation: "ito_integral".to_string(),
        result: Some(serde_json::json!({
            "integral_value": result,
            "integrand_values": integrand_values,
            "brownian_increments": brownian_increments,
            "dt": dt
        })),
        error: None,
    }
}

pub fn handle_stratonovich_integral(params: &HashMap<String, Value>) -> StochasticCalculusResult {
    let integrand_values: Vec<f64> = params
        .get("integrand_values")
        .and_then(|v| v.as_array())
        .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect())
        .unwrap_or_else(|| vec![1.0, 1.0, 1.0]);

    let dt = params.get("dt").and_then(|v| v.as_f64()).unwrap_or(0.01);

    // Generate random Brownian increments for demo
    let mut rng = thread_rng();
    let normal = Normal::new(0.0, dt.sqrt()).unwrap();
    let brownian_increments: Vec<f64> = (0..integrand_values.len())
        .map(|_| normal.sample(&mut rng))
        .collect();

    let result = stratonovich_integral(&integrand_values, &brownian_increments, dt);

    StochasticCalculusResult {
        success: true,
        operation: "stratonovich_integral".to_string(),
        result: Some(serde_json::json!({
            "integral_value": result,
            "integrand_values": integrand_values,
            "brownian_increments": brownian_increments,
            "dt": dt
        })),
        error: None,
    }
}

pub fn handle_sde_solution(params: &HashMap<String, Value>) -> StochasticCalculusResult {
    let sde_type = params
        .get("sde_type")
        .and_then(|v| v.as_str())
        .unwrap_or("geometric_brownian");
    let initial_value = params
        .get("initial_value")
        .and_then(|v| v.as_f64())
        .unwrap_or(1.0);
    let t_max = params.get("t_max").and_then(|v| v.as_f64()).unwrap_or(1.0);
    let n_steps = params
        .get("n_steps")
        .and_then(|v| v.as_u64())
        .unwrap_or(1000) as usize;
    let method = params
        .get("method")
        .and_then(|v| v.as_str())
        .unwrap_or("euler_maruyama");

    let solution = match (sde_type, method) {
        ("geometric_brownian", "euler_maruyama") => solve_sde_euler_maruyama(
            geometric_drift,
            geometric_diffusion,
            initial_value,
            t_max,
            n_steps,
        ),
        ("geometric_brownian", "milstein") => solve_sde_milstein(
            geometric_drift,
            geometric_diffusion,
            geometric_diffusion_derivative,
            initial_value,
            t_max,
            n_steps,
        ),
        ("ornstein_uhlenbeck", _) => {
            let theta = params.get("theta").and_then(|v| v.as_f64()).unwrap_or(0.1);
            let mu = params.get("mu").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let sigma = params.get("sigma").and_then(|v| v.as_f64()).unwrap_or(0.2);
            ornstein_uhlenbeck_process(theta, mu, sigma, initial_value, t_max, n_steps)
        }
        ("mean_reverting", "euler_maruyama") => solve_sde_euler_maruyama(
            linear_drift,
            constant_diffusion,
            initial_value,
            t_max,
            n_steps,
        ),
        _ => solve_sde_euler_maruyama(
            geometric_drift,
            geometric_diffusion,
            initial_value,
            t_max,
            n_steps,
        ),
    };

    StochasticCalculusResult {
        success: true,
        operation: "sde_solution".to_string(),
        result: Some(serde_json::json!({
            "solution": solution,
            "sde_type": sde_type,
            "method": method,
            "initial_value": initial_value,
            "t_max": t_max,
            "n_steps": n_steps
        })),
        error: None,
    }
}

pub fn handle_stochastic_calculus_ops(params: &HashMap<String, Value>) -> StochasticCalculusResult {
    let operation = params
        .get("operation")
        .and_then(|v| v.as_str())
        .unwrap_or("brownian_motion");

    match operation {
        "brownian_motion" => {
            let t_max = params.get("t_max").and_then(|v| v.as_f64()).unwrap_or(1.0);
            let n_steps = params
                .get("n_steps")
                .and_then(|v| v.as_u64())
                .unwrap_or(1000) as usize;
            let initial_value = params
                .get("initial_value")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let drift = params.get("drift").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let volatility = params
                .get("volatility")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);

            let path = generate_brownian_motion(t_max, n_steps, initial_value, drift, volatility);

            StochasticCalculusResult {
                success: true,
                operation: "brownian_motion".to_string(),
                result: Some(serde_json::json!({
                    "path": path,
                    "t_max": t_max,
                    "n_steps": n_steps,
                    "parameters": {
                        "initial_value": initial_value,
                        "drift": drift,
                        "volatility": volatility
                    }
                })),
                error: None,
            }
        }
        "black_scholes_monte_carlo" => {
            let spot_price = params
                .get("spot_price")
                .and_then(|v| v.as_f64())
                .unwrap_or(100.0);
            let strike_price = params
                .get("strike_price")
                .and_then(|v| v.as_f64())
                .unwrap_or(100.0);
            let risk_free_rate = params
                .get("risk_free_rate")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.05);
            let volatility = params
                .get("volatility")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.2);
            let time_to_expiry = params
                .get("time_to_expiry")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);
            let n_simulations = params
                .get("n_simulations")
                .and_then(|v| v.as_u64())
                .unwrap_or(10000) as usize;
            let option_type = params
                .get("option_type")
                .and_then(|v| v.as_str())
                .unwrap_or("call");

            let option_value = black_scholes_monte_carlo(
                spot_price,
                strike_price,
                risk_free_rate,
                volatility,
                time_to_expiry,
                n_simulations,
                option_type,
            );

            StochasticCalculusResult {
                success: true,
                operation: "black_scholes_monte_carlo".to_string(),
                result: Some(serde_json::json!({
                    "option_value": option_value,
                    "parameters": {
                        "spot_price": spot_price,
                        "strike_price": strike_price,
                        "risk_free_rate": risk_free_rate,
                        "volatility": volatility,
                        "time_to_expiry": time_to_expiry,
                        "n_simulations": n_simulations,
                        "option_type": option_type
                    }
                })),
                error: None,
            }
        }
        _ => StochasticCalculusResult {
            success: false,
            operation: operation.to_string(),
            result: None,
            error: Some(format!(
                "Unknown stochastic calculus operation: {}",
                operation
            )),
        },
    }
}
