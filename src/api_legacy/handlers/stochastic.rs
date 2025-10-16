use super::super::{ComputationRequest, ComputationResponse};
use serde_json::{json, Value};

pub fn process_stochastic_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::stochastic_processes::*;

    let result = match request.operation.as_str() {
        "brownian_motion" => {
            let params: BrownianMotionParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(p) => p,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match generate_brownian_motion(params) {
                Ok(path) => Ok(json!({ "path": path })),
                Err(e) => Err(e),
            }
        },
        "markov_chain" => {
            let params: MarkovChainParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(p) => p,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match simulate_markov_chain(params) {
                Ok(states) => Ok(json!({ "states": states })),
                Err(e) => Err(e),
            }
        },
        "poisson_process" => {
            let time_horizon = request.parameters.get("time_horizon")
                .and_then(|v| v.as_f64()).unwrap_or(10.0);
            let rate = request.parameters.get("rate")
                .and_then(|v| v.as_f64()).unwrap_or(1.0);

            match simulate_poisson_process(time_horizon, rate) {
                Ok(events) => Ok(json!({ "events": events })),
                Err(e) => Err(e),
            }
        },
        "stochastic_integral" => {
            let params: StochasticIntegralParams = match serde_json::from_value(serde_json::Value::Object(
                request.parameters.clone().into_iter().collect()
            )) {
                Ok(p) => p,
                Err(e) => return ComputationResponse {
                    success: false,
                    module: request.module.clone(),
                    operation: request.operation.clone(),
                    result: None,
                    error: Some(format!("Invalid parameters: {}", e)),
                },
            };

            match compute_stochastic_integral_monte_carlo(params) {
                Ok(integral_value) => Ok(json!({ "integral": integral_value })),
                Err(e) => Err(e),
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
                Ok(path) => Ok(json!({ "path": path })),
                Err(e) => Err(e),
            }
        },
        _ => Err(format!("Unknown stochastic operation: {}", request.operation)),
    };

    match result {
        Ok(value) => ComputationResponse {
            success: true,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: Some(value),
            error: None,
        },
        Err(e) => ComputationResponse {
            success: false,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: None,
            error: Some(e),
        },
    }
}
