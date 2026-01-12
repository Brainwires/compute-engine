//! Financial Model Simulations
//!
//! Black-Scholes, Heston, SABR, and stochastic volatility models

use crate::engine::*;

/// Simulate finance models
pub fn simulate_finance(
    finance_model: &FinanceModel,
    input: &SimulateInput,
) -> ToolResult<SimulateOutput> {
    let range = input.range.ok_or("range [start, end] required")?;
    let steps = input.steps.unwrap_or(1000);
    let dt = (range[1] - range[0]) / steps as f64;

    match finance_model {
        FinanceModel::Heston => {
            let s0 = input.parameters.get("initial_price").unwrap_or(&100.0);
            let v0 = input.parameters.get("initial_variance").unwrap_or(&0.04);
            let kappa = input.parameters.get("kappa").unwrap_or(&2.0);
            let theta = input.parameters.get("theta").unwrap_or(&0.04);
            let sigma = input.parameters.get("sigma").unwrap_or(&0.3);
            let rho = input.parameters.get("rho").unwrap_or(&-0.7);
            let r = input.parameters.get("risk_free_rate").unwrap_or(&0.05);

            let mut times = Vec::with_capacity(steps + 1);
            let mut prices = Vec::with_capacity(steps + 1);
            let mut variances = Vec::with_capacity(steps + 1);

            let mut s = *s0;
            let mut v = *v0;
            let mut t = range[0];

            times.push(t);
            prices.push(s);
            variances.push(v);

            for _ in 0..steps {
                let dw1 = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let dw2_indep = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let dw2 = rho * dw1 + (1.0 - rho * rho).sqrt() * dw2_indep;

                s = s * (1.0 + r * dt + v.abs().sqrt() * dw1);
                v = v + kappa * (theta - v) * dt + sigma * v.abs().sqrt() * dw2;
                v = v.max(0.0);

                t += dt;
                times.push(t);
                prices.push(s);
                variances.push(v);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("price".to_string(), prices);
            results.insert("variance".to_string(), variances);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "model": "heston",
                    "kappa": kappa,
                    "theta": theta,
                    "sigma": sigma,
                    "rho": rho
                })),
            })
        }

        FinanceModel::SABR => {
            let f0 = input.parameters.get("forward_rate").unwrap_or(&0.05);
            let alpha = input.parameters.get("alpha").unwrap_or(&0.3);
            let beta = input.parameters.get("beta").unwrap_or(&0.5);
            let rho = input.parameters.get("rho").unwrap_or(&-0.3);
            let nu = input.parameters.get("nu").unwrap_or(&0.4);

            let mut times = Vec::with_capacity(steps + 1);
            let mut forward_rates = Vec::with_capacity(steps + 1);
            let mut volatilities = Vec::with_capacity(steps + 1);

            let mut f = *f0;
            let mut alpha_t = *alpha;
            let mut t = range[0];

            times.push(t);
            forward_rates.push(f);
            volatilities.push(alpha_t);

            for _ in 0..steps {
                let dw1 = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let dw2_indep = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let dw2 = rho * dw1 + (1.0 - rho * rho).sqrt() * dw2_indep;

                f = f + alpha_t * f.powf(*beta) * dw1;
                alpha_t = alpha_t * (1.0 + nu * dw2);

                t += dt;
                times.push(t);
                forward_rates.push(f);
                volatilities.push(alpha_t);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("forward_rate".to_string(), forward_rates);
            results.insert("volatility".to_string(), volatilities);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "model": "SABR",
                    "beta": beta,
                    "rho": rho,
                    "nu": nu
                })),
            })
        }

        FinanceModel::StochasticVolatility => {
            let s0 = input.parameters.get("initial_price").unwrap_or(&100.0);
            let v0 = input.parameters.get("initial_volatility").unwrap_or(&0.2);
            let mu = input.parameters.get("drift").unwrap_or(&0.05);
            let kappa = input.parameters.get("mean_reversion").unwrap_or(&1.0);
            let theta = input.parameters.get("long_term_vol").unwrap_or(&0.2);
            let sigma_v = input.parameters.get("vol_of_vol").unwrap_or(&0.3);

            let mut times = Vec::with_capacity(steps + 1);
            let mut prices = Vec::with_capacity(steps + 1);
            let mut vols = Vec::with_capacity(steps + 1);

            let mut s = *s0;
            let mut v = *v0;
            let mut t = range[0];

            times.push(t);
            prices.push(s);
            vols.push(v);

            for _ in 0..steps {
                let dw1 = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                let dw2 = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;

                s = s * (1.0 + mu * dt + v * dw1);
                v = v + kappa * (theta - v) * dt + sigma_v * v * dw2;
                v = v.max(0.01);

                t += dt;
                times.push(t);
                prices.push(s);
                vols.push(v);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("price".to_string(), prices);
            results.insert("volatility".to_string(), vols);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "model": "stochastic_volatility",
                    "mean_reversion": kappa,
                    "long_term_vol": theta
                })),
            })
        }

        FinanceModel::BlackScholes => {
            let s0 = input.parameters.get("initial_price").unwrap_or(&100.0);
            let mu = input.parameters.get("drift").unwrap_or(&0.05);
            let sigma = input.parameters.get("volatility").unwrap_or(&0.2);

            let mut times = Vec::with_capacity(steps + 1);
            let mut prices = Vec::with_capacity(steps + 1);

            let mut s = *s0;
            let mut t = range[0];

            times.push(t);
            prices.push(s);

            for _ in 0..steps {
                let dw = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0;
                s = s * (1.0 + mu * dt + sigma * dw);

                t += dt;
                times.push(t);
                prices.push(s);
            }

            let mut results = std::collections::HashMap::new();
            results.insert("price".to_string(), prices);

            Ok(SimulateOutput {
                results,
                time: Some(times),
                moments: None,
                plots: None,
                metadata: Some(serde_json::json!({
                    "model": "black_scholes",
                    "drift": mu,
                    "volatility": sigma,
                    "sde": "dS = μS dt + σS dW"
                })),
            })
        }
    }
}
