//! Unified Simulator implementation
//!
//! Routes simulation requests to ODE/PDE solvers and stochastic process generators

use crate::engine::*;

pub struct UnifiedSimulator;

impl UnifiedSimulator {
    pub fn new() -> Self {
        Self
    }

    /// Simulate time evolution using ODE/PDE solvers
    fn simulate_time_evolution(
        &self,
        method: &TimeEvolutionMethod,
        input: &SimulateInput,
    ) -> ToolResult<SimulateOutput> {
        use crate::tools::numerical_methods;

        let initial_conditions = input
            .initial_conditions
            .as_ref()
            .ok_or("initial_conditions required for time evolution")?;

        let range = input.range.ok_or("range [start, end] required")?;
        let steps = input.steps.unwrap_or(100);

        match method {
            TimeEvolutionMethod::Euler | TimeEvolutionMethod::RungeKutta4 => {
                // ODE solver
                let method_str = match method {
                    TimeEvolutionMethod::Euler => "euler",
                    TimeEvolutionMethod::RungeKutta4 => "rk4",
                    _ => "euler",
                };

                let initial_value = initial_conditions
                    .values()
                    .next()
                    .ok_or("At least one initial condition required")?;

                let step_size = (range[1] - range[0]) / steps as f64;

                let result = numerical_methods::solve_ode(numerical_methods::ODESolverRequest {
                    method: method_str.to_string(),
                    initial_value: *initial_value,
                    t_start: range[0],
                    t_end: range[1],
                    step_size,
                    derivative_expression: None,
                })
                .map_err(|e| e.to_string())?;

                let mut results = std::collections::HashMap::new();
                results.insert(
                    input.variables.first().unwrap_or(&"y".to_string()).clone(),
                    result.y_values,
                );

                Ok(SimulateOutput {
                    results,
                    time: Some(result.t_values),
                    moments: None,
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": result.method_used,
                        "steps": result.steps_taken
                    })),
                })
            }

            TimeEvolutionMethod::AdaptiveStep => {
                // Adaptive step size using error estimation
                let initial_value = initial_conditions
                    .values()
                    .next()
                    .ok_or("At least one initial condition required")?;

                let mut t = range[0];
                let mut y = *initial_value;
                let mut times = vec![t];
                let mut values = vec![y];

                let mut h = (range[1] - range[0]) / steps as f64;
                let tolerance = input.parameters.get("tolerance").unwrap_or(&1e-6);

                while t < range[1] {
                    // Simple adaptive step: if error estimate high, reduce step
                    let k1 = 0.1 * y; // Simplified derivative
                    let y_euler = y + h * k1;
                    let y_rk2 = y + h * (k1 + 0.1 * y_euler) / 2.0;

                    let error_estimate = (y_rk2 - y_euler).abs();

                    if error_estimate > *tolerance {
                        h *= 0.5; // Reduce step size
                        continue;
                    }

                    y = y_rk2;
                    t += h;
                    times.push(t);
                    values.push(y);

                    if error_estimate < *tolerance / 10.0 {
                        h *= 1.5; // Increase step size
                    }
                }

                let mut results = std::collections::HashMap::new();
                results.insert(
                    input.variables.first().unwrap_or(&"y".to_string()).clone(),
                    values,
                );

                Ok(SimulateOutput {
                    results,
                    time: Some(times),
                    moments: None,
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": "adaptive_step",
                        "tolerance": tolerance
                    })),
                })
            }

            TimeEvolutionMethod::ImplicitEuler => {
                // Implicit Euler: y_{n+1} = y_n + h*f(t_{n+1}, y_{n+1})
                let initial_value = initial_conditions
                    .values()
                    .next()
                    .ok_or("At least one initial condition required")?;

                let h = (range[1] - range[0]) / steps as f64;
                let mut times = Vec::with_capacity(steps + 1);
                let mut values = Vec::with_capacity(steps + 1);

                let mut t = range[0];
                let mut y = *initial_value;
                times.push(t);
                values.push(y);

                for _ in 0..steps {
                    // Implicit Euler with fixed-point iteration
                    let max_iter = 10;
                    let mut y_next = y;

                    for _ in 0..max_iter {
                        // f(t, y) approximation
                        let f = -0.1 * y_next; // Simplified derivative
                        let y_new = y + h * f;

                        if (y_new - y_next).abs() < 1e-8 {
                            y_next = y_new;
                            break;
                        }
                        y_next = y_new;
                    }

                    t += h;
                    y = y_next;
                    times.push(t);
                    values.push(y);
                }

                let mut results = std::collections::HashMap::new();
                results.insert(
                    input.variables.first().unwrap_or(&"y".to_string()).clone(),
                    values,
                );

                Ok(SimulateOutput {
                    results,
                    time: Some(times),
                    moments: None,
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": "implicit_euler",
                        "note": "Stable for stiff equations"
                    })),
                })
            }
        }
    }

    /// Simulate stochastic processes
    fn simulate_stochastic(
        &self,
        process: &StochasticProcess,
        input: &SimulateInput,
    ) -> ToolResult<SimulateOutput> {
        use crate::mathematics::calculus::stochastic;

        let range = input.range.ok_or("range [start, end] required")?;
        let steps = input.steps.unwrap_or(1000);
        let num_paths = input.num_paths.unwrap_or(1);

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
                // Geometric Brownian Motion: dS = μS dt + σS dW
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
                    let dw = (dt.sqrt()) * (rand::random::<f64>() - 0.5) * 2.0; // Wiener increment
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
                // Poisson process with rate λ
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
                // Lévy process (simplified as stable process)
                let alpha = input.parameters.get("alpha").unwrap_or(&1.5); // Stability parameter
                let beta = input.parameters.get("beta").unwrap_or(&0.0); // Skewness

                let dt = (range[1] - range[0]) / steps as f64;
                let mut times = Vec::with_capacity(steps + 1);
                let mut values = Vec::with_capacity(steps + 1);

                let mut x = 0.0;
                let mut t = range[0];
                times.push(t);
                values.push(x);

                for _ in 0..steps {
                    // Simplified Lévy increment
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
                // Merton jump-diffusion: dS = μS dt + σS dW + S dJ
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
                // Fractional Brownian Motion with Hurst parameter H
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
                    // Simplified fBM increment with Hurst scaling
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
                // Mean-reverting process (Cox-Ingersoll-Ross type)
                let kappa = input.parameters.get("kappa").unwrap_or(&0.5); // Mean reversion speed
                let theta = input.parameters.get("theta").unwrap_or(&1.0); // Long-term mean
                let sigma = input.parameters.get("sigma").unwrap_or(&0.2); // Volatility
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
                    // dx = κ(θ - x)dt + σ√x dW (CIR process)
                    let dx = kappa * (theta - x) * dt + sigma * x.abs().sqrt() * dw;
                    x += dx;
                    x = x.max(0.0); // Keep positive
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
                // Variance Gamma process
                let nu = input.parameters.get("nu").unwrap_or(&0.2); // Variance rate
                let theta_vg = input.parameters.get("theta").unwrap_or(&0.0); // Drift
                let sigma = input.parameters.get("sigma").unwrap_or(&0.3); // Volatility

                let dt = (range[1] - range[0]) / steps as f64;
                let mut times = Vec::with_capacity(steps + 1);
                let mut values = Vec::with_capacity(steps + 1);

                let mut x = 0.0;
                let mut t = range[0];
                times.push(t);
                values.push(x);

                for _ in 0..steps {
                    // VG = θ*γ + σ*W(γ) where γ is gamma time change
                    let gamma_increment = nu * dt; // Simplified
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

    /// Simulate fluid dynamics
    fn simulate_fluid(
        &self,
        fluid_sim: &FluidSim,
        input: &SimulateInput,
    ) -> ToolResult<SimulateOutput> {
        match fluid_sim {
            FluidSim::LatticeBotzmann => {
                // Lattice Boltzmann Method for fluid simulation
                let nx = input
                    .parameters
                    .get("nx")
                    .map(|v| *v as usize)
                    .unwrap_or(100);
                let ny = input
                    .parameters
                    .get("ny")
                    .map(|v| *v as usize)
                    .unwrap_or(100);
                let steps = input.steps.unwrap_or(1000);
                let tau = input.parameters.get("tau").unwrap_or(&0.6); // Relaxation time
                let u_lid = input.parameters.get("lid_velocity").unwrap_or(&0.1); // Lid-driven cavity velocity

                // D2Q9 model (2D, 9 velocities)
                // Initialize distribution functions
                let mut f = vec![vec![vec![0.0; 9]; ny]; nx];
                let mut f_eq = vec![vec![vec![0.0; 9]; ny]; nx];

                // Lattice velocities for D2Q9
                let cx = vec![0, 1, 0, -1, 0, 1, -1, -1, 1];
                let cy = vec![0, 0, 1, 0, -1, 1, 1, -1, -1];
                let w = vec![
                    4.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 36.0,
                    1.0 / 36.0,
                    1.0 / 36.0,
                    1.0 / 36.0,
                ];

                // Initialize with equilibrium distribution
                for i in 0..nx {
                    for j in 0..ny {
                        for k in 0..9 {
                            f[i][j][k] = w[k];
                        }
                    }
                }

                // Main LBM loop
                for _step in 0..steps {
                    // Compute macroscopic quantities (density and velocity)
                    let mut rho = vec![vec![0.0; ny]; nx];
                    let mut ux = vec![vec![0.0; ny]; nx];
                    let mut uy = vec![vec![0.0; ny]; nx];

                    for i in 0..nx {
                        for j in 0..ny {
                            for k in 0..9 {
                                rho[i][j] += f[i][j][k];
                                ux[i][j] += f[i][j][k] * cx[k] as f64;
                                uy[i][j] += f[i][j][k] * cy[k] as f64;
                            }
                            ux[i][j] /= rho[i][j];
                            uy[i][j] /= rho[i][j];
                        }
                    }

                    // Apply lid boundary condition (top boundary moves with velocity u_lid)
                    for i in 0..nx {
                        ux[i][ny - 1] = *u_lid;
                        uy[i][ny - 1] = 0.0;
                    }

                    // Compute equilibrium distribution
                    for i in 0..nx {
                        for j in 0..ny {
                            let u2 = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
                            for k in 0..9 {
                                let cu = cx[k] as f64 * ux[i][j] + cy[k] as f64 * uy[i][j];
                                f_eq[i][j][k] =
                                    rho[i][j] * w[k] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
                            }
                        }
                    }

                    // Collision step (BGK approximation)
                    for i in 0..nx {
                        for j in 0..ny {
                            for k in 0..9 {
                                f[i][j][k] += -(f[i][j][k] - f_eq[i][j][k]) / tau;
                            }
                        }
                    }

                    // Streaming step
                    let mut f_new = f.clone();
                    for i in 0..nx {
                        for j in 0..ny {
                            for k in 0..9 {
                                let i_new = ((i as i32 + cx[k]) + nx as i32) % nx as i32;
                                let j_new = ((j as i32 + cy[k]) + ny as i32) % ny as i32;
                                f_new[i_new as usize][j_new as usize][k] = f[i][j][k];
                            }
                        }
                    }
                    f = f_new;
                }

                // Extract final velocity field
                let mut final_ux = Vec::new();
                let mut final_uy = Vec::new();
                let mut final_rho = Vec::new();

                for i in 0..nx {
                    for j in 0..ny {
                        let mut rho_ij = 0.0;
                        let mut ux_ij = 0.0;
                        let mut uy_ij = 0.0;
                        for k in 0..9 {
                            rho_ij += f[i][j][k];
                            ux_ij += f[i][j][k] * cx[k] as f64;
                            uy_ij += f[i][j][k] * cy[k] as f64;
                        }
                        ux_ij /= rho_ij;
                        uy_ij /= rho_ij;

                        final_ux.push(ux_ij);
                        final_uy.push(uy_ij);
                        final_rho.push(rho_ij);
                    }
                }

                let mut results = std::collections::HashMap::new();
                results.insert("velocity_x".to_string(), final_ux);
                results.insert("velocity_y".to_string(), final_uy);
                results.insert("density".to_string(), final_rho);

                Ok(SimulateOutput {
                    results,
                    time: None,
                    moments: None,
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": "lattice_boltzmann_d2q9",
                        "grid_size": [nx, ny],
                        "steps": steps,
                        "tau": tau,
                        "lid_velocity": u_lid,
                        "note": "D2Q9 model with BGK collision operator"
                    })),
                })
            }

            FluidSim::NavierStokes | FluidSim::Euler => Err(format!(
                "Fluid simulation {:?} should use Solve tool with FluidEquation type",
                fluid_sim
            )),

            FluidSim::QuantumNavierStokes1D => {
                use crate::physics::fluid_dynamics::quantum_navier_stokes_1d::{QNS1DConfig, QNS1DSolver};

                // Extract parameters
                let nx = input
                    .parameters
                    .get("nx")
                    .map(|v| *v as usize)
                    .unwrap_or(200);
                let domain_length = *input.parameters.get("domain_length").unwrap_or(&1.0);
                let viscosity = *input.parameters.get("viscosity").unwrap_or(&0.01);
                let particle_mass = *input.parameters.get("mass").unwrap_or(&4.8e-26); // N2 molecule
                let gamma = *input.parameters.get("gamma").unwrap_or(&1.4);
                let enable_quantum = input
                    .parameters
                    .get("enable_quantum")
                    .map(|v| *v != 0.0)
                    .unwrap_or(true);
                let cfl = *input.parameters.get("cfl").unwrap_or(&0.5);

                // Time parameters
                let range = input.range.ok_or("range [start, end] required for QNS simulation")?;
                let end_time = range[1];

                // Create config
                let config = QNS1DConfig {
                    nx,
                    length: domain_length,
                    viscosity,
                    particle_mass,
                    gamma,
                    enable_quantum,
                    cfl,
                };

                // Create solver
                let mut solver = QNS1DSolver::new(config);

                // Check for initial condition type
                let init_type = input
                    .parameters
                    .get("init_type")
                    .map(|v| *v as i32)
                    .unwrap_or(0);

                if init_type == 1 {
                    // Sod shock tube (uses standard Sod parameters)
                    solver.init_sod_shock_tube();
                }

                // Run simulation until end time
                solver.run_until(end_time);

                // Extract results
                let state = solver.get_state();
                let mut results = std::collections::HashMap::new();
                results.insert("density".to_string(), state.density.clone());
                results.insert("velocity".to_string(), state.velocity.clone());
                results.insert("pressure".to_string(), state.pressure.clone());
                results.insert("x".to_string(), solver.get_x_coords());

                Ok(SimulateOutput {
                    results,
                    time: Some(vec![state.time]),
                    moments: None,
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": "quantum_navier_stokes_1d",
                        "grid_size": nx,
                        "domain_length": domain_length,
                        "viscosity": viscosity,
                        "particle_mass": particle_mass,
                        "gamma": gamma,
                        "enable_quantum": enable_quantum,
                        "cfl": cfl,
                        "steps_taken": state.step_count,
                        "note": "1D Quantum Navier-Stokes with Bohm potential correction"
                    })),
                })
            }

            FluidSim::QuantumNavierStokes2D => {
                use crate::physics::fluid_dynamics::quantum_navier_stokes_2d::{QNS2DConfig, QNS2DSolver};

                // Extract parameters
                let nx = input
                    .parameters
                    .get("nx")
                    .map(|v| *v as usize)
                    .unwrap_or(64);
                let ny = input
                    .parameters
                    .get("ny")
                    .map(|v| *v as usize)
                    .unwrap_or(64);
                let lx = *input.parameters.get("lx").unwrap_or(&(2.0 * std::f64::consts::PI));
                let ly = *input.parameters.get("ly").unwrap_or(&(2.0 * std::f64::consts::PI));
                let viscosity = *input.parameters.get("viscosity").unwrap_or(&0.01);
                let particle_mass = *input.parameters.get("mass").unwrap_or(&4.8e-26);
                let rho_ref = *input.parameters.get("rho_ref").unwrap_or(&1.0);
                let gamma = *input.parameters.get("gamma").unwrap_or(&1.4);
                let enable_quantum = input
                    .parameters
                    .get("enable_quantum")
                    .map(|v| *v != 0.0)
                    .unwrap_or(true);
                let cfl = *input.parameters.get("cfl").unwrap_or(&0.3);
                let sound_speed = *input.parameters.get("sound_speed").unwrap_or(&10.0);

                // Time parameters
                let range = input.range.ok_or("range [start, end] required for QNS simulation")?;
                let end_time = range[1];

                // Create config
                let config = QNS2DConfig {
                    nx,
                    ny,
                    lx,
                    ly,
                    viscosity,
                    particle_mass,
                    rho_ref,
                    gamma,
                    enable_quantum,
                    cfl,
                    sound_speed,
                };

                // Create solver
                let mut solver = QNS2DSolver::new(config);

                // Check for initial condition type
                let init_type = input
                    .parameters
                    .get("init_type")
                    .map(|v| *v as i32)
                    .unwrap_or(0);

                match init_type {
                    1 => {
                        // Taylor-Green vortex
                        let u0 = *input.parameters.get("u0").unwrap_or(&1.0);
                        let k = *input.parameters.get("wavenumber").unwrap_or(&1.0);
                        solver.init_taylor_green(u0, k);
                    }
                    2 => {
                        // Decaying turbulence
                        let energy = *input.parameters.get("energy").unwrap_or(&1.0);
                        let k_peak = *input.parameters.get("k_peak").unwrap_or(&4.0);
                        solver.init_decaying_turbulence(energy, k_peak);
                    }
                    _ => {
                        // Default: uniform flow (already initialized)
                    }
                }

                // Run simulation until end time
                solver.run_until(end_time);

                // Extract results - flatten 2D arrays to 1D for JSON
                let state = solver.get_state();
                let mut results = std::collections::HashMap::new();

                // Flatten density, u, v, pressure
                let density_flat: Vec<f64> = state.density.iter().cloned().collect();
                let u_flat: Vec<f64> = state.u.iter().cloned().collect();
                let v_flat: Vec<f64> = state.v.iter().cloned().collect();
                let pressure_flat: Vec<f64> = state.pressure.iter().cloned().collect();

                results.insert("density".to_string(), density_flat);
                results.insert("velocity_x".to_string(), u_flat);
                results.insert("velocity_y".to_string(), v_flat);
                results.insert("pressure".to_string(), pressure_flat);
                results.insert("x".to_string(), solver.get_x_coords());
                results.insert("y".to_string(), solver.get_y_coords());

                Ok(SimulateOutput {
                    results,
                    time: Some(vec![state.time]),
                    moments: Some(std::collections::HashMap::from([
                        ("kinetic_energy".to_string(), solver.kinetic_energy()),
                        ("enstrophy".to_string(), solver.enstrophy()),
                        ("max_vorticity".to_string(), solver.max_vorticity()),
                        ("max_divergence".to_string(), solver.max_divergence()),
                    ])),
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": "quantum_navier_stokes_2d",
                        "grid_size": [nx, ny],
                        "domain_size": [lx, ly],
                        "viscosity": viscosity,
                        "particle_mass": particle_mass,
                        "gamma": gamma,
                        "enable_quantum": enable_quantum,
                        "cfl": cfl,
                        "sound_speed": sound_speed,
                        "steps_taken": state.step_count,
                        "note": "2D Quantum Navier-Stokes with Bohm potential correction"
                    })),
                })
            }

            FluidSim::NavierStokes3D => {
                use crate::physics::fluid_dynamics::navier_stokes_3d::{NS3DConfig, NS3DSolver};

                // Extract parameters
                let nx = input
                    .parameters
                    .get("nx")
                    .map(|v| *v as usize)
                    .unwrap_or(32);
                let ny = input
                    .parameters
                    .get("ny")
                    .map(|v| *v as usize)
                    .unwrap_or(32);
                let nz = input
                    .parameters
                    .get("nz")
                    .map(|v| *v as usize)
                    .unwrap_or(32);
                let lx = *input.parameters.get("lx").unwrap_or(&(2.0 * std::f64::consts::PI));
                let ly = *input.parameters.get("ly").unwrap_or(&(2.0 * std::f64::consts::PI));
                let lz = *input.parameters.get("lz").unwrap_or(&(2.0 * std::f64::consts::PI));
                let viscosity = *input.parameters.get("viscosity").unwrap_or(&0.01);
                let rho = *input.parameters.get("rho").unwrap_or(&1.0);
                let cfl = *input.parameters.get("cfl").unwrap_or(&0.5);

                // Time parameters
                let range = input.range.ok_or("range [start, end] required for NS3D simulation")?;
                let end_time = range[1];

                // Create config
                let config = NS3DConfig {
                    nx,
                    ny,
                    nz,
                    lx,
                    ly,
                    lz,
                    viscosity,
                    rho,
                    cfl,
                    ..Default::default()
                };

                // Create solver
                let mut solver = NS3DSolver::new(config);

                // Check for initial condition type
                let init_type = input
                    .parameters
                    .get("init_type")
                    .map(|v| *v as i32)
                    .unwrap_or(0);

                match init_type {
                    1 => {
                        // Taylor-Green vortex
                        let u0 = *input.parameters.get("u0").unwrap_or(&1.0);
                        solver.init_taylor_green(u0);
                    }
                    _ => {
                        // Default: uniform flow
                        let u0 = *input.parameters.get("u0").unwrap_or(&0.0);
                        let v0 = *input.parameters.get("v0").unwrap_or(&0.0);
                        let w0 = *input.parameters.get("w0").unwrap_or(&0.0);
                        solver.init_uniform(u0, v0, w0);
                    }
                }

                // Run simulation until end time
                solver.run_until(end_time);

                // Extract results - flatten 3D arrays to 1D
                let state = solver.get_state();
                let mut results = std::collections::HashMap::new();

                let u_flat: Vec<f64> = state.u.iter().cloned().collect();
                let v_flat: Vec<f64> = state.v.iter().cloned().collect();
                let w_flat: Vec<f64> = state.w.iter().cloned().collect();
                let p_flat: Vec<f64> = state.p.iter().cloned().collect();

                results.insert("velocity_x".to_string(), u_flat);
                results.insert("velocity_y".to_string(), v_flat);
                results.insert("velocity_z".to_string(), w_flat);
                results.insert("pressure".to_string(), p_flat);
                results.insert("x".to_string(), solver.get_x_coords());
                results.insert("y".to_string(), solver.get_y_coords());
                results.insert("z".to_string(), solver.get_z_coords());

                Ok(SimulateOutput {
                    results,
                    time: Some(vec![state.time]),
                    moments: Some(std::collections::HashMap::from([
                        ("kinetic_energy".to_string(), solver.kinetic_energy()),
                        ("enstrophy".to_string(), solver.enstrophy()),
                        ("max_vorticity".to_string(), solver.max_vorticity()),
                        ("max_divergence".to_string(), solver.max_divergence()),
                        ("reynolds_number".to_string(), solver.reynolds_number()),
                    ])),
                    plots: None,
                    metadata: Some(serde_json::json!({
                        "method": "navier_stokes_3d",
                        "grid_size": [nx, ny, nz],
                        "domain_size": [lx, ly, lz],
                        "viscosity": viscosity,
                        "rho": rho,
                        "cfl": cfl,
                        "steps_taken": state.step_count,
                        "note": "3D incompressible Navier-Stokes with fractional step method"
                    })),
                })
            }
        }
    }

    /// Simulate finance models
    fn simulate_finance(
        &self,
        finance_model: &FinanceModel,
        input: &SimulateInput,
    ) -> ToolResult<SimulateOutput> {
        let range = input.range.ok_or("range [start, end] required")?;
        let steps = input.steps.unwrap_or(1000);
        let dt = (range[1] - range[0]) / steps as f64;

        match finance_model {
            FinanceModel::Heston => {
                // Heston stochastic volatility model
                let s0 = input.parameters.get("initial_price").unwrap_or(&100.0);
                let v0 = input.parameters.get("initial_variance").unwrap_or(&0.04);
                let kappa = input.parameters.get("kappa").unwrap_or(&2.0); // Mean reversion
                let theta = input.parameters.get("theta").unwrap_or(&0.04); // Long-term variance
                let sigma = input.parameters.get("sigma").unwrap_or(&0.3); // Vol of vol
                let rho = input.parameters.get("rho").unwrap_or(&-0.7); // Correlation
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

                    // dS = rS dt + √v S dW1
                    // dv = κ(θ - v)dt + σ√v dW2
                    s = s * (1.0 + r * dt + v.abs().sqrt() * dw1);
                    v = v + kappa * (theta - v) * dt + sigma * v.abs().sqrt() * dw2;
                    v = v.max(0.0); // Keep variance positive

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
                // SABR (Stochastic Alpha Beta Rho) model
                let f0 = input.parameters.get("forward_rate").unwrap_or(&0.05);
                let alpha = input.parameters.get("alpha").unwrap_or(&0.3); // Initial vol
                let beta = input.parameters.get("beta").unwrap_or(&0.5); // CEV exponent
                let rho = input.parameters.get("rho").unwrap_or(&-0.3); // Correlation
                let nu = input.parameters.get("nu").unwrap_or(&0.4); // Vol of vol

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

                    // dF = α F^β dW1
                    // dα = ν α dW2
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
                // General stochastic volatility model (similar to Heston but simplified)
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
                    v = v.max(0.01); // Keep volatility positive

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
                // Black-Scholes model (Geometric Brownian Motion for asset price)
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
                    // dS = μS dt + σS dW
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
}

impl Simulate for UnifiedSimulator {
    fn simulate(&self, input: &SimulateInput) -> ToolResult<SimulateOutput> {
        match &input.model {
            SimulationModel::TimeEvolution(method) => self.simulate_time_evolution(method, input),

            SimulationModel::Stochastic(process) => self.simulate_stochastic(process, input),

            SimulationModel::FluidDynamics(fluid_sim) => self.simulate_fluid(fluid_sim, input),

            SimulationModel::Finance(finance_model) => self.simulate_finance(finance_model, input),
        }
    }
}
