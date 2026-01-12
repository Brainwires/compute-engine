//! Fluid Dynamics Module
//!
//! Provides computational fluid dynamics capabilities:
//! - Navier-Stokes equation solvers (2D classical, 1D/2D quantum, 3D)
//! - Quantum Navier-Stokes with Bohm potential corrections
//! - Analytical flow solutions
//! - Boundary condition handling
//! - Cavity flow simulations
//! - Flow field analysis
//! - Spectral analysis (energy spectrum E(k))
//! - Lattice Boltzmann method

pub mod analytical;
pub mod boundary_conditions;
pub mod cavity_flow;
pub mod flow_analysis;
pub mod grid;
pub mod navier_stokes;
pub mod navier_stokes_3d;
pub mod quantum_navier_stokes_1d;
pub mod quantum_navier_stokes_2d;
pub mod spectral_analysis;

// Re-export main types and functions
pub use analytical::*;
pub use boundary_conditions::*;
pub use cavity_flow::*;
pub use flow_analysis::*;
pub use grid::*;
pub use navier_stokes::*;
pub use navier_stokes_3d::*;
pub use quantum_navier_stokes_1d::*;
pub use quantum_navier_stokes_2d::*;
pub use spectral_analysis::*;

use crate::engine::*;

/// Simulate fluid dynamics
pub fn simulate_fluid(fluid_sim: &FluidSim, input: &SimulateInput) -> ToolResult<SimulateOutput> {
    match fluid_sim {
        FluidSim::LatticeBotzmann => simulate_lattice_boltzmann(input),
        FluidSim::NavierStokes | FluidSim::Euler => Err(format!(
            "Fluid simulation {:?} should use Solve tool with FluidEquation type",
            fluid_sim
        )),
        FluidSim::QuantumNavierStokes1D => simulate_qns_1d(input),
        FluidSim::QuantumNavierStokes2D => simulate_qns_2d(input),
        FluidSim::NavierStokes3D => simulate_ns_3d(input),
    }
}

fn simulate_lattice_boltzmann(input: &SimulateInput) -> ToolResult<SimulateOutput> {
    let nx = input.parameters.get("nx").map(|v| *v as usize).unwrap_or(100);
    let ny = input.parameters.get("ny").map(|v| *v as usize).unwrap_or(100);
    let steps = input.steps.unwrap_or(1000);
    let tau = input.parameters.get("tau").unwrap_or(&0.6);
    let u_lid = input.parameters.get("lid_velocity").unwrap_or(&0.1);

    // D2Q9 model
    let mut f = vec![vec![vec![0.0; 9]; ny]; nx];
    let mut f_eq = vec![vec![vec![0.0; 9]; ny]; nx];

    let cx = vec![0, 1, 0, -1, 0, 1, -1, -1, 1];
    let cy = vec![0, 0, 1, 0, -1, 1, 1, -1, -1];
    let w = vec![
        4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
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

        // Lid boundary condition
        for i in 0..nx {
            ux[i][ny - 1] = *u_lid;
            uy[i][ny - 1] = 0.0;
        }

        // Equilibrium distribution
        for i in 0..nx {
            for j in 0..ny {
                let u2 = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
                for k in 0..9 {
                    let cu = cx[k] as f64 * ux[i][j] + cy[k] as f64 * uy[i][j];
                    f_eq[i][j][k] = rho[i][j] * w[k] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
                }
            }
        }

        // Collision step (BGK)
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

fn simulate_qns_1d(input: &SimulateInput) -> ToolResult<SimulateOutput> {
    use quantum_navier_stokes_1d::{QNS1DConfig, QNS1DSolver};

    let nx = input.parameters.get("nx").map(|v| *v as usize).unwrap_or(200);
    let domain_length = *input.parameters.get("domain_length").unwrap_or(&1.0);
    let viscosity = *input.parameters.get("viscosity").unwrap_or(&0.01);
    let particle_mass = *input.parameters.get("mass").unwrap_or(&4.8e-26);
    let gamma = *input.parameters.get("gamma").unwrap_or(&1.4);
    let enable_quantum = input.parameters.get("enable_quantum").map(|v| *v != 0.0).unwrap_or(true);
    let cfl = *input.parameters.get("cfl").unwrap_or(&0.5);

    let range = input.range.ok_or("range [start, end] required for QNS simulation")?;
    let end_time = range[1];

    let config = QNS1DConfig {
        nx,
        length: domain_length,
        viscosity,
        particle_mass,
        gamma,
        enable_quantum,
        cfl,
    };

    let mut solver = QNS1DSolver::new(config);

    let init_type = input.parameters.get("init_type").map(|v| *v as i32).unwrap_or(0);
    if init_type == 1 {
        solver.init_sod_shock_tube();
    }

    solver.run_until(end_time);

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

fn simulate_qns_2d(input: &SimulateInput) -> ToolResult<SimulateOutput> {
    use quantum_navier_stokes_2d::{QNS2DConfig, QNS2DSolver};

    let nx = input.parameters.get("nx").map(|v| *v as usize).unwrap_or(64);
    let ny = input.parameters.get("ny").map(|v| *v as usize).unwrap_or(64);
    let lx = *input.parameters.get("lx").unwrap_or(&(2.0 * std::f64::consts::PI));
    let ly = *input.parameters.get("ly").unwrap_or(&(2.0 * std::f64::consts::PI));
    let viscosity = *input.parameters.get("viscosity").unwrap_or(&0.01);
    let particle_mass = *input.parameters.get("mass").unwrap_or(&4.8e-26);
    let rho_ref = *input.parameters.get("rho_ref").unwrap_or(&1.0);
    let gamma = *input.parameters.get("gamma").unwrap_or(&1.4);
    let enable_quantum = input.parameters.get("enable_quantum").map(|v| *v != 0.0).unwrap_or(true);
    let cfl = *input.parameters.get("cfl").unwrap_or(&0.3);
    let sound_speed = *input.parameters.get("sound_speed").unwrap_or(&10.0);

    let range = input.range.ok_or("range [start, end] required for QNS simulation")?;
    let end_time = range[1];

    let config = QNS2DConfig {
        nx, ny, lx, ly, viscosity, particle_mass, rho_ref, gamma, enable_quantum, cfl, sound_speed,
    };

    let mut solver = QNS2DSolver::new(config);

    let init_type = input.parameters.get("init_type").map(|v| *v as i32).unwrap_or(0);
    match init_type {
        1 => {
            let u0 = *input.parameters.get("u0").unwrap_or(&1.0);
            let k = *input.parameters.get("wavenumber").unwrap_or(&1.0);
            solver.init_taylor_green(u0, k);
        }
        2 => {
            let energy = *input.parameters.get("energy").unwrap_or(&1.0);
            let k_peak = *input.parameters.get("k_peak").unwrap_or(&4.0);
            solver.init_decaying_turbulence(energy, k_peak);
        }
        _ => {}
    }

    solver.run_until(end_time);

    let state = solver.get_state();
    let mut results = std::collections::HashMap::new();

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

fn simulate_ns_3d(input: &SimulateInput) -> ToolResult<SimulateOutput> {
    use navier_stokes_3d::{NS3DConfig, NS3DSolver};

    let nx = input.parameters.get("nx").map(|v| *v as usize).unwrap_or(32);
    let ny = input.parameters.get("ny").map(|v| *v as usize).unwrap_or(32);
    let nz = input.parameters.get("nz").map(|v| *v as usize).unwrap_or(32);
    let lx = *input.parameters.get("lx").unwrap_or(&(2.0 * std::f64::consts::PI));
    let ly = *input.parameters.get("ly").unwrap_or(&(2.0 * std::f64::consts::PI));
    let lz = *input.parameters.get("lz").unwrap_or(&(2.0 * std::f64::consts::PI));
    let viscosity = *input.parameters.get("viscosity").unwrap_or(&0.01);
    let rho = *input.parameters.get("rho").unwrap_or(&1.0);
    let cfl = *input.parameters.get("cfl").unwrap_or(&0.5);

    let range = input.range.ok_or("range [start, end] required for NS3D simulation")?;
    let end_time = range[1];

    let config = NS3DConfig {
        nx, ny, nz, lx, ly, lz, viscosity, rho, cfl,
        ..Default::default()
    };

    let mut solver = NS3DSolver::new(config);

    let init_type = input.parameters.get("init_type").map(|v| *v as i32).unwrap_or(0);
    match init_type {
        1 => {
            let u0 = *input.parameters.get("u0").unwrap_or(&1.0);
            solver.init_taylor_green(u0);
        }
        _ => {
            let u0 = *input.parameters.get("u0").unwrap_or(&0.0);
            let v0 = *input.parameters.get("v0").unwrap_or(&0.0);
            let w0 = *input.parameters.get("w0").unwrap_or(&0.0);
            solver.init_uniform(u0, v0, w0);
        }
    }

    solver.run_until(end_time);

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
