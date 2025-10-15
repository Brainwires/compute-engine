use ndarray::Array2;
use serde_json::{json, Value};
use std::error::Error;

pub struct NavierStokesSolver {
    config: Value,
}

impl NavierStokesSolver {
    pub fn new(config: Value) -> Self {
        Self { config }
    }

    pub fn solve(&self) -> Result<Value, Box<dyn Error>> {
        // Extract configuration parameters
        let domain = self.config["domain"].as_object()
            .ok_or("Missing domain configuration")?;
        let fluid_props = self.config["fluid_properties"].as_object()
            .ok_or("Missing fluid properties")?;
        let boundary_conds = self.config["boundary_conditions"].as_object()
            .ok_or("Missing boundary conditions")?;
        let initial_conds = self.config["initial_conditions"].as_object()
            .ok_or("Missing initial conditions")?;
        let sim_params = self.config["simulation_parameters"].as_object()
            .ok_or("Missing simulation parameters")?;

        // Parse domain parameters
        let width = domain["width"].as_f64().ok_or("Invalid domain width")?;
        let height = domain["height"].as_f64().ok_or("Invalid domain height")?;
        let nx = domain["grid_points_x"].as_u64().ok_or("Invalid grid points x")? as usize;
        let ny = domain["grid_points_y"].as_u64().ok_or("Invalid grid points y")? as usize;

        // Parse fluid properties
        let viscosity = fluid_props["viscosity"].as_f64().ok_or("Invalid viscosity")?;
        let density = fluid_props["density"].as_f64().ok_or("Invalid density")?;

        // Parse simulation parameters
        let dt = sim_params["time_step"].as_f64().ok_or("Invalid time step")?;
        let total_time = sim_params["total_time"].as_f64().ok_or("Invalid total time")?;
        let time_steps = (total_time / dt) as usize;

        // Initialize grid
        let dx = width / (nx - 1) as f64;
        let dy = height / (ny - 1) as f64;

        // Initialize velocity and pressure fields
        let mut u = Array2::<f64>::zeros((nx, ny));
        let mut v = Array2::<f64>::zeros((nx, ny));
        let mut p = Array2::<f64>::zeros((nx, ny));

        // Set initial conditions
        let initial_u = initial_conds["velocity_x"].as_f64().unwrap_or(0.0);
        let initial_v = initial_conds["velocity_y"].as_f64().unwrap_or(0.0);
        let initial_p = initial_conds["pressure"].as_f64().unwrap_or(0.0);

        u.fill(initial_u);
        v.fill(initial_v);
        p.fill(initial_p);

        // Apply boundary conditions
        self.apply_boundary_conditions(&mut u, &mut v, boundary_conds, nx, ny)?;

        // Time stepping loop
        for step in 0..time_steps {
            // Solve using fractional step method
            self.solve_momentum_step(&mut u, &mut v, &p, dt, dx, dy, viscosity, density)?;
            self.apply_boundary_conditions(&mut u, &mut v, boundary_conds, nx, ny)?;
            self.solve_pressure_step(&mut u, &mut v, &mut p, dt, dx, dy)?;
            self.correct_velocity(&mut u, &mut v, &p, dt, dx, dy, density);

            // Progress reporting
            if step % (time_steps / 10).max(1) == 0 {
                eprintln!("2D Navier-Stokes simulation: {:.1}% complete", 
                         100.0 * step as f64 / time_steps as f64);
            }
        }

        // Calculate additional flow properties
        let vorticity = self.calculate_vorticity(&u, &v, dx, dy);
        let divergence = self.calculate_divergence(&u, &v, dx, dy);

        // Convert ndarray to nested vectors for JSON serialization
        let u_vec: Vec<Vec<f64>> = (0..nx).map(|i| (0..ny).map(|j| u[[i, j]]).collect()).collect();
        let v_vec: Vec<Vec<f64>> = (0..nx).map(|i| (0..ny).map(|j| v[[i, j]]).collect()).collect();
        let p_vec: Vec<Vec<f64>> = (0..nx).map(|i| (0..ny).map(|j| p[[i, j]]).collect()).collect();
        let vort_vec: Vec<Vec<f64>> = (0..nx).map(|i| (0..ny).map(|j| vorticity[[i, j]]).collect()).collect();
        let div_vec: Vec<Vec<f64>> = (0..nx).map(|i| (0..ny).map(|j| divergence[[i, j]]).collect()).collect();

        Ok(json!({
            "solution": {
                "velocity_field": {
                    "u": u_vec,
                    "v": v_vec
                },
                "pressure_field": p_vec,
                "derived_quantities": {
                    "vorticity": vort_vec,
                    "divergence": div_vec
                }
            },
            "simulation_parameters": {
                "domain": { "width": width, "height": height },
                "grid": { "nx": nx, "ny": ny, "dx": dx, "dy": dy },
                "time": { "dt": dt, "total_time": total_time, "time_steps": time_steps },
                "fluid": { "viscosity": viscosity, "density": density }
            }
        }))
    }

    fn apply_boundary_conditions(&self, u: &mut Array2<f64>, v: &mut Array2<f64>, 
                                boundary_conds: &serde_json::Map<String, Value>,
                                nx: usize, ny: usize) -> Result<(), Box<dyn Error>> {
        // Apply boundary conditions based on configuration
        let top = boundary_conds["top"].as_str().unwrap_or("no_slip");
        let bottom = boundary_conds["bottom"].as_str().unwrap_or("no_slip");  
        let left = boundary_conds["left"].as_str().unwrap_or("no_slip");
        let right = boundary_conds["right"].as_str().unwrap_or("no_slip");

        // Bottom boundary (j = 0)
        for i in 0..nx {
            match bottom {
                "no_slip" => { u[[i, 0]] = 0.0; v[[i, 0]] = 0.0; },
                "free_slip" => { v[[i, 0]] = 0.0; }, // ∂u/∂n = 0 handled in solver
                "pressure_outlet" => { /* handled in pressure solver */ },
                _ => return Err("Unknown bottom boundary condition".into()),
            }
        }

        // Top boundary (j = ny-1)
        for i in 0..nx {
            match top {
                "no_slip" => { u[[i, ny-1]] = 0.0; v[[i, ny-1]] = 0.0; },
                "free_slip" => { v[[i, ny-1]] = 0.0; },
                "pressure_outlet" => { /* handled in pressure solver */ },
                _ => return Err("Unknown top boundary condition".into()),
            }
        }

        // Left boundary (i = 0)
        for j in 0..ny {
            match left {
                "no_slip" => { u[[0, j]] = 0.0; v[[0, j]] = 0.0; },
                "free_slip" => { u[[0, j]] = 0.0; },
                "velocity_inlet" => { u[[0, j]] = 1.0; v[[0, j]] = 0.0; }, // default inlet velocity
                _ => return Err("Unknown left boundary condition".into()),
            }
        }

        // Right boundary (i = nx-1)  
        for j in 0..ny {
            match right {
                "no_slip" => { u[[nx-1, j]] = 0.0; v[[nx-1, j]] = 0.0; },
                "free_slip" => { u[[nx-1, j]] = 0.0; },
                "pressure_outlet" => { /* handled in pressure solver */ },
                _ => return Err("Unknown right boundary condition".into()),
            }
        }

        Ok(())
    }

    fn solve_momentum_step(&self, u: &mut Array2<f64>, v: &mut Array2<f64>, p: &Array2<f64>,
                          dt: f64, dx: f64, dy: f64, nu: f64, rho: f64) -> Result<(), Box<dyn Error>> {
        let nx = u.shape()[0];
        let ny = u.shape()[1];
        let mut u_new = u.clone();
        let mut v_new = v.clone();

        // Solve u-momentum equation
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                // Convective terms (upwind scheme)
                let u_dudx = if u[[i,j]] > 0.0 {
                    u[[i,j]] * (u[[i,j]] - u[[i-1,j]]) / dx
                } else {
                    u[[i,j]] * (u[[i+1,j]] - u[[i,j]]) / dx
                };

                let v_dudy = if v[[i,j]] > 0.0 {
                    v[[i,j]] * (u[[i,j]] - u[[i,j-1]]) / dy
                } else {
                    v[[i,j]] * (u[[i,j+1]] - u[[i,j]]) / dy
                };

                // Pressure gradient
                let dpdx = (p[[i+1,j]] - p[[i-1,j]]) / (2.0 * dx);

                // Viscous terms
                let d2udx2 = (u[[i+1,j]] - 2.0*u[[i,j]] + u[[i-1,j]]) / (dx*dx);
                let d2udy2 = (u[[i,j+1]] - 2.0*u[[i,j]] + u[[i,j-1]]) / (dy*dy);

                u_new[[i,j]] = u[[i,j]] + dt * (-u_dudx - v_dudy - dpdx/rho + nu * (d2udx2 + d2udy2));
            }
        }

        // Solve v-momentum equation
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                // Convective terms (upwind scheme)
                let u_dvdx = if u[[i,j]] > 0.0 {
                    u[[i,j]] * (v[[i,j]] - v[[i-1,j]]) / dx
                } else {
                    u[[i,j]] * (v[[i+1,j]] - v[[i,j]]) / dx
                };

                let v_dvdy = if v[[i,j]] > 0.0 {
                    v[[i,j]] * (v[[i,j]] - v[[i,j-1]]) / dy
                } else {
                    v[[i,j]] * (v[[i,j+1]] - v[[i,j]]) / dy
                };

                // Pressure gradient
                let dpdy = (p[[i,j+1]] - p[[i,j-1]]) / (2.0 * dy);

                // Viscous terms
                let d2vdx2 = (v[[i+1,j]] - 2.0*v[[i,j]] + v[[i-1,j]]) / (dx*dx);
                let d2vdy2 = (v[[i,j+1]] - 2.0*v[[i,j]] + v[[i,j-1]]) / (dy*dy);

                v_new[[i,j]] = v[[i,j]] + dt * (-u_dvdx - v_dvdy - dpdy/rho + nu * (d2vdx2 + d2vdy2));
            }
        }

        *u = u_new;
        *v = v_new;
        Ok(())
    }

    fn solve_pressure_step(&self, u: &Array2<f64>, v: &Array2<f64>, p: &mut Array2<f64>,
                          dt: f64, dx: f64, dy: f64) -> Result<(), Box<dyn Error>> {
        let nx = u.shape()[0];
        let ny = u.shape()[1];
        
        // Calculate RHS for pressure Poisson equation
        let mut rhs = Array2::<f64>::zeros((nx, ny));
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let dudx = (u[[i+1,j]] - u[[i-1,j]]) / (2.0 * dx);
                let dvdy = (v[[i,j+1]] - v[[i,j-1]]) / (2.0 * dy);
                rhs[[i,j]] = (dudx + dvdy) / dt;
            }
        }

        // Solve Poisson equation using Jacobi iterations
        let mut p_old = p.clone();
        let max_iterations = 1000;
        let tolerance = 1e-6;

        for _iter in 0..max_iterations {
            let mut max_residual = 0.0_f64;

            for i in 1..nx-1 {
                for j in 1..ny-1 {
                    let p_new = 0.25 * (
                        p_old[[i+1,j]] + p_old[[i-1,j]] + 
                        p_old[[i,j+1]] + p_old[[i,j-1]] - 
                        dx*dx * rhs[[i,j]]
                    );
                    
                    let residual = (p_new - p_old[[i,j]]).abs();
                    max_residual = max_residual.max(residual);
                    
                    p[[i,j]] = p_new;
                }
            }

            // Apply pressure boundary conditions (zero gradient)
            for i in 0..nx {
                p[[i, 0]] = p[[i, 1]];
                p[[i, ny-1]] = p[[i, ny-2]];
            }
            for j in 0..ny {
                p[[0, j]] = p[[1, j]];
                p[[nx-1, j]] = p[[nx-2, j]];
            }

            if max_residual < tolerance {
                break;
            }

            p_old = p.clone();
        }

        Ok(())
    }

    fn correct_velocity(&self, u: &mut Array2<f64>, v: &mut Array2<f64>, p: &Array2<f64>,
                       dt: f64, dx: f64, dy: f64, rho: f64) {
        let nx = u.shape()[0];
        let ny = u.shape()[1];

        // Correct velocities using pressure gradients
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let dpdx = (p[[i+1,j]] - p[[i-1,j]]) / (2.0 * dx);
                let dpdy = (p[[i,j+1]] - p[[i,j-1]]) / (2.0 * dy);

                u[[i,j]] -= dt * dpdx / rho;
                v[[i,j]] -= dt * dpdy / rho;
            }
        }
    }

    fn calculate_vorticity(&self, u: &Array2<f64>, v: &Array2<f64>, dx: f64, dy: f64) -> Array2<f64> {
        let nx = u.shape()[0];
        let ny = u.shape()[1];
        let mut vorticity = Array2::<f64>::zeros((nx, ny));

        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let dvdx = (v[[i+1,j]] - v[[i-1,j]]) / (2.0 * dx);
                let dudy = (u[[i,j+1]] - u[[i,j-1]]) / (2.0 * dy);
                vorticity[[i,j]] = dvdx - dudy;
            }
        }

        vorticity
    }

    fn calculate_divergence(&self, u: &Array2<f64>, v: &Array2<f64>, dx: f64, dy: f64) -> Array2<f64> {
        let nx = u.shape()[0];
        let ny = u.shape()[1];
        let mut divergence = Array2::<f64>::zeros((nx, ny));

        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let dudx = (u[[i+1,j]] - u[[i-1,j]]) / (2.0 * dx);
                let dvdy = (v[[i,j+1]] - v[[i,j-1]]) / (2.0 * dy);
                divergence[[i,j]] = dudx + dvdy;
            }
        }

        divergence
    }
}