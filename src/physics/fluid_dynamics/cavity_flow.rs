use ndarray::Array2;
use serde_json::{json, Value};
use std::error::Error;

pub struct CavityFlowSolver {
    cavity_size: f64,
    grid_resolution: usize,
    reynolds: f64,
    lid_velocity: f64,
    dx: f64,
    dy: f64,
}

impl CavityFlowSolver {
    pub fn new(cavity_size: f64, grid_resolution: usize, reynolds: f64, lid_velocity: f64) -> Self {
        let dx = cavity_size / (grid_resolution - 1) as f64;
        let dy = dx;
        
        Self {
            cavity_size,
            grid_resolution,
            reynolds,
            lid_velocity,
            dx,
            dy,
        }
    }

    pub fn solve(&self, time_steps: usize, dt: f64) -> Result<Value, Box<dyn Error>> {
        let n = self.grid_resolution;
        
        // Initialize velocity and pressure fields
        let mut u = Array2::<f64>::zeros((n, n)); // x-velocity
        let mut v = Array2::<f64>::zeros((n, n)); // y-velocity  
        let mut p = Array2::<f64>::zeros((n, n)); // pressure
        
        // Temporary arrays for calculations
        let mut u_new = Array2::<f64>::zeros((n, n));
        let mut v_new = Array2::<f64>::zeros((n, n));
        let mut p_new = Array2::<f64>::zeros((n, n));

        // Set boundary conditions
        self.apply_boundary_conditions(&mut u, &mut v);

        // Time stepping loop
        for step in 0..time_steps {
            // Solve momentum equations (simplified explicit scheme)
            self.solve_momentum(&u, &v, &p, &mut u_new, &mut v_new, dt)?;
            
            // Apply boundary conditions
            self.apply_boundary_conditions(&mut u_new, &mut v_new);
            
            // Solve pressure Poisson equation
            self.solve_pressure(&u_new, &v_new, &mut p_new, dt)?;
            
            // Correct velocities
            self.correct_velocities(&mut u_new, &mut v_new, &p_new, dt);
            
            // Update arrays
            u = u_new.clone();
            v = v_new.clone();
            p = p_new.clone();

            // Progress reporting (every 10% of simulation)
            if step % (time_steps / 10).max(1) == 0 {
                eprintln!("Cavity flow simulation: {:.1}% complete", 
                         100.0 * step as f64 / time_steps as f64);
            }
        }

        // Calculate derived quantities
        let vorticity = self.calculate_vorticity(&u, &v);
        let stream_function = self.calculate_stream_function(&u, &v)?;
        let velocity_magnitude = self.calculate_velocity_magnitude(&u, &v);

        // Convert ndarray to nested vectors for JSON serialization
        let u_vec: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| u[[i, j]]).collect()).collect();
        let v_vec: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| v[[i, j]]).collect()).collect();
        let p_vec: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| p[[i, j]]).collect()).collect();
        let vort_vec: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| vorticity[[i, j]]).collect()).collect();
        let stream_vec: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| stream_function[[i, j]]).collect()).collect();
        let vel_mag_vec: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| velocity_magnitude[[i, j]]).collect()).collect();

        Ok(json!({
            "cavity_size": self.cavity_size,
            "grid_resolution": self.grid_resolution,
            "reynolds_number": self.reynolds,
            "lid_velocity": self.lid_velocity,
            "time_steps": time_steps,
            "dt": dt,
            "velocity_field": {
                "u": u_vec,
                "v": v_vec
            },
            "pressure_field": p_vec,
            "derived_quantities": {
                "vorticity": vort_vec,
                "stream_function": stream_vec,
                "velocity_magnitude": vel_mag_vec
            },
            "grid_parameters": {
                "dx": self.dx,
                "dy": self.dy,
                "x_coords": (0..n).map(|i| i as f64 * self.dx).collect::<Vec<_>>(),
                "y_coords": (0..n).map(|i| i as f64 * self.dy).collect::<Vec<_>>()
            }
        }))
    }

    fn apply_boundary_conditions(&self, u: &mut Array2<f64>, v: &mut Array2<f64>) {
        let n = self.grid_resolution;
        
        // No-slip boundary conditions on walls
        // Bottom wall (y = 0)
        for i in 0..n {
            u[[i, 0]] = 0.0;
            v[[i, 0]] = 0.0;
        }
        
        // Left wall (x = 0)  
        for j in 0..n {
            u[[0, j]] = 0.0;
            v[[0, j]] = 0.0;
        }
        
        // Right wall (x = L)
        for j in 0..n {
            u[[n-1, j]] = 0.0;
            v[[n-1, j]] = 0.0;
        }
        
        // Top wall (moving lid at y = L)
        for i in 0..n {
            u[[i, n-1]] = self.lid_velocity;
            v[[i, n-1]] = 0.0;
        }
    }

    fn solve_momentum(&self, u: &Array2<f64>, v: &Array2<f64>, p: &Array2<f64>, 
                     u_new: &mut Array2<f64>, v_new: &mut Array2<f64>, dt: f64) -> Result<(), Box<dyn Error>> {
        let n = self.grid_resolution;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let nu = self.lid_velocity * self.cavity_size / self.reynolds; // kinematic viscosity

        // Solve u-momentum equation
        for i in 1..n-1 {
            for j in 1..n-1 {
                // Convective terms
                let u_dudx = u[[i,j]] * (u[[i+1,j]] - u[[i-1,j]]) / (2.0 * self.dx);
                let v_dudy = v[[i,j]] * (u[[i,j+1]] - u[[i,j-1]]) / (2.0 * self.dy);
                
                // Pressure gradient
                let dpdx = (p[[i+1,j]] - p[[i-1,j]]) / (2.0 * self.dx);
                
                // Viscous terms (Laplacian)
                let d2udx2 = (u[[i+1,j]] - 2.0*u[[i,j]] + u[[i-1,j]]) / dx2;
                let d2udy2 = (u[[i,j+1]] - 2.0*u[[i,j]] + u[[i,j-1]]) / dy2;
                
                u_new[[i,j]] = u[[i,j]] + dt * (-u_dudx - v_dudy - dpdx + nu * (d2udx2 + d2udy2));
            }
        }

        // Solve v-momentum equation  
        for i in 1..n-1 {
            for j in 1..n-1 {
                // Convective terms
                let u_dvdx = u[[i,j]] * (v[[i+1,j]] - v[[i-1,j]]) / (2.0 * self.dx);
                let v_dvdy = v[[i,j]] * (v[[i,j+1]] - v[[i,j-1]]) / (2.0 * self.dy);
                
                // Pressure gradient
                let dpdy = (p[[i,j+1]] - p[[i,j-1]]) / (2.0 * self.dy);
                
                // Viscous terms (Laplacian)
                let d2vdx2 = (v[[i+1,j]] - 2.0*v[[i,j]] + v[[i-1,j]]) / dx2;
                let d2vdy2 = (v[[i,j+1]] - 2.0*v[[i,j]] + v[[i,j-1]]) / dy2;
                
                v_new[[i,j]] = v[[i,j]] + dt * (-u_dvdx - v_dvdy - dpdy + nu * (d2vdx2 + d2vdy2));
            }
        }

        Ok(())
    }

    fn solve_pressure(&self, u: &Array2<f64>, v: &Array2<f64>, p: &mut Array2<f64>, dt: f64) -> Result<(), Box<dyn Error>> {
        let n = self.grid_resolution;
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        
        // Calculate divergence of velocity field
        let mut div = Array2::<f64>::zeros((n, n));
        for i in 1..n-1 {
            for j in 1..n-1 {
                let dudx = (u[[i+1,j]] - u[[i-1,j]]) / (2.0 * self.dx);
                let dvdy = (v[[i,j+1]] - v[[i,j-1]]) / (2.0 * self.dy);
                div[[i,j]] = (dudx + dvdy) / dt;
            }
        }

        // Solve Poisson equation for pressure using Jacobi iterations
        let mut p_old = p.clone();
        let max_iterations = 1000;
        let tolerance = 1e-6;
        
        for iter in 0..max_iterations {
            let mut max_residual = 0.0_f64;
            
            for i in 1..n-1 {
                for j in 1..n-1 {
                    let p_new_val = 0.25 * (
                        p_old[[i+1,j]] + p_old[[i-1,j]] + 
                        p_old[[i,j+1]] + p_old[[i,j-1]] - 
                        dx2 * div[[i,j]]
                    );
                    
                    let residual = (p_new_val - p_old[[i,j]]).abs();
                    max_residual = max_residual.max(residual);
                    
                    p[[i,j]] = p_new_val;
                }
            }
            
            // Pressure boundary conditions (Neumann - zero gradient)
            for i in 0..n {
                p[[i, 0]] = p[[i, 1]];      // bottom
                p[[i, n-1]] = p[[i, n-2]];  // top
            }
            for j in 0..n {
                p[[0, j]] = p[[1, j]];      // left
                p[[n-1, j]] = p[[n-2, j]];  // right
            }
            
            if max_residual < tolerance {
                break;
            }
            
            p_old = p.clone();
        }

        Ok(())
    }

    fn correct_velocities(&self, u: &mut Array2<f64>, v: &mut Array2<f64>, p: &Array2<f64>, dt: f64) {
        let n = self.grid_resolution;
        
        // Correct velocities using pressure gradients
        for i in 1..n-1 {
            for j in 1..n-1 {
                let dpdx = (p[[i+1,j]] - p[[i-1,j]]) / (2.0 * self.dx);
                let dpdy = (p[[i,j+1]] - p[[i,j-1]]) / (2.0 * self.dy);
                
                u[[i,j]] -= dt * dpdx;
                v[[i,j]] -= dt * dpdy;
            }
        }
    }

    fn calculate_vorticity(&self, u: &Array2<f64>, v: &Array2<f64>) -> Array2<f64> {
        let n = self.grid_resolution;
        let mut vorticity = Array2::<f64>::zeros((n, n));
        
        for i in 1..n-1 {
            for j in 1..n-1 {
                let dvdx = (v[[i+1,j]] - v[[i-1,j]]) / (2.0 * self.dx);
                let dudy = (u[[i,j+1]] - u[[i,j-1]]) / (2.0 * self.dy);
                vorticity[[i,j]] = dvdx - dudy;
            }
        }
        
        vorticity
    }

    fn calculate_stream_function(&self, u: &Array2<f64>, _v: &Array2<f64>) -> Result<Array2<f64>, Box<dyn Error>> {
        let n = self.grid_resolution;
        let mut psi = Array2::<f64>::zeros((n, n));
        
        // Simple integration method for stream function
        // ∂ψ/∂y = u, ∂ψ/∂x = -v
        for i in 1..n {
            for j in 1..n {
                psi[[i,j]] = psi[[i,j-1]] + u[[i,j]] * self.dy;
            }
        }
        
        Ok(psi)
    }

    fn calculate_velocity_magnitude(&self, u: &Array2<f64>, v: &Array2<f64>) -> Array2<f64> {
        let n = self.grid_resolution;
        let mut magnitude = Array2::<f64>::zeros((n, n));
        
        for i in 0..n {
            for j in 0..n {
                magnitude[[i,j]] = (u[[i,j]].powi(2) + v[[i,j]].powi(2)).sqrt();
            }
        }
        
        magnitude
    }
}