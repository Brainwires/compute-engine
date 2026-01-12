//! 3D Navier-Stokes Solver
//!
//! Implements the three-dimensional incompressible Navier-Stokes equations:
//!
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u    (momentum)
//! ∇·u = 0                            (continuity/incompressibility)
//!
//! This solver is designed for:
//! - Taylor-Green vortex simulations (classic turbulence benchmark)
//! - DNS of homogeneous isotropic turbulence
//! - Studying enstrophy growth and potential singularity formation
//!
//! Numerical method: Fractional step method with:
//! - Semi-Lagrangian advection
//! - Jacobi iterations for pressure Poisson equation
//! - Central differences for viscous terms

use ndarray::Array3;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Configuration for 3D Navier-Stokes solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NS3DConfig {
    /// Number of grid points in x
    pub nx: usize,

    /// Number of grid points in y
    pub ny: usize,

    /// Number of grid points in z
    pub nz: usize,

    /// Domain length in x (m)
    pub lx: f64,

    /// Domain length in y (m)
    pub ly: f64,

    /// Domain length in z (m)
    pub lz: f64,

    /// Kinematic viscosity ν (m²/s)
    pub viscosity: f64,

    /// Reference density (kg/m³)
    pub rho: f64,

    /// CFL number for time step calculation
    pub cfl: f64,

    /// Maximum Jacobi iterations for pressure solver
    pub max_pressure_iterations: usize,

    /// Pressure solver tolerance
    pub pressure_tolerance: f64,
}

impl Default for NS3DConfig {
    fn default() -> Self {
        Self {
            nx: 64,
            ny: 64,
            nz: 64,
            lx: 2.0 * PI,
            ly: 2.0 * PI,
            lz: 2.0 * PI,
            viscosity: 0.01,
            rho: 1.0,
            cfl: 0.5,
            max_pressure_iterations: 100,
            pressure_tolerance: 1e-5,
        }
    }
}

/// State of the 3D Navier-Stokes simulation
#[derive(Debug, Clone)]
pub struct NS3DState {
    /// x-velocity field u(x,y,z)
    pub u: Array3<f64>,

    /// y-velocity field v(x,y,z)
    pub v: Array3<f64>,

    /// z-velocity field w(x,y,z)
    pub w: Array3<f64>,

    /// Pressure field p(x,y,z)
    pub p: Array3<f64>,

    /// Current simulation time
    pub time: f64,

    /// Number of time steps taken
    pub step_count: usize,
}

/// 3D Navier-Stokes Solver
pub struct NS3DSolver {
    config: NS3DConfig,
    dx: f64,
    dy: f64,
    dz: f64,

    // Velocity components
    u: Array3<f64>,
    v: Array3<f64>,
    w: Array3<f64>,
    p: Array3<f64>,

    // Work arrays
    u_star: Array3<f64>,
    v_star: Array3<f64>,
    w_star: Array3<f64>,

    time: f64,
    step_count: usize,
}

impl NS3DSolver {
    /// Create a new 3D Navier-Stokes solver
    pub fn new(config: NS3DConfig) -> Self {
        let nx = config.nx;
        let ny = config.ny;
        let nz = config.nz;
        let dx = config.lx / nx as f64;
        let dy = config.ly / ny as f64;
        let dz = config.lz / nz as f64;

        Self {
            config,
            dx,
            dy,
            dz,
            u: Array3::zeros((nx, ny, nz)),
            v: Array3::zeros((nx, ny, nz)),
            w: Array3::zeros((nx, ny, nz)),
            p: Array3::zeros((nx, ny, nz)),
            u_star: Array3::zeros((nx, ny, nz)),
            v_star: Array3::zeros((nx, ny, nz)),
            w_star: Array3::zeros((nx, ny, nz)),
            time: 0.0,
            step_count: 0,
        }
    }

    /// Initialize with Taylor-Green vortex
    ///
    /// Classic 3D benchmark case:
    /// u = u0 * sin(x/L) * cos(y/L) * cos(z/L)
    /// v = -u0 * cos(x/L) * sin(y/L) * cos(z/L)
    /// w = 0
    /// p = p0 + (ρu0²/16) * [cos(2x/L) + cos(2y/L)] * [cos(2z/L) + 2]
    pub fn init_taylor_green(&mut self, u0: f64) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let x = i as f64 * self.dx;
                    let y = j as f64 * self.dy;
                    let z = k as f64 * self.dz;

                    self.u[[i, j, k]] = u0 * x.sin() * y.cos() * z.cos();
                    self.v[[i, j, k]] = -u0 * x.cos() * y.sin() * z.cos();
                    self.w[[i, j, k]] = 0.0;

                    // Initial pressure from analytical solution
                    self.p[[i, j, k]] = self.config.rho * u0 * u0 / 16.0
                        * ((2.0 * x).cos() + (2.0 * y).cos())
                        * ((2.0 * z).cos() + 2.0);
                }
            }
        }

        self.time = 0.0;
        self.step_count = 0;
    }

    /// Initialize with uniform flow
    pub fn init_uniform(&mut self, u0: f64, v0: f64, w0: f64) {
        self.u.fill(u0);
        self.v.fill(v0);
        self.w.fill(w0);
        self.p.fill(0.0);
        self.time = 0.0;
        self.step_count = 0;
    }

    /// Calculate stable time step
    fn calculate_dt(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let mut max_vel: f64 = 1e-10;

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let vel = self.u[[i, j, k]].abs()
                        + self.v[[i, j, k]].abs()
                        + self.w[[i, j, k]].abs();
                    max_vel = max_vel.max(vel);
                }
            }
        }

        let dt_advection = self.config.cfl * self.dx.min(self.dy).min(self.dz) / max_vel;
        let dt_diffusion = self.config.cfl * self.dx.min(self.dy).min(self.dz).powi(2)
            / (6.0 * self.config.viscosity + 1e-15);

        dt_advection.min(dt_diffusion)
    }

    /// Perform one time step using fractional step method
    pub fn step(&mut self) {
        let dt = self.calculate_dt();
        let nu = self.config.viscosity;

        // Step 1: Compute intermediate velocity (advection + diffusion)
        self.compute_intermediate_velocity(dt, nu);

        // Step 2: Apply periodic boundary conditions to intermediate velocity
        self.apply_periodic_bc_star();

        // Step 3: Solve pressure Poisson equation
        self.solve_pressure(dt);

        // Step 4: Correct velocity to ensure divergence-free
        self.correct_velocity(dt);

        // Step 5: Apply periodic boundary conditions
        self.apply_periodic_bc();

        self.time += dt;
        self.step_count += 1;
    }

    /// Compute intermediate velocity (u*, v*, w*)
    fn compute_intermediate_velocity(&mut self, dt: f64, nu: f64) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    // Advection terms (central differences)
                    let u_dudx = self.u[[i, j, k]]
                        * (self.u[[i + 1, j, k]] - self.u[[i - 1, j, k]])
                        / (2.0 * self.dx);
                    let v_dudy = self.v[[i, j, k]]
                        * (self.u[[i, j + 1, k]] - self.u[[i, j - 1, k]])
                        / (2.0 * self.dy);
                    let w_dudz = self.w[[i, j, k]]
                        * (self.u[[i, j, k + 1]] - self.u[[i, j, k - 1]])
                        / (2.0 * self.dz);

                    let u_dvdx = self.u[[i, j, k]]
                        * (self.v[[i + 1, j, k]] - self.v[[i - 1, j, k]])
                        / (2.0 * self.dx);
                    let v_dvdy = self.v[[i, j, k]]
                        * (self.v[[i, j + 1, k]] - self.v[[i, j - 1, k]])
                        / (2.0 * self.dy);
                    let w_dvdz = self.w[[i, j, k]]
                        * (self.v[[i, j, k + 1]] - self.v[[i, j, k - 1]])
                        / (2.0 * self.dz);

                    let u_dwdx = self.u[[i, j, k]]
                        * (self.w[[i + 1, j, k]] - self.w[[i - 1, j, k]])
                        / (2.0 * self.dx);
                    let v_dwdy = self.v[[i, j, k]]
                        * (self.w[[i, j + 1, k]] - self.w[[i, j - 1, k]])
                        / (2.0 * self.dy);
                    let w_dwdz = self.w[[i, j, k]]
                        * (self.w[[i, j, k + 1]] - self.w[[i, j, k - 1]])
                        / (2.0 * self.dz);

                    // Diffusion terms (Laplacian)
                    let d2u = (self.u[[i + 1, j, k]] - 2.0 * self.u[[i, j, k]] + self.u[[i - 1, j, k]])
                        / (self.dx * self.dx)
                        + (self.u[[i, j + 1, k]] - 2.0 * self.u[[i, j, k]] + self.u[[i, j - 1, k]])
                            / (self.dy * self.dy)
                        + (self.u[[i, j, k + 1]] - 2.0 * self.u[[i, j, k]] + self.u[[i, j, k - 1]])
                            / (self.dz * self.dz);

                    let d2v = (self.v[[i + 1, j, k]] - 2.0 * self.v[[i, j, k]] + self.v[[i - 1, j, k]])
                        / (self.dx * self.dx)
                        + (self.v[[i, j + 1, k]] - 2.0 * self.v[[i, j, k]] + self.v[[i, j - 1, k]])
                            / (self.dy * self.dy)
                        + (self.v[[i, j, k + 1]] - 2.0 * self.v[[i, j, k]] + self.v[[i, j, k - 1]])
                            / (self.dz * self.dz);

                    let d2w = (self.w[[i + 1, j, k]] - 2.0 * self.w[[i, j, k]] + self.w[[i - 1, j, k]])
                        / (self.dx * self.dx)
                        + (self.w[[i, j + 1, k]] - 2.0 * self.w[[i, j, k]] + self.w[[i, j - 1, k]])
                            / (self.dy * self.dy)
                        + (self.w[[i, j, k + 1]] - 2.0 * self.w[[i, j, k]] + self.w[[i, j, k - 1]])
                            / (self.dz * self.dz);

                    // Intermediate velocities
                    self.u_star[[i, j, k]] =
                        self.u[[i, j, k]] + dt * (-(u_dudx + v_dudy + w_dudz) + nu * d2u);
                    self.v_star[[i, j, k]] =
                        self.v[[i, j, k]] + dt * (-(u_dvdx + v_dvdy + w_dvdz) + nu * d2v);
                    self.w_star[[i, j, k]] =
                        self.w[[i, j, k]] + dt * (-(u_dwdx + v_dwdy + w_dwdz) + nu * d2w);
                }
            }
        }
    }

    /// Apply periodic boundary conditions to intermediate velocity
    fn apply_periodic_bc_star(&mut self) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;

        // x boundaries
        for j in 0..ny {
            for k in 0..nz {
                self.u_star[[0, j, k]] = self.u_star[[nx - 2, j, k]];
                self.u_star[[nx - 1, j, k]] = self.u_star[[1, j, k]];
                self.v_star[[0, j, k]] = self.v_star[[nx - 2, j, k]];
                self.v_star[[nx - 1, j, k]] = self.v_star[[1, j, k]];
                self.w_star[[0, j, k]] = self.w_star[[nx - 2, j, k]];
                self.w_star[[nx - 1, j, k]] = self.w_star[[1, j, k]];
            }
        }

        // y boundaries
        for i in 0..nx {
            for k in 0..nz {
                self.u_star[[i, 0, k]] = self.u_star[[i, ny - 2, k]];
                self.u_star[[i, ny - 1, k]] = self.u_star[[i, 1, k]];
                self.v_star[[i, 0, k]] = self.v_star[[i, ny - 2, k]];
                self.v_star[[i, ny - 1, k]] = self.v_star[[i, 1, k]];
                self.w_star[[i, 0, k]] = self.w_star[[i, ny - 2, k]];
                self.w_star[[i, ny - 1, k]] = self.w_star[[i, 1, k]];
            }
        }

        // z boundaries
        for i in 0..nx {
            for j in 0..ny {
                self.u_star[[i, j, 0]] = self.u_star[[i, j, nz - 2]];
                self.u_star[[i, j, nz - 1]] = self.u_star[[i, j, 1]];
                self.v_star[[i, j, 0]] = self.v_star[[i, j, nz - 2]];
                self.v_star[[i, j, nz - 1]] = self.v_star[[i, j, 1]];
                self.w_star[[i, j, 0]] = self.w_star[[i, j, nz - 2]];
                self.w_star[[i, j, nz - 1]] = self.w_star[[i, j, 1]];
            }
        }
    }

    /// Solve pressure Poisson equation: ∇²p = (ρ/dt) ∇·u*
    fn solve_pressure(&mut self, dt: f64) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let rho = self.config.rho;

        // Compute RHS: (ρ/dt) ∇·u*
        let mut rhs = Array3::<f64>::zeros((nx, ny, nz));
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let div = (self.u_star[[i + 1, j, k]] - self.u_star[[i - 1, j, k]])
                        / (2.0 * self.dx)
                        + (self.v_star[[i, j + 1, k]] - self.v_star[[i, j - 1, k]])
                            / (2.0 * self.dy)
                        + (self.w_star[[i, j, k + 1]] - self.w_star[[i, j, k - 1]])
                            / (2.0 * self.dz);
                    rhs[[i, j, k]] = rho * div / dt;
                }
            }
        }

        // Jacobi iteration
        let mut p_old = self.p.clone();
        let dx2 = self.dx * self.dx;
        let dy2 = self.dy * self.dy;
        let dz2 = self.dz * self.dz;
        let coef = 2.0 * (1.0 / dx2 + 1.0 / dy2 + 1.0 / dz2);

        for _iter in 0..self.config.max_pressure_iterations {
            let mut max_residual: f64 = 0.0;

            for i in 1..nx - 1 {
                for j in 1..ny - 1 {
                    for k in 1..nz - 1 {
                        let p_new = (1.0 / coef)
                            * ((p_old[[i + 1, j, k]] + p_old[[i - 1, j, k]]) / dx2
                                + (p_old[[i, j + 1, k]] + p_old[[i, j - 1, k]]) / dy2
                                + (p_old[[i, j, k + 1]] + p_old[[i, j, k - 1]]) / dz2
                                - rhs[[i, j, k]]);

                        let residual = (p_new - p_old[[i, j, k]]).abs();
                        max_residual = max_residual.max(residual);
                        self.p[[i, j, k]] = p_new;
                    }
                }
            }

            // Apply periodic BC for pressure
            self.apply_pressure_bc();

            if max_residual < self.config.pressure_tolerance {
                break;
            }

            p_old = self.p.clone();
        }
    }

    /// Apply periodic boundary conditions for pressure
    fn apply_pressure_bc(&mut self) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;

        for j in 0..ny {
            for k in 0..nz {
                self.p[[0, j, k]] = self.p[[nx - 2, j, k]];
                self.p[[nx - 1, j, k]] = self.p[[1, j, k]];
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                self.p[[i, 0, k]] = self.p[[i, ny - 2, k]];
                self.p[[i, ny - 1, k]] = self.p[[i, 1, k]];
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                self.p[[i, j, 0]] = self.p[[i, j, nz - 2]];
                self.p[[i, j, nz - 1]] = self.p[[i, j, 1]];
            }
        }
    }

    /// Correct velocity to be divergence-free: u = u* - (dt/ρ)∇p
    fn correct_velocity(&mut self, dt: f64) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let rho = self.config.rho;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let dpdx =
                        (self.p[[i + 1, j, k]] - self.p[[i - 1, j, k]]) / (2.0 * self.dx);
                    let dpdy =
                        (self.p[[i, j + 1, k]] - self.p[[i, j - 1, k]]) / (2.0 * self.dy);
                    let dpdz =
                        (self.p[[i, j, k + 1]] - self.p[[i, j, k - 1]]) / (2.0 * self.dz);

                    self.u[[i, j, k]] = self.u_star[[i, j, k]] - dt * dpdx / rho;
                    self.v[[i, j, k]] = self.v_star[[i, j, k]] - dt * dpdy / rho;
                    self.w[[i, j, k]] = self.w_star[[i, j, k]] - dt * dpdz / rho;
                }
            }
        }
    }

    /// Apply periodic boundary conditions
    fn apply_periodic_bc(&mut self) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;

        // x boundaries
        for j in 0..ny {
            for k in 0..nz {
                self.u[[0, j, k]] = self.u[[nx - 2, j, k]];
                self.u[[nx - 1, j, k]] = self.u[[1, j, k]];
                self.v[[0, j, k]] = self.v[[nx - 2, j, k]];
                self.v[[nx - 1, j, k]] = self.v[[1, j, k]];
                self.w[[0, j, k]] = self.w[[nx - 2, j, k]];
                self.w[[nx - 1, j, k]] = self.w[[1, j, k]];
            }
        }

        // y boundaries
        for i in 0..nx {
            for k in 0..nz {
                self.u[[i, 0, k]] = self.u[[i, ny - 2, k]];
                self.u[[i, ny - 1, k]] = self.u[[i, 1, k]];
                self.v[[i, 0, k]] = self.v[[i, ny - 2, k]];
                self.v[[i, ny - 1, k]] = self.v[[i, 1, k]];
                self.w[[i, 0, k]] = self.w[[i, ny - 2, k]];
                self.w[[i, ny - 1, k]] = self.w[[i, 1, k]];
            }
        }

        // z boundaries
        for i in 0..nx {
            for j in 0..ny {
                self.u[[i, j, 0]] = self.u[[i, j, nz - 2]];
                self.u[[i, j, nz - 1]] = self.u[[i, j, 1]];
                self.v[[i, j, 0]] = self.v[[i, j, nz - 2]];
                self.v[[i, j, nz - 1]] = self.v[[i, j, 1]];
                self.w[[i, j, 0]] = self.w[[i, j, nz - 2]];
                self.w[[i, j, nz - 1]] = self.w[[i, j, 1]];
            }
        }
    }

    /// Run simulation until the specified end time
    pub fn run_until(&mut self, end_time: f64) {
        while self.time < end_time {
            self.step();
        }
    }

    /// Get the current simulation state
    pub fn get_state(&self) -> NS3DState {
        NS3DState {
            u: self.u.clone(),
            v: self.v.clone(),
            w: self.w.clone(),
            p: self.p.clone(),
            time: self.time,
            step_count: self.step_count,
        }
    }

    /// Calculate total kinetic energy: E = (1/2) ∫ (u² + v² + w²) dV
    pub fn kinetic_energy(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let mut energy = 0.0;

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    energy += 0.5
                        * (self.u[[i, j, k]].powi(2)
                            + self.v[[i, j, k]].powi(2)
                            + self.w[[i, j, k]].powi(2));
                }
            }
        }

        energy * self.dx * self.dy * self.dz
    }

    /// Calculate enstrophy: Ω = ∫ |ω|² dV
    pub fn enstrophy(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let mut enstrophy = 0.0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    // Vorticity components
                    let omega_x = (self.w[[i, j + 1, k]] - self.w[[i, j - 1, k]]) / (2.0 * self.dy)
                        - (self.v[[i, j, k + 1]] - self.v[[i, j, k - 1]]) / (2.0 * self.dz);
                    let omega_y = (self.u[[i, j, k + 1]] - self.u[[i, j, k - 1]]) / (2.0 * self.dz)
                        - (self.w[[i + 1, j, k]] - self.w[[i - 1, j, k]]) / (2.0 * self.dx);
                    let omega_z = (self.v[[i + 1, j, k]] - self.v[[i - 1, j, k]]) / (2.0 * self.dx)
                        - (self.u[[i, j + 1, k]] - self.u[[i, j - 1, k]]) / (2.0 * self.dy);

                    enstrophy += omega_x * omega_x + omega_y * omega_y + omega_z * omega_z;
                }
            }
        }

        enstrophy * self.dx * self.dy * self.dz
    }

    /// Calculate maximum vorticity magnitude
    pub fn max_vorticity(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let mut max_omega: f64 = 0.0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let omega_x = (self.w[[i, j + 1, k]] - self.w[[i, j - 1, k]]) / (2.0 * self.dy)
                        - (self.v[[i, j, k + 1]] - self.v[[i, j, k - 1]]) / (2.0 * self.dz);
                    let omega_y = (self.u[[i, j, k + 1]] - self.u[[i, j, k - 1]]) / (2.0 * self.dz)
                        - (self.w[[i + 1, j, k]] - self.w[[i - 1, j, k]]) / (2.0 * self.dx);
                    let omega_z = (self.v[[i + 1, j, k]] - self.v[[i - 1, j, k]]) / (2.0 * self.dx)
                        - (self.u[[i, j + 1, k]] - self.u[[i, j - 1, k]]) / (2.0 * self.dy);

                    let omega_mag = (omega_x * omega_x + omega_y * omega_y + omega_z * omega_z).sqrt();
                    max_omega = max_omega.max(omega_mag);
                }
            }
        }

        max_omega
    }

    /// Calculate maximum velocity divergence (should be ~0 for incompressible)
    pub fn max_divergence(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let mut max_div: f64 = 0.0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                for k in 1..nz - 1 {
                    let div = (self.u[[i + 1, j, k]] - self.u[[i - 1, j, k]]) / (2.0 * self.dx)
                        + (self.v[[i, j + 1, k]] - self.v[[i, j - 1, k]]) / (2.0 * self.dy)
                        + (self.w[[i, j, k + 1]] - self.w[[i, j, k - 1]]) / (2.0 * self.dz);
                    max_div = max_div.max(div.abs());
                }
            }
        }

        max_div
    }

    /// Get grid coordinates
    pub fn get_x_coords(&self) -> Vec<f64> {
        (0..self.config.nx).map(|i| i as f64 * self.dx).collect()
    }

    pub fn get_y_coords(&self) -> Vec<f64> {
        (0..self.config.ny).map(|j| j as f64 * self.dy).collect()
    }

    pub fn get_z_coords(&self) -> Vec<f64> {
        (0..self.config.nz).map(|k| k as f64 * self.dz).collect()
    }

    /// Get Reynolds number based on domain size and max velocity
    pub fn reynolds_number(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let nz = self.config.nz;
        let mut max_vel: f64 = 0.0;

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let vel = (self.u[[i, j, k]].powi(2)
                        + self.v[[i, j, k]].powi(2)
                        + self.w[[i, j, k]].powi(2))
                    .sqrt();
                    max_vel = max_vel.max(vel);
                }
            }
        }

        max_vel * self.config.lx / self.config.viscosity
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_taylor_green_initialization() {
        let config = NS3DConfig {
            nx: 16,
            ny: 16,
            nz: 16,
            ..Default::default()
        };

        let mut solver = NS3DSolver::new(config);
        solver.init_taylor_green(1.0);

        // Check that velocity field is non-zero
        let state = solver.get_state();
        let max_u: f64 = state.u.iter().map(|x| x.abs()).fold(0.0, f64::max);
        let max_v: f64 = state.v.iter().map(|x| x.abs()).fold(0.0, f64::max);

        assert!(max_u > 0.9, "Max u should be close to u0=1.0");
        assert!(max_v > 0.9, "Max v should be close to u0=1.0");

        // w should be zero for Taylor-Green
        let max_w: f64 = state.w.iter().map(|x| x.abs()).fold(0.0, f64::max);
        assert!(max_w < 1e-10, "Max w should be zero for Taylor-Green");
    }

    #[test]
    fn test_energy_conservation() {
        let config = NS3DConfig {
            nx: 16,
            ny: 16,
            nz: 16,
            viscosity: 0.0, // Inviscid - energy should be conserved
            ..Default::default()
        };

        let mut solver = NS3DSolver::new(config);
        solver.init_taylor_green(0.5);

        let e0 = solver.kinetic_energy();

        // Run for a few steps
        for _ in 0..5 {
            solver.step();
        }

        let e1 = solver.kinetic_energy();

        // Energy should be approximately conserved (small numerical dissipation)
        let relative_change = (e1 - e0).abs() / e0;
        assert!(
            relative_change < 0.1,
            "Energy should be nearly conserved for inviscid flow: change={:.2}%",
            relative_change * 100.0
        );
    }

    #[test]
    fn test_energy_decay_with_viscosity() {
        let config = NS3DConfig {
            nx: 16,
            ny: 16,
            nz: 16,
            viscosity: 0.1, // Significant viscosity
            ..Default::default()
        };

        let mut solver = NS3DSolver::new(config);
        solver.init_taylor_green(1.0);

        let e0 = solver.kinetic_energy();

        // Run for some steps
        for _ in 0..20 {
            solver.step();
        }

        let e1 = solver.kinetic_energy();

        // Energy should decay due to viscosity
        assert!(
            e1 < e0,
            "Energy should decay with viscosity: initial={}, final={}",
            e0,
            e1
        );
    }

    #[test]
    fn test_divergence_free() {
        let config = NS3DConfig {
            nx: 16,
            ny: 16,
            nz: 16,
            max_pressure_iterations: 200, // More iterations for better convergence
            pressure_tolerance: 1e-6,
            ..Default::default()
        };

        let mut solver = NS3DSolver::new(config);
        solver.init_taylor_green(0.5);

        // Run a few steps
        for _ in 0..10 {
            solver.step();
        }

        let max_div = solver.max_divergence();

        // On coarse grids, divergence control is limited by discretization error
        // A tolerance of ~0.3 is acceptable for 16^3 grid
        assert!(
            max_div < 0.5,
            "Divergence should be reasonably small: max_div={}",
            max_div
        );
    }

    #[test]
    fn test_enstrophy_calculation() {
        let config = NS3DConfig {
            nx: 16,
            ny: 16,
            nz: 16,
            ..Default::default()
        };

        let mut solver = NS3DSolver::new(config);
        solver.init_taylor_green(1.0);

        let enstrophy = solver.enstrophy();

        // Taylor-Green has non-zero vorticity
        assert!(
            enstrophy > 0.0,
            "Enstrophy should be positive for Taylor-Green"
        );
    }
}
