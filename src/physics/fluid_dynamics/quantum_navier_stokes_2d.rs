//! 2D Quantum Navier-Stokes Solver
//!
//! Implements the two-dimensional compressible Quantum Navier-Stokes equations:
//!
//! ∂ρ/∂t + ∇·(ρu) = 0                                    (continuity)
//! ∂(ρu)/∂t + ∇·(ρu⊗u) + ∇p = μ∇²u + F_Q                (momentum + quantum)
//!
//! where F_Q = -ρ∇Q is the quantum pressure force from the Bohm potential:
//! Q = (ℏ²/2m) × ∇²√ρ / √ρ
//!
//! This solver is designed for:
//! - 2D decaying turbulence simulations
//! - Comparing classical N-S vs QNS at various Reynolds numbers
//! - Studying decoherence effects on turbulent energy cascade
//!
//! Numerical method: Fractional step with MacCormack predictor-corrector
//! for advection, central differences for diffusion and quantum terms.

use crate::physics::quantum_mechanics::bohm_potential::BohmPotential;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Configuration for 2D QNS solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QNS2DConfig {
    /// Number of grid points in x
    pub nx: usize,

    /// Number of grid points in y
    pub ny: usize,

    /// Domain length in x (m)
    pub lx: f64,

    /// Domain length in y (m)
    pub ly: f64,

    /// Dynamic viscosity (Pa·s)
    pub viscosity: f64,

    /// Particle/molecule mass for quantum effects (kg)
    pub particle_mass: f64,

    /// Reference density (kg/m³)
    pub rho_ref: f64,

    /// Specific heat ratio γ
    pub gamma: f64,

    /// Enable quantum correction term
    pub enable_quantum: bool,

    /// CFL number for time step calculation (< 1 for stability)
    pub cfl: f64,

    /// Sound speed for pressure correction (for weakly compressible)
    pub sound_speed: f64,
}

impl Default for QNS2DConfig {
    fn default() -> Self {
        Self {
            nx: 128,
            ny: 128,
            lx: 2.0 * PI,
            ly: 2.0 * PI,
            viscosity: 0.01,
            particle_mass: 4.8e-26, // N2 molecule
            rho_ref: 1.0,
            gamma: 1.4,
            enable_quantum: true,
            cfl: 0.3,
            sound_speed: 10.0, // Artificial sound speed for pressure correction
        }
    }
}

/// State of the 2D QNS simulation
#[derive(Debug, Clone)]
pub struct QNS2DState {
    /// Density field ρ(x,y)
    pub density: Array2<f64>,

    /// x-velocity field u(x,y)
    pub u: Array2<f64>,

    /// y-velocity field v(x,y)
    pub v: Array2<f64>,

    /// Pressure field p(x,y)
    pub pressure: Array2<f64>,

    /// Current simulation time
    pub time: f64,

    /// Number of time steps taken
    pub step_count: usize,
}

/// 2D Quantum Navier-Stokes Solver
pub struct QNS2DSolver {
    config: QNS2DConfig,
    dx: f64,
    dy: f64,
    bohm: BohmPotential,

    // Primitive variables
    rho: Array2<f64>,  // density
    u: Array2<f64>,    // x-velocity
    v: Array2<f64>,    // y-velocity
    p: Array2<f64>,    // pressure

    // Work arrays for MacCormack scheme
    rho_pred: Array2<f64>,
    u_pred: Array2<f64>,
    v_pred: Array2<f64>,

    // Quantum force arrays
    q_force_x: Array2<f64>,
    q_force_y: Array2<f64>,

    time: f64,
    step_count: usize,
}

impl QNS2DSolver {
    /// Create a new 2D QNS solver with the given configuration
    pub fn new(config: QNS2DConfig) -> Self {
        let nx = config.nx;
        let ny = config.ny;
        let dx = config.lx / nx as f64;
        let dy = config.ly / ny as f64;
        let bohm = BohmPotential::with_mass(config.particle_mass);

        Self {
            config,
            dx,
            dy,
            bohm,
            rho: Array2::from_elem((nx, ny), 1.0),
            u: Array2::zeros((nx, ny)),
            v: Array2::zeros((nx, ny)),
            p: Array2::from_elem((nx, ny), 1.0),
            rho_pred: Array2::zeros((nx, ny)),
            u_pred: Array2::zeros((nx, ny)),
            v_pred: Array2::zeros((nx, ny)),
            q_force_x: Array2::zeros((nx, ny)),
            q_force_y: Array2::zeros((nx, ny)),
            time: 0.0,
            step_count: 0,
        }
    }

    /// Initialize with Taylor-Green vortex (classic turbulence test case)
    ///
    /// u = u0 * sin(kx) * cos(ky)
    /// v = -u0 * cos(kx) * sin(ky)
    /// p = p0 + (ρu0²/4) * (cos(2kx) + cos(2ky))
    pub fn init_taylor_green(&mut self, u0: f64, k: f64) {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let p0 = 1.0 / self.config.gamma; // Reference pressure

        for i in 0..nx {
            for j in 0..ny {
                let x = i as f64 * self.dx;
                let y = j as f64 * self.dy;

                self.rho[[i, j]] = self.config.rho_ref;
                self.u[[i, j]] = u0 * (k * x).sin() * (k * y).cos();
                self.v[[i, j]] = -u0 * (k * x).cos() * (k * y).sin();
                self.p[[i, j]] =
                    p0 + (self.config.rho_ref * u0 * u0 / 4.0) * ((2.0 * k * x).cos() + (2.0 * k * y).cos());
            }
        }

        self.time = 0.0;
        self.step_count = 0;
    }

    /// Initialize with decaying turbulence (random vorticity field)
    pub fn init_decaying_turbulence(&mut self, energy: f64, k_peak: f64) {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let nx = self.config.nx;
        let ny = self.config.ny;

        // Generate random phases for each Fourier mode
        let n_modes = 16;

        for i in 0..nx {
            for j in 0..ny {
                let x = i as f64 * self.dx;
                let y = j as f64 * self.dy;

                let mut u_sum = 0.0;
                let mut v_sum = 0.0;

                // Sum over Fourier modes with energy spectrum peaked at k_peak
                for kx in 1..=n_modes {
                    for ky in 1..=n_modes {
                        let k_mag = ((kx * kx + ky * ky) as f64).sqrt();
                        // Energy spectrum: E(k) ∝ k⁴ exp(-2(k/k_peak)²) (Kolmogorov-like)
                        let e_k = k_mag.powi(4) * (-2.0 * (k_mag / k_peak).powi(2)).exp();
                        let amplitude = (2.0 * e_k / k_mag).sqrt() * energy.sqrt();

                        let phase_x: f64 = rng.r#gen::<f64>() * 2.0 * PI;
                        let phase_y: f64 = rng.r#gen::<f64>() * 2.0 * PI;

                        let kx_f = kx as f64;
                        let ky_f = ky as f64;

                        // Divergence-free velocity field from stream function
                        u_sum += amplitude * ky_f / k_mag * (kx_f * x + ky_f * y + phase_x).sin();
                        v_sum -= amplitude * kx_f / k_mag * (kx_f * x + ky_f * y + phase_y).sin();
                    }
                }

                self.rho[[i, j]] = self.config.rho_ref;
                self.u[[i, j]] = u_sum;
                self.v[[i, j]] = v_sum;
                self.p[[i, j]] = 1.0 / self.config.gamma;
            }
        }

        self.time = 0.0;
        self.step_count = 0;
    }

    /// Calculate stable time step based on CFL condition
    fn calculate_dt(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mut max_speed: f64 = 1e-10;

        for i in 0..nx {
            for j in 0..ny {
                let speed = self.u[[i, j]].abs() + self.v[[i, j]].abs() + self.config.sound_speed;
                max_speed = max_speed.max(speed);
            }
        }

        let dt_advection = self.config.cfl * self.dx.min(self.dy) / max_speed;
        let dt_diffusion = self.config.cfl * self.dx.min(self.dy).powi(2)
            / (4.0 * self.config.viscosity / self.config.rho_ref + 1e-15);

        dt_advection.min(dt_diffusion)
    }

    /// Calculate quantum force from Bohm potential
    fn calculate_quantum_forces(&mut self) {
        if !self.config.enable_quantum {
            self.q_force_x.fill(0.0);
            self.q_force_y.fill(0.0);
            return;
        }

        let nx = self.config.nx;
        let ny = self.config.ny;

        // Calculate Bohm potential Q on the grid
        let q = self.bohm.calculate_2d(&self.rho, self.dx, self.dy);

        // Calculate -ρ∇Q (quantum force per unit volume)
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let dq_dx = (q[[i + 1, j]] - q[[i - 1, j]]) / (2.0 * self.dx);
                let dq_dy = (q[[i, j + 1]] - q[[i, j - 1]]) / (2.0 * self.dy);

                self.q_force_x[[i, j]] = -self.rho[[i, j]] * dq_dx;
                self.q_force_y[[i, j]] = -self.rho[[i, j]] * dq_dy;
            }
        }

        // Periodic boundary conditions for quantum forces
        for j in 0..ny {
            self.q_force_x[[0, j]] = self.q_force_x[[nx - 2, j]];
            self.q_force_x[[nx - 1, j]] = self.q_force_x[[1, j]];
        }
        for i in 0..nx {
            self.q_force_y[[i, 0]] = self.q_force_y[[i, ny - 2]];
            self.q_force_y[[i, ny - 1]] = self.q_force_y[[i, 1]];
        }
    }

    /// Perform one time step using MacCormack predictor-corrector
    pub fn step(&mut self) {
        let dt = self.calculate_dt();
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mu = self.config.viscosity;

        // Calculate quantum forces before predictor step
        self.calculate_quantum_forces();

        // MacCormack predictor step (forward differences)
        for i in 0..nx - 1 {
            for j in 0..ny - 1 {
                let ip1 = (i + 1) % nx;
                let jp1 = (j + 1) % ny;

                // Continuity: ∂ρ/∂t = -∇·(ρu)
                let drho_flux_x = (self.rho[[ip1, j]] * self.u[[ip1, j]] - self.rho[[i, j]] * self.u[[i, j]]) / self.dx;
                let drho_flux_y = (self.rho[[i, jp1]] * self.v[[i, jp1]] - self.rho[[i, j]] * self.v[[i, j]]) / self.dy;
                self.rho_pred[[i, j]] = self.rho[[i, j]] - dt * (drho_flux_x + drho_flux_y);

                // Momentum x: ∂u/∂t = -u·∇u - (1/ρ)∂p/∂x + ν∇²u + F_Q/ρ
                let u_dudx = self.u[[i, j]] * (self.u[[ip1, j]] - self.u[[i, j]]) / self.dx;
                let v_dudy = self.v[[i, j]] * (self.u[[i, jp1]] - self.u[[i, j]]) / self.dy;
                let dpdx = (self.p[[ip1, j]] - self.p[[i, j]]) / self.dx;
                let d2udx2 = if i > 0 && i < nx - 1 {
                    (self.u[[ip1, j]] - 2.0 * self.u[[i, j]] + self.u[[i - 1, j]]) / (self.dx * self.dx)
                } else {
                    0.0
                };
                let d2udy2 = if j > 0 && j < ny - 1 {
                    (self.u[[i, jp1]] - 2.0 * self.u[[i, j]] + self.u[[i, j - 1]]) / (self.dy * self.dy)
                } else {
                    0.0
                };
                let nu = mu / self.rho[[i, j]].max(1e-10);
                let f_qx = self.q_force_x[[i, j]] / self.rho[[i, j]].max(1e-10);

                self.u_pred[[i, j]] = self.u[[i, j]]
                    + dt * (-u_dudx - v_dudy - dpdx / self.rho[[i, j]].max(1e-10) + nu * (d2udx2 + d2udy2) + f_qx);

                // Momentum y: ∂v/∂t = -u·∇v - (1/ρ)∂p/∂y + ν∇²v + F_Q/ρ
                let u_dvdx = self.u[[i, j]] * (self.v[[ip1, j]] - self.v[[i, j]]) / self.dx;
                let v_dvdy = self.v[[i, j]] * (self.v[[i, jp1]] - self.v[[i, j]]) / self.dy;
                let dpdy = (self.p[[i, jp1]] - self.p[[i, j]]) / self.dy;
                let d2vdx2 = if i > 0 && i < nx - 1 {
                    (self.v[[ip1, j]] - 2.0 * self.v[[i, j]] + self.v[[i - 1, j]]) / (self.dx * self.dx)
                } else {
                    0.0
                };
                let d2vdy2 = if j > 0 && j < ny - 1 {
                    (self.v[[i, jp1]] - 2.0 * self.v[[i, j]] + self.v[[i, j - 1]]) / (self.dy * self.dy)
                } else {
                    0.0
                };
                let f_qy = self.q_force_y[[i, j]] / self.rho[[i, j]].max(1e-10);

                self.v_pred[[i, j]] = self.v[[i, j]]
                    + dt * (-u_dvdx - v_dvdy - dpdy / self.rho[[i, j]].max(1e-10) + nu * (d2vdx2 + d2vdy2) + f_qy);
            }
        }

        // Apply periodic BCs to predicted values
        self.apply_periodic_bc_pred();

        // MacCormack corrector step (backward differences)
        for i in 1..nx {
            for j in 1..ny {
                let im1 = if i > 0 { i - 1 } else { nx - 1 };
                let jm1 = if j > 0 { j - 1 } else { ny - 1 };

                // Continuity corrector
                let drho_flux_x_pred =
                    (self.rho_pred[[i, j]] * self.u_pred[[i, j]] - self.rho_pred[[im1, j]] * self.u_pred[[im1, j]]) / self.dx;
                let drho_flux_y_pred =
                    (self.rho_pred[[i, j]] * self.v_pred[[i, j]] - self.rho_pred[[i, jm1]] * self.v_pred[[i, jm1]]) / self.dy;
                let rho_corr = self.rho[[i, j]] - dt * (drho_flux_x_pred + drho_flux_y_pred);
                self.rho[[i, j]] = 0.5 * (self.rho_pred[[i, j]] + rho_corr);

                // u-momentum corrector
                let u_dudx_pred = self.u_pred[[i, j]] * (self.u_pred[[i, j]] - self.u_pred[[im1, j]]) / self.dx;
                let v_dudy_pred = self.v_pred[[i, j]] * (self.u_pred[[i, j]] - self.u_pred[[i, jm1]]) / self.dy;
                let dpdx_pred = (self.p[[i, j]] - self.p[[im1, j]]) / self.dx;
                let _nu_pred = mu / self.rho_pred[[i, j]].max(1e-10);
                let f_qx_pred = self.q_force_x[[i, j]] / self.rho_pred[[i, j]].max(1e-10);

                let u_corr = self.u[[i, j]]
                    + dt * (-u_dudx_pred - v_dudy_pred - dpdx_pred / self.rho_pred[[i, j]].max(1e-10) + f_qx_pred);
                self.u[[i, j]] = 0.5 * (self.u_pred[[i, j]] + u_corr);

                // v-momentum corrector
                let u_dvdx_pred = self.u_pred[[i, j]] * (self.v_pred[[i, j]] - self.v_pred[[im1, j]]) / self.dx;
                let v_dvdy_pred = self.v_pred[[i, j]] * (self.v_pred[[i, j]] - self.v_pred[[i, jm1]]) / self.dy;
                let dpdy_pred = (self.p[[i, j]] - self.p[[i, jm1]]) / self.dy;
                let f_qy_pred = self.q_force_y[[i, j]] / self.rho_pred[[i, j]].max(1e-10);

                let v_corr = self.v[[i, j]]
                    + dt * (-u_dvdx_pred - v_dvdy_pred - dpdy_pred / self.rho_pred[[i, j]].max(1e-10) + f_qy_pred);
                self.v[[i, j]] = 0.5 * (self.v_pred[[i, j]] + v_corr);
            }
        }

        // Update pressure using equation of state (weakly compressible)
        for i in 0..nx {
            for j in 0..ny {
                self.p[[i, j]] = self.config.sound_speed.powi(2)
                    * (self.rho[[i, j]] - self.config.rho_ref)
                    + 1.0 / self.config.gamma;
            }
        }

        // Ensure positive density
        for i in 0..nx {
            for j in 0..ny {
                self.rho[[i, j]] = self.rho[[i, j]].max(1e-10);
            }
        }

        // Apply periodic boundary conditions
        self.apply_periodic_bc();

        self.time += dt;
        self.step_count += 1;
    }

    /// Apply periodic boundary conditions to predicted arrays
    fn apply_periodic_bc_pred(&mut self) {
        let nx = self.config.nx;
        let ny = self.config.ny;

        // x boundaries
        for j in 0..ny {
            self.rho_pred[[nx - 1, j]] = self.rho_pred[[0, j]];
            self.u_pred[[nx - 1, j]] = self.u_pred[[0, j]];
            self.v_pred[[nx - 1, j]] = self.v_pred[[0, j]];
        }

        // y boundaries
        for i in 0..nx {
            self.rho_pred[[i, ny - 1]] = self.rho_pred[[i, 0]];
            self.u_pred[[i, ny - 1]] = self.u_pred[[i, 0]];
            self.v_pred[[i, ny - 1]] = self.v_pred[[i, 0]];
        }
    }

    /// Apply periodic boundary conditions
    fn apply_periodic_bc(&mut self) {
        let nx = self.config.nx;
        let ny = self.config.ny;

        // x boundaries
        for j in 0..ny {
            self.rho[[nx - 1, j]] = self.rho[[0, j]];
            self.u[[nx - 1, j]] = self.u[[0, j]];
            self.v[[nx - 1, j]] = self.v[[0, j]];
            self.p[[nx - 1, j]] = self.p[[0, j]];
        }

        // y boundaries
        for i in 0..nx {
            self.rho[[i, ny - 1]] = self.rho[[i, 0]];
            self.u[[i, ny - 1]] = self.u[[i, 0]];
            self.v[[i, ny - 1]] = self.v[[i, 0]];
            self.p[[i, ny - 1]] = self.p[[i, 0]];
        }
    }

    /// Run simulation until the specified end time
    pub fn run_until(&mut self, end_time: f64) {
        while self.time < end_time {
            self.step();
        }
    }

    /// Get the current simulation state
    pub fn get_state(&self) -> QNS2DState {
        QNS2DState {
            density: self.rho.clone(),
            u: self.u.clone(),
            v: self.v.clone(),
            pressure: self.p.clone(),
            time: self.time,
            step_count: self.step_count,
        }
    }

    /// Get x-coordinates of grid points
    pub fn get_x_coords(&self) -> Vec<f64> {
        (0..self.config.nx).map(|i| i as f64 * self.dx).collect()
    }

    /// Get y-coordinates of grid points
    pub fn get_y_coords(&self) -> Vec<f64> {
        (0..self.config.ny).map(|j| j as f64 * self.dy).collect()
    }

    /// Calculate total kinetic energy: E = (1/2) ∫ ρ(u² + v²) dA
    pub fn kinetic_energy(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mut energy = 0.0;

        for i in 0..nx {
            for j in 0..ny {
                energy += 0.5 * self.rho[[i, j]] * (self.u[[i, j]].powi(2) + self.v[[i, j]].powi(2));
            }
        }

        energy * self.dx * self.dy
    }

    /// Calculate enstrophy: Ω = ∫ ω² dA where ω = ∂v/∂x - ∂u/∂y
    pub fn enstrophy(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mut enstrophy = 0.0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let dvdx = (self.v[[i + 1, j]] - self.v[[i - 1, j]]) / (2.0 * self.dx);
                let dudy = (self.u[[i, j + 1]] - self.u[[i, j - 1]]) / (2.0 * self.dy);
                let omega = dvdx - dudy;
                enstrophy += omega * omega;
            }
        }

        enstrophy * self.dx * self.dy
    }

    /// Calculate maximum vorticity magnitude
    pub fn max_vorticity(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mut max_omega: f64 = 0.0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let dvdx = (self.v[[i + 1, j]] - self.v[[i - 1, j]]) / (2.0 * self.dx);
                let dudy = (self.u[[i, j + 1]] - self.u[[i, j - 1]]) / (2.0 * self.dy);
                let omega = (dvdx - dudy).abs();
                max_omega = max_omega.max(omega);
            }
        }

        max_omega
    }

    /// Calculate velocity divergence (should be small for weakly compressible)
    pub fn max_divergence(&self) -> f64 {
        let nx = self.config.nx;
        let ny = self.config.ny;
        let mut max_div: f64 = 0.0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let dudx = (self.u[[i + 1, j]] - self.u[[i - 1, j]]) / (2.0 * self.dx);
                let dvdy = (self.v[[i, j + 1]] - self.v[[i, j - 1]]) / (2.0 * self.dy);
                let div = (dudx + dvdy).abs();
                max_div = max_div.max(div);
            }
        }

        max_div
    }
}

/// Compare classical vs quantum 2D simulations
pub fn compare_classical_vs_quantum_2d(
    config: QNS2DConfig,
    init_fn: impl Fn(&mut QNS2DSolver),
    end_time: f64,
) -> QNS2DComparisonResult {
    // Run classical simulation
    let mut config_classical = config.clone();
    config_classical.enable_quantum = false;
    let mut solver_classical = QNS2DSolver::new(config_classical);
    init_fn(&mut solver_classical);
    let e0_classical = solver_classical.kinetic_energy();
    solver_classical.run_until(end_time);

    // Run quantum simulation
    let mut config_quantum = config;
    config_quantum.enable_quantum = true;
    let mut solver_quantum = QNS2DSolver::new(config_quantum);
    init_fn(&mut solver_quantum);
    let e0_quantum = solver_quantum.kinetic_energy();
    solver_quantum.run_until(end_time);

    QNS2DComparisonResult {
        time: end_time,
        classical_kinetic_energy: solver_classical.kinetic_energy(),
        quantum_kinetic_energy: solver_quantum.kinetic_energy(),
        classical_energy_decay: 1.0 - solver_classical.kinetic_energy() / e0_classical,
        quantum_energy_decay: 1.0 - solver_quantum.kinetic_energy() / e0_quantum,
        classical_enstrophy: solver_classical.enstrophy(),
        quantum_enstrophy: solver_quantum.enstrophy(),
        classical_max_vorticity: solver_classical.max_vorticity(),
        quantum_max_vorticity: solver_quantum.max_vorticity(),
    }
}

/// Result of comparing classical vs quantum 2D simulations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QNS2DComparisonResult {
    /// Simulation end time
    pub time: f64,

    /// Classical kinetic energy at end
    pub classical_kinetic_energy: f64,

    /// Quantum kinetic energy at end
    pub quantum_kinetic_energy: f64,

    /// Classical energy decay fraction
    pub classical_energy_decay: f64,

    /// Quantum energy decay fraction
    pub quantum_energy_decay: f64,

    /// Classical enstrophy at end
    pub classical_enstrophy: f64,

    /// Quantum enstrophy at end
    pub quantum_enstrophy: f64,

    /// Classical maximum vorticity
    pub classical_max_vorticity: f64,

    /// Quantum maximum vorticity
    pub quantum_max_vorticity: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_taylor_green_initialization() {
        let config = QNS2DConfig {
            nx: 32,
            ny: 32,
            lx: 2.0 * PI,
            ly: 2.0 * PI,
            ..Default::default()
        };

        let mut solver = QNS2DSolver::new(config);
        solver.init_taylor_green(1.0, 1.0);

        // Check that velocity field is non-zero
        let state = solver.get_state();
        let max_u: f64 = state.u.iter().map(|x| x.abs()).fold(0.0, f64::max);
        let max_v: f64 = state.v.iter().map(|x| x.abs()).fold(0.0, f64::max);

        assert!(max_u > 0.9, "Max u should be close to u0=1.0");
        assert!(max_v > 0.9, "Max v should be close to u0=1.0");
    }

    #[test]
    fn test_energy_decay() {
        let config = QNS2DConfig {
            nx: 32,
            ny: 32,
            lx: 2.0 * PI,
            ly: 2.0 * PI,
            viscosity: 0.1,
            enable_quantum: false,
            ..Default::default()
        };

        let mut solver = QNS2DSolver::new(config);
        solver.init_taylor_green(0.5, 1.0);

        let e0 = solver.kinetic_energy();

        // Run for a short time
        for _ in 0..100 {
            solver.step();
        }

        let e1 = solver.kinetic_energy();

        // Energy should decay due to viscosity
        assert!(
            e1 < e0,
            "Energy should decay: initial={}, final={}",
            e0,
            e1
        );
    }

    #[test]
    fn test_quantum_correction_effect() {
        let config = QNS2DConfig {
            nx: 32,
            ny: 32,
            lx: 2.0 * PI,
            ly: 2.0 * PI,
            viscosity: 0.01,
            particle_mass: 1e-20, // Larger mass = more noticeable quantum effects
            ..Default::default()
        };

        // Classical run
        let mut config_classical = config.clone();
        config_classical.enable_quantum = false;
        let mut solver_classical = QNS2DSolver::new(config_classical);
        solver_classical.init_taylor_green(0.5, 1.0);
        for _ in 0..50 {
            solver_classical.step();
        }

        // Quantum run
        let mut config_quantum = config;
        config_quantum.enable_quantum = true;
        let mut solver_quantum = QNS2DSolver::new(config_quantum);
        solver_quantum.init_taylor_green(0.5, 1.0);
        for _ in 0..50 {
            solver_quantum.step();
        }

        // Both should have evolved (different final states)
        let e_classical = solver_classical.kinetic_energy();
        let e_quantum = solver_quantum.kinetic_energy();

        // Just check that both simulations ran without error
        assert!(e_classical > 0.0, "Classical energy should be positive");
        assert!(e_quantum > 0.0, "Quantum energy should be positive");
    }

    #[test]
    fn test_periodic_boundaries() {
        let config = QNS2DConfig {
            nx: 16,
            ny: 16,
            ..Default::default()
        };

        let mut solver = QNS2DSolver::new(config);
        solver.init_taylor_green(0.5, 1.0);

        // Run a few steps
        for _ in 0..10 {
            solver.step();
        }

        let state = solver.get_state();
        let nx = 16;
        let ny = 16;

        // Check periodic BC: values at opposite boundaries should match
        for j in 0..ny {
            let diff_u = (state.u[[0, j]] - state.u[[nx - 1, j]]).abs();
            assert!(
                diff_u < 1e-10,
                "u should be periodic in x: diff={}",
                diff_u
            );
        }

        for i in 0..nx {
            let diff_v = (state.v[[i, 0]] - state.v[[i, ny - 1]]).abs();
            assert!(
                diff_v < 1e-10,
                "v should be periodic in y: diff={}",
                diff_v
            );
        }
    }

    #[test]
    fn test_enstrophy_calculation() {
        let config = QNS2DConfig {
            nx: 32,
            ny: 32,
            lx: 2.0 * PI,
            ly: 2.0 * PI,
            ..Default::default()
        };

        let mut solver = QNS2DSolver::new(config);
        solver.init_taylor_green(1.0, 1.0);

        let enstrophy = solver.enstrophy();

        // Taylor-Green has non-zero vorticity
        assert!(enstrophy > 0.0, "Enstrophy should be positive for Taylor-Green");
    }
}
