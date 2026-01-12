//! Numerical Integration Methods for N-Body Simulation

use super::{Body, IntegrationMethod, NBodyConfig, Vec3};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationResult {
    pub trajectories: Vec<Vec<Vec3>>, // [body_index][time_step]
    pub energies: Vec<f64>,            // Total energy at each step
    pub times: Vec<f64>,               // Time values
}

/// Run N-body simulation
pub fn simulate(config: &NBodyConfig) -> SimulationResult {
    match config.method {
        IntegrationMethod::Euler => simulate_euler(config),
        IntegrationMethod::Verlet => simulate_verlet(config),
        IntegrationMethod::RungeKutta4 => simulate_rk4(config),
        IntegrationMethod::BarnesHut { theta: _ } => simulate_verlet(config), // TODO: Barnes-Hut
    }
}

/// Euler method (simple but inaccurate)
fn simulate_euler(config: &NBodyConfig) -> SimulationResult {
    let mut bodies = config.bodies.clone();
    let n = bodies.len();
    let mut trajectories = vec![vec![]; n];
    let mut energies = vec![];
    let mut times = vec![];

    for _ in 0..=config.steps {
        // Record state
        let current_config = NBodyConfig {
            bodies: bodies.clone(),
            dt: config.dt,
            steps: config.steps,
            method: config.method,
            softening: config.softening,
        };
        for (i, body) in bodies.iter().enumerate() {
            trajectories[i].push(body.position);
        }
        energies.push(current_config.system_properties().total_energy);
        times.push(times.len() as f64 * config.dt);

        // Update velocities and positions
        let forces: Vec<Vec3> = (0..n)
            .map(|i| current_config.total_force(i))
            .collect();

        for (i, body) in bodies.iter_mut().enumerate() {
            let accel = forces[i] / body.mass;
            body.velocity = body.velocity + accel * config.dt;
            body.position = body.position + body.velocity * config.dt;
        }
    }

    SimulationResult { trajectories, energies, times }
}

/// Verlet integration (symplectic, energy-conserving)
fn simulate_verlet(config: &NBodyConfig) -> SimulationResult {
    let mut bodies = config.bodies.clone();
    let n = bodies.len();
    let mut trajectories = vec![vec![]; n];
    let mut energies = vec![];
    let mut times = vec![];

    // Calculate initial accelerations
    let current_config = NBodyConfig {
        bodies: bodies.clone(),
        dt: config.dt,
        steps: config.steps,
        method: config.method,
        softening: config.softening,
    };
    let mut accels: Vec<Vec3> = (0..n)
        .map(|i| current_config.total_force(i) / bodies[i].mass)
        .collect();

    for step in 0..=config.steps {
        // Record state
        let current_config = NBodyConfig {
            bodies: bodies.clone(),
            dt: config.dt,
            steps: config.steps,
            method: config.method,
            softening: config.softening,
        };
        for (i, body) in bodies.iter().enumerate() {
            trajectories[i].push(body.position);
        }
        energies.push(current_config.system_properties().total_energy);
        times.push(step as f64 * config.dt);

        // Update positions: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
        for (i, body) in bodies.iter_mut().enumerate() {
            body.position = body.position + body.velocity * config.dt + accels[i] * (0.5 * config.dt * config.dt);
        }

        // Calculate new accelerations
        let current_config = NBodyConfig {
            bodies: bodies.clone(),
            dt: config.dt,
            steps: config.steps,
            method: config.method,
            softening: config.softening,
        };
        let new_accels: Vec<Vec3> = (0..n)
            .map(|i| current_config.total_force(i) / bodies[i].mass)
            .collect();

        // Update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
        for (i, body) in bodies.iter_mut().enumerate() {
            body.velocity = body.velocity + (accels[i] + new_accels[i]) * (0.5 * config.dt);
        }

        accels = new_accels;
    }

    SimulationResult { trajectories, energies, times }
}

/// Runge-Kutta 4th order
fn simulate_rk4(config: &NBodyConfig) -> SimulationResult {
    let mut bodies = config.bodies.clone();
    let n = bodies.len();
    let mut trajectories = vec![vec![]; n];
    let mut energies = vec![];
    let mut times = vec![];

    for step in 0..=config.steps {
        // Record state
        let current_config = NBodyConfig {
            bodies: bodies.clone(),
            dt: config.dt,
            steps: config.steps,
            method: config.method,
            softening: config.softening,
        };
        for (i, body) in bodies.iter().enumerate() {
            trajectories[i].push(body.position);
        }
        energies.push(current_config.system_properties().total_energy);
        times.push(step as f64 * config.dt);

        // RK4 step
        bodies = rk4_step(&bodies, config.dt, config.softening);
    }

    SimulationResult { trajectories, energies, times }
}

fn rk4_step(bodies: &[Body], dt: f64, softening: f64) -> Vec<Body> {
    let n = bodies.len();

    // k1 = f(y)
    let k1_v: Vec<Vec3> = (0..n).map(|i| {
        let config = NBodyConfig {
            bodies: bodies.to_vec(),
            dt,
            steps: 1,
            method: IntegrationMethod::RungeKutta4,
            softening,
        };
        config.total_force(i) / bodies[i].mass
    }).collect();
    let k1_x: Vec<Vec3> = bodies.iter().map(|b| b.velocity).collect();

    // k2 = f(y + 0.5*dt*k1)
    let bodies_k2: Vec<Body> = bodies.iter().enumerate().map(|(i, b)| {
        Body {
            mass: b.mass,
            position: b.position + k1_x[i] * (0.5 * dt),
            velocity: b.velocity + k1_v[i] * (0.5 * dt),
            name: b.name.clone(),
        }
    }).collect();
    let k2_v: Vec<Vec3> = (0..n).map(|i| {
        let config = NBodyConfig {
            bodies: bodies_k2.clone(),
            dt,
            steps: 1,
            method: IntegrationMethod::RungeKutta4,
            softening,
        };
        config.total_force(i) / bodies_k2[i].mass
    }).collect();
    let k2_x: Vec<Vec3> = bodies_k2.iter().map(|b| b.velocity).collect();

    // k3 = f(y + 0.5*dt*k2)
    let bodies_k3: Vec<Body> = bodies.iter().enumerate().map(|(i, b)| {
        Body {
            mass: b.mass,
            position: b.position + k2_x[i] * (0.5 * dt),
            velocity: b.velocity + k2_v[i] * (0.5 * dt),
            name: b.name.clone(),
        }
    }).collect();
    let k3_v: Vec<Vec3> = (0..n).map(|i| {
        let config = NBodyConfig {
            bodies: bodies_k3.clone(),
            dt,
            steps: 1,
            method: IntegrationMethod::RungeKutta4,
            softening,
        };
        config.total_force(i) / bodies_k3[i].mass
    }).collect();
    let k3_x: Vec<Vec3> = bodies_k3.iter().map(|b| b.velocity).collect();

    // k4 = f(y + dt*k3)
    let bodies_k4: Vec<Body> = bodies.iter().enumerate().map(|(i, b)| {
        Body {
            mass: b.mass,
            position: b.position + k3_x[i] * dt,
            velocity: b.velocity + k3_v[i] * dt,
            name: b.name.clone(),
        }
    }).collect();
    let k4_v: Vec<Vec3> = (0..n).map(|i| {
        let config = NBodyConfig {
            bodies: bodies_k4.clone(),
            dt,
            steps: 1,
            method: IntegrationMethod::RungeKutta4,
            softening,
        };
        config.total_force(i) / bodies_k4[i].mass
    }).collect();
    let k4_x: Vec<Vec3> = bodies_k4.iter().map(|b| b.velocity).collect();

    // y(t+dt) = y(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    bodies.iter().enumerate().map(|(i, b)| {
        Body {
            mass: b.mass,
            position: b.position + (k1_x[i] + k2_x[i] * 2.0 + k3_x[i] * 2.0 + k4_x[i]) * (dt / 6.0),
            velocity: b.velocity + (k1_v[i] + k2_v[i] * 2.0 + k3_v[i] * 2.0 + k4_v[i]) * (dt / 6.0),
            name: b.name.clone(),
        }
    }).collect()
}

