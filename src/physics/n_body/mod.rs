//! N-Body Gravitational Simulation
//!
//! Simulates gravitational interactions between N bodies using numerical integration.
//!
//! **Methods:**
//! - Euler method (fast, low accuracy)
//! - Verlet integration (energy-conserving)
//! - Runge-Kutta 4 (higher accuracy)
//! - Barnes-Hut algorithm (O(N log N) for large systems)
//!
//! **Applications:**
//! - Solar system evolution
//! - Galaxy collisions
//! - Star cluster dynamics
//! - Orbital mechanics

use serde::{Deserialize, Serialize};

pub mod integrator;
pub mod barnes_hut;

pub use integrator::*;
pub use barnes_hut::*;

/// Physical constants
pub const G: f64 = 6.67430e-11; // m³/(kg·s²)

/// 3D vector
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0 }
    }

    pub fn mag(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn mag_sq(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn normalize(&self) -> Self {
        let m = self.mag();
        if m > 1e-10 {
            Self { x: self.x / m, y: self.y / m, z: self.z / m }
        } else {
            Self::zero()
        }
    }

    pub fn dot(&self, other: &Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}

impl std::ops::Add for Vec3 {
    type Output = Vec3;
    fn add(self, other: Vec3) -> Vec3 {
        Vec3 { x: self.x + other.x, y: self.y + other.y, z: self.z + other.z }
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, other: Vec3) -> Vec3 {
        Vec3 { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

impl std::ops::Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, scalar: f64) -> Vec3 {
        Vec3 { x: self.x * scalar, y: self.y * scalar, z: self.z * scalar }
    }
}

impl std::ops::Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, scalar: f64) -> Vec3 {
        Vec3 { x: self.x / scalar, y: self.y / scalar, z: self.z / scalar }
    }
}

/// Body in N-body simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Body {
    pub mass: f64,
    pub position: Vec3,
    pub velocity: Vec3,
    pub name: String,
}

impl Body {
    pub fn new(mass: f64, position: Vec3, velocity: Vec3, name: &str) -> Self {
        Self { mass, position, velocity, name: name.to_string() }
    }

    /// Calculate kinetic energy: KE = 0.5 * m * v²
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity.mag_sq()
    }

    /// Calculate momentum: p = m * v
    pub fn momentum(&self) -> Vec3 {
        self.velocity * self.mass
    }

    /// Calculate angular momentum: L = r × p
    pub fn angular_momentum(&self) -> Vec3 {
        self.position.cross(&self.momentum())
    }
}

/// N-body simulation configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NBodyConfig {
    pub bodies: Vec<Body>,
    pub dt: f64,            // Time step (seconds)
    pub steps: usize,       // Number of steps
    pub method: IntegrationMethod,
    pub softening: f64,     // Softening length (to avoid singularities)
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum IntegrationMethod {
    Euler,
    Verlet,
    RungeKutta4,
    BarnesHut { theta: f64 }, // Opening angle for tree
}

/// System-wide properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemProperties {
    pub total_energy: f64,
    pub kinetic_energy: f64,
    pub potential_energy: f64,
    pub total_momentum: Vec3,
    pub total_angular_momentum: Vec3,
    pub center_of_mass: Vec3,
}

impl NBodyConfig {
    /// Calculate gravitational force on body i from body j
    pub fn force_between(&self, i: usize, j: usize) -> Vec3 {
        let bi = &self.bodies[i];
        let bj = &self.bodies[j];

        let r = bj.position - bi.position;
        let dist_sq = r.mag_sq() + self.softening * self.softening;
        let dist = dist_sq.sqrt();

        // F = G * m1 * m2 / r² * r̂
        let f_mag = G * bi.mass * bj.mass / dist_sq;
        r.normalize() * f_mag
    }

    /// Calculate total gravitational force on body i
    pub fn total_force(&self, i: usize) -> Vec3 {
        let mut force = Vec3::zero();
        for j in 0..self.bodies.len() {
            if i != j {
                force = force + self.force_between(i, j);
            }
        }
        force
    }

    /// Calculate gravitational potential energy of the system
    pub fn potential_energy(&self) -> f64 {
        let mut pe = 0.0;
        for i in 0..self.bodies.len() {
            for j in (i+1)..self.bodies.len() {
                let r = (self.bodies[i].position - self.bodies[j].position).mag();
                if r > 1e-10 {
                    pe -= G * self.bodies[i].mass * self.bodies[j].mass / r;
                }
            }
        }
        pe
    }

    /// Calculate system properties
    pub fn system_properties(&self) -> SystemProperties {
        let ke: f64 = self.bodies.iter().map(|b| b.kinetic_energy()).sum();
        let pe = self.potential_energy();

        let total_momentum: Vec3 = self.bodies.iter()
            .map(|b| b.momentum())
            .fold(Vec3::zero(), |acc, p| acc + p);

        let total_angular_momentum: Vec3 = self.bodies.iter()
            .map(|b| b.angular_momentum())
            .fold(Vec3::zero(), |acc, l| acc + l);

        let total_mass: f64 = self.bodies.iter().map(|b| b.mass).sum();
        let center_of_mass = self.bodies.iter()
            .map(|b| b.position * b.mass)
            .fold(Vec3::zero(), |acc, p| acc + p) / total_mass;

        SystemProperties {
            total_energy: ke + pe,
            kinetic_energy: ke,
            potential_energy: pe,
            total_momentum,
            total_angular_momentum,
            center_of_mass,
        }
    }
}

