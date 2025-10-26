//! Chaos Theory and Fractals
//!
//! **Chaos Theory:**
//! - Lorenz attractor (weather prediction)
//! - Lyapunov exponents (measure of chaos)
//! - Strange attractors
//! - Bifurcation diagrams
//!
//! **Fractals:**
//! - Mandelbrot set
//! - Julia sets
//! - Fractal dimension (box-counting)
//! - Iterated function systems

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub mod fractals;
pub mod attractors;
pub mod lyapunov;

pub use fractals::*;
pub use attractors::*;
pub use lyapunov::*;

/// Complex number for fractal calculations
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }

    pub fn mag_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    pub fn mag(&self) -> f64 {
        self.mag_sq().sqrt()
    }
}

impl std::ops::Add for Complex {
    type Output = Complex;
    fn add(self, other: Complex) -> Complex {
        Complex { re: self.re + other.re, im: self.im + other.im }
    }
}

impl std::ops::Mul for Complex {
    type Output = Complex;
    fn mul(self, other: Complex) -> Complex {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

/// 3D point for attractor visualization
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn distance(&self, other: &Point3D) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/mathematics/chaos_tests.rs"]
mod tests;

