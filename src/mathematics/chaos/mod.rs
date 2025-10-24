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
mod tests {
    use super::*;

    #[test]
    fn test_complex_operations() {
        let c1 = Complex::new(1.0, 2.0);
        let c2 = Complex::new(3.0, 4.0);

        let sum = c1 + c2;
        assert_eq!(sum.re, 4.0);
        assert_eq!(sum.im, 6.0);

        let prod = c1 * c2;
        assert_eq!(prod.re, -5.0); // (1*3 - 2*4)
        assert_eq!(prod.im, 10.0);  // (1*4 + 2*3)
    }

    #[test]
    fn test_complex_magnitude() {
        let c = Complex::new(3.0, 4.0);
        assert_eq!(c.mag(), 5.0);
        assert_eq!(c.mag_sq(), 25.0);
    }

    #[test]
    fn test_point3d_distance() {
        let p1 = Point3D::new(0.0, 0.0, 0.0);
        let p2 = Point3D::new(3.0, 4.0, 0.0);
        assert_eq!(p1.distance(&p2), 5.0);
    }
}
