//! Materials Science Module
//!
//! Comprehensive materials science calculations including:
//! - Crystal structures and crystallography
//! - Band theory and electronic structure
//! - Mechanical properties (stress, strain, elasticity)
//! - Thermal properties (heat capacity, thermal expansion)
//! - Diffusion and phase transformations
//! - X-ray diffraction analysis
//! - Surface science and interfaces

pub mod crystal_structures;
pub mod band_theory;
pub mod mechanical_properties;
pub mod thermal_properties;
pub mod diffusion;
pub mod xrd_analysis;

pub use crystal_structures::*;
pub use band_theory::*;
pub use mechanical_properties::*;
pub use thermal_properties::*;
pub use diffusion::*;
pub use xrd_analysis::*;

use std::f64::consts::PI;

/// Common physical constants for materials science
pub mod constants {
    /// Boltzmann constant (J/K)
    pub const K_B: f64 = 1.380649e-23;

    /// Planck constant (J·s)
    pub const H: f64 = 6.62607015e-34;

    /// Reduced Planck constant (ℏ = h/2π) (J·s)
    pub const HBAR: f64 = 1.054571817e-34;

    /// Elementary charge (C)
    pub const E: f64 = 1.602176634e-19;

    /// Electron mass (kg)
    pub const M_E: f64 = 9.1093837015e-31;

    /// Avogadro constant (1/mol)
    pub const N_A: f64 = 6.02214076e23;

    /// Gas constant (J/(mol·K))
    pub const R: f64 = 8.314462618;

    /// Vacuum permittivity (F/m)
    pub const EPSILON_0: f64 = 8.8541878128e-12;

    /// Bohr radius (m)
    pub const A_0: f64 = 5.29177210903e-11;
}

/// Crystal system classification
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CrystalSystem {
    Cubic,
    Tetragonal,
    Orthorhombic,
    Hexagonal,
    Trigonal,
    Monoclinic,
    Triclinic,
}

/// Bravais lattice types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BravaisLattice {
    // Cubic (3)
    SimpleCubic,
    BodyCenteredCubic,
    FaceCenteredCubic,

    // Tetragonal (2)
    SimpleTetragonal,
    BodyCenteredTetragonal,

    // Orthorhombic (4)
    SimpleOrthorhombic,
    BodyCenteredOrthorhombic,
    BaseCenteredOrthorhombic,
    FaceCenteredOrthorhombic,

    // Hexagonal (1)
    Hexagonal,

    // Trigonal/Rhombohedral (1)
    Rhombohedral,

    // Monoclinic (2)
    SimpleMonoclinic,
    BaseCenteredMonoclinic,

    // Triclinic (1)
    Triclinic,
}

/// Material class for electronic properties
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialClass {
    Insulator,
    Semiconductor,
    Conductor,
    Superconductor,
}

/// Vector in 3D space (used for lattice vectors, positions)
pub type Vector3D = [f64; 3];

/// Matrix 3x3 (used for strain tensors, rotation matrices)
pub type Matrix3x3 = [[f64; 3]; 3];
