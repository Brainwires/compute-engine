//! Crystal Structures and Crystallography
//!
//! Calculations for crystal lattices, unit cells, Miller indices,
//! interplanar spacing, and crystallographic analysis.

use super::{BravaisLattice, CrystalSystem, Vector3D};
use std::f64::consts::PI;

/// Unit cell parameters
#[derive(Debug, Clone)]
pub struct UnitCell {
    /// Lattice parameter a (Å)
    pub a: f64,
    /// Lattice parameter b (Å)
    pub b: f64,
    /// Lattice parameter c (Å)
    pub c: f64,
    /// Angle α (degrees)
    pub alpha: f64,
    /// Angle β (degrees)
    pub beta: f64,
    /// Angle γ (degrees)
    pub gamma: f64,
}

impl UnitCell {
    /// Create cubic unit cell
    pub fn cubic(a: f64) -> Self {
        Self {
            a,
            b: a,
            c: a,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }

    /// Create tetragonal unit cell
    pub fn tetragonal(a: f64, c: f64) -> Self {
        Self {
            a,
            b: a,
            c,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }

    /// Create orthorhombic unit cell
    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Self {
        Self {
            a,
            b,
            c,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }

    /// Create hexagonal unit cell
    pub fn hexagonal(a: f64, c: f64) -> Self {
        Self {
            a,
            b: a,
            c,
            alpha: 90.0,
            beta: 90.0,
            gamma: 120.0,
        }
    }
}

/// Calculate volume of unit cell
///
/// V = abc√(1 - cos²α - cos²β - cos²γ + 2cosαcosβcosγ)
///
/// # Example
/// ```
/// use computational_engine::materials_science::*;
/// let cell = UnitCell::cubic(3.5); // 3.5 Å
/// let volume = unit_cell_volume(&cell);
/// assert!((volume - 42.875).abs() < 0.001);
/// ```
pub fn unit_cell_volume(cell: &UnitCell) -> f64 {
    let alpha_rad = cell.alpha.to_radians();
    let beta_rad = cell.beta.to_radians();
    let gamma_rad = cell.gamma.to_radians();

    let cos_alpha = alpha_rad.cos();
    let cos_beta = beta_rad.cos();
    let cos_gamma = gamma_rad.cos();

    let discriminant = 1.0 - cos_alpha.powi(2) - cos_beta.powi(2) - cos_gamma.powi(2)
        + 2.0 * cos_alpha * cos_beta * cos_gamma;

    cell.a * cell.b * cell.c * discriminant.sqrt()
}

/// Calculate interplanar spacing d_{hkl} for cubic crystals
///
/// d = a / √(h² + k² + l²)
///
/// # Arguments
/// * `a` - Lattice parameter (Å)
/// * `h`, `k`, `l` - Miller indices
///
/// # Returns
/// Interplanar spacing (Å)
pub fn interplanar_spacing_cubic(a: f64, h: i32, k: i32, l: i32) -> f64 {
    a / ((h * h + k * k + l * l) as f64).sqrt()
}

/// Calculate interplanar spacing for tetragonal crystals
///
/// 1/d² = (h² + k²)/a² + l²/c²
pub fn interplanar_spacing_tetragonal(a: f64, c: f64, h: i32, k: i32, l: i32) -> f64 {
    let inv_d_sq = ((h * h + k * k) as f64) / (a * a) + (l * l) as f64 / (c * c);
    1.0 / inv_d_sq.sqrt()
}

/// Calculate interplanar spacing for orthorhombic crystals
///
/// 1/d² = h²/a² + k²/b² + l²/c²
pub fn interplanar_spacing_orthorhombic(
    a: f64,
    b: f64,
    c: f64,
    h: i32,
    k: i32,
    l: i32,
) -> f64 {
    let inv_d_sq =
        (h * h) as f64 / (a * a) + (k * k) as f64 / (b * b) + (l * l) as f64 / (c * c);
    1.0 / inv_d_sq.sqrt()
}

/// Calculate interplanar spacing for hexagonal crystals
///
/// 1/d² = 4/3 * (h² + hk + k²)/a² + l²/c²
pub fn interplanar_spacing_hexagonal(a: f64, c: f64, h: i32, k: i32, l: i32) -> f64 {
    let inv_d_sq = (4.0 / 3.0) * ((h * h + h * k + k * k) as f64) / (a * a)
        + (l * l) as f64 / (c * c);
    1.0 / inv_d_sq.sqrt()
}

/// Calculate Bragg angle for X-ray diffraction
///
/// 2d sinθ = nλ → θ = arcsin(nλ/(2d))
///
/// # Arguments
/// * `d_spacing` - Interplanar spacing (Å)
/// * `wavelength` - X-ray wavelength (Å), typically Cu Kα = 1.5418 Å
/// * `n` - Order of diffraction (usually 1)
///
/// # Returns
/// Bragg angle θ in degrees, or None if condition cannot be satisfied
pub fn bragg_angle(d_spacing: f64, wavelength: f64, n: u32) -> Option<f64> {
    let arg = (n as f64) * wavelength / (2.0 * d_spacing);
    if arg > 1.0 {
        None // Cannot satisfy Bragg condition
    } else {
        Some(arg.asin().to_degrees())
    }
}

/// Calculate 2θ angle for XRD pattern
pub fn two_theta(d_spacing: f64, wavelength: f64, n: u32) -> Option<f64> {
    bragg_angle(d_spacing, wavelength, n).map(|theta| 2.0 * theta)
}

/// Calculate atomic packing factor for different structures
///
/// APF = (Volume of atoms in unit cell) / (Volume of unit cell)
pub fn atomic_packing_factor(lattice: BravaisLattice) -> f64 {
    match lattice {
        BravaisLattice::SimpleCubic => PI / 6.0,                   // 0.524
        BravaisLattice::BodyCenteredCubic => 3.0_f64.sqrt() * PI / 8.0, // 0.680
        BravaisLattice::FaceCenteredCubic => PI * 2.0_f64.sqrt() / 6.0, // 0.740
        BravaisLattice::Hexagonal => PI * 2.0_f64.sqrt() / 6.0,  // 0.740 (ideal c/a)
        _ => 0.0, // Other structures need more parameters
    }
}

/// Calculate coordination number for common structures
pub fn coordination_number(lattice: BravaisLattice) -> u32 {
    match lattice {
        BravaisLattice::SimpleCubic => 6,
        BravaisLattice::BodyCenteredCubic => 8,
        BravaisLattice::FaceCenteredCubic => 12,
        BravaisLattice::Hexagonal => 12,
        _ => 0,
    }
}

/// Calculate number of atoms per unit cell
pub fn atoms_per_unit_cell(lattice: BravaisLattice) -> f64 {
    match lattice {
        BravaisLattice::SimpleCubic => 1.0,
        BravaisLattice::BodyCenteredCubic => 2.0,
        BravaisLattice::FaceCenteredCubic => 4.0,
        BravaisLattice::Hexagonal => 2.0,
        BravaisLattice::SimpleTetragonal => 1.0,
        BravaisLattice::BodyCenteredTetragonal => 2.0,
        _ => 0.0,
    }
}

/// Calculate nearest neighbor distance for cubic structures
///
/// # Arguments
/// * `lattice_param` - Lattice parameter a (Å)
/// * `lattice` - Bravais lattice type
pub fn nearest_neighbor_distance(lattice_param: f64, lattice: BravaisLattice) -> f64 {
    match lattice {
        BravaisLattice::SimpleCubic => lattice_param,
        BravaisLattice::BodyCenteredCubic => 3.0_f64.sqrt() * lattice_param / 2.0,
        BravaisLattice::FaceCenteredCubic => 2.0_f64.sqrt() * lattice_param / 2.0,
        _ => lattice_param, // Simplified
    }
}

/// Calculate density of a crystal
///
/// ρ = (n * M) / (N_A * V)
///
/// # Arguments
/// * `atoms_per_cell` - Number of atoms per unit cell
/// * `atomic_mass` - Atomic mass (g/mol)
/// * `cell_volume` - Unit cell volume (Å³)
///
/// # Returns
/// Density in g/cm³
pub fn crystal_density(atoms_per_cell: f64, atomic_mass: f64, cell_volume_angstrom3: f64) -> f64 {
    const N_A: f64 = 6.02214076e23;
    let cell_volume_cm3 = cell_volume_angstrom3 * 1e-24; // Å³ to cm³
    (atoms_per_cell * atomic_mass) / (N_A * cell_volume_cm3)
}

/// Calculate reciprocal lattice vector magnitude
///
/// For cubic: b* = 2π/a
pub fn reciprocal_lattice_cubic(a: f64) -> f64 {
    2.0 * PI / a
}

/// Calculate structure factor for simple cubic (real valued for simplicity)
///
/// F_{hkl} = Σ f_j exp[2πi(hx_j + ky_j + lz_j)]
///
/// For simple structures, returns magnitude
pub fn structure_factor_simple_cubic(h: i32, k: i32, l: i32) -> f64 {
    // For simple cubic, all reflections allowed
    // |F|² = f² for all hkl
    1.0
}

/// Check systematic absences for BCC
///
/// BCC: Reflections present only when h + k + l = even
pub fn is_reflection_allowed_bcc(h: i32, k: i32, l: i32) -> bool {
    (h + k + l) % 2 == 0
}

/// Check systematic absences for FCC
///
/// FCC: Reflections present only when h, k, l are all even or all odd
pub fn is_reflection_allowed_fcc(h: i32, k: i32, l: i32) -> bool {
    let h_even = h % 2 == 0;
    let k_even = k % 2 == 0;
    let l_even = l % 2 == 0;

    (h_even && k_even && l_even) || (!h_even && !k_even && !l_even)
}

/// Convert Miller indices to direction vector
pub fn miller_to_direction(h: i32, k: i32, l: i32, cell: &UnitCell) -> Vector3D {
    [h as f64 * cell.a, k as f64 * cell.b, l as f64 * cell.c]
}

/// Calculate angle between two crystallographic planes
///
/// For cubic systems: cos(φ) = (h₁h₂ + k₁k₂ + l₁l₂) / √[(h₁²+k₁²+l₁²)(h₂²+k₂²+l₂²)]
pub fn angle_between_planes_cubic(
    h1: i32,
    k1: i32,
    l1: i32,
    h2: i32,
    k2: i32,
    l2: i32,
) -> f64 {
    let dot = (h1 * h2 + k1 * k2 + l1 * l2) as f64;
    let mag1 = ((h1 * h1 + k1 * k1 + l1 * l1) as f64).sqrt();
    let mag2 = ((h2 * h2 + k2 * k2 + l2 * l2) as f64).sqrt();

    (dot / (mag1 * mag2)).acos().to_degrees()
}

/// Calculate linear density (atoms per unit length) along a direction
///
/// For [100] in SC: ρ_L = 1/a
/// For [110] in SC: ρ_L = 1/(√2 * a)
/// For [111] in SC: ρ_L = 1/(√3 * a)
pub fn linear_density_cubic(a: f64, h: i32, k: i32, l: i32) -> f64 {
    let direction_length = a * ((h * h + k * k + l * l) as f64).sqrt();
    // For simple cubic along any direction
    1.0 / direction_length
}

/// Calculate planar density (atoms per unit area) for a plane
///
/// Returns atoms/Ų
pub fn planar_density_cubic_100(a: f64) -> f64 {
    // For (100) plane in simple cubic: 1 atom per a² area
    1.0 / (a * a)
}

pub fn planar_density_cubic_110(a: f64) -> f64 {
    // For (110) plane in simple cubic
    1.0 / (2.0_f64.sqrt() * a * a)
}

pub fn planar_density_cubic_111(a: f64) -> f64 {
    // For (111) plane in simple cubic
    1.0 / (3.0_f64.sqrt() * a * a)
}

