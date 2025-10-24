//! X-Ray Diffraction Analysis
//!
//! XRD pattern simulation, peak analysis, lattice parameter refinement,
//! crystallite size determination, and phase identification.

use super::constants::*;
use super::{BravaisLattice, UnitCell};
use std::f64::consts::PI;

/// Standard X-ray wavelengths (Å)
pub mod wavelengths {
    /// Cu Kα (most common)
    pub const CU_KA: f64 = 1.5418;
    /// Cu Kα1
    pub const CU_KA1: f64 = 1.54056;
    /// Cu Kα2
    pub const CU_KA2: f64 = 1.54439;
    /// Mo Kα
    pub const MO_KA: f64 = 0.71073;
    /// Co Kα
    pub const CO_KA: f64 = 1.79021;
    /// Fe Kα
    pub const FE_KA: f64 = 1.93604;
    /// Cr Kα
    pub const CR_KA: f64 = 2.29100;
}

/// XRD peak information
#[derive(Debug, Clone)]
pub struct XRDPeak {
    /// Miller indices (h, k, l)
    pub hkl: (i32, i32, i32),
    /// 2θ angle (degrees)
    pub two_theta: f64,
    /// d-spacing (Å)
    pub d_spacing: f64,
    /// Relative intensity (0-100)
    pub intensity: f64,
}

/// Calculate 2θ for a given d-spacing and wavelength
///
/// λ = 2d sinθ → θ = arcsin(λ/(2d))
///
/// # Returns
/// 2θ in degrees, or None if Bragg condition cannot be satisfied
pub fn calculate_two_theta(d_spacing: f64, wavelength: f64) -> Option<f64> {
    let sin_theta = wavelength / (2.0 * d_spacing);
    if sin_theta > 1.0 {
        None
    } else {
        Some(2.0 * sin_theta.asin().to_degrees())
    }
}

/// Calculate d-spacing from 2θ angle
///
/// d = λ/(2 sinθ)
pub fn calculate_d_spacing(two_theta_deg: f64, wavelength: f64) -> f64 {
    let theta_rad = (two_theta_deg / 2.0).to_radians();
    wavelength / (2.0 * theta_rad.sin())
}

/// Calculate structure factor magnitude for simple structures
///
/// For more complex structures, would need atomic form factors
///
/// Returns relative intensity based on multiplicity and structure factor
pub fn calculate_intensity_simple(
    hkl: (i32, i32, i32),
    lattice: BravaisLattice,
) -> f64 {
    let (h, k, l) = hkl;

    // Check systematic absences
    let allowed = match lattice {
        BravaisLattice::SimpleCubic => true,
        BravaisLattice::BodyCenteredCubic => (h + k + l) % 2 == 0,
        BravaisLattice::FaceCenteredCubic => {
            let h_even = h % 2 == 0;
            let k_even = k % 2 == 0;
            let l_even = l % 2 == 0;
            (h_even && k_even && l_even) || (!h_even && !k_even && !l_even)
        }
        _ => true, // Simplified
    };

    if !allowed {
        return 0.0;
    }

    // Calculate multiplicity (number of equivalent planes)
    let multiplicity = calculate_multiplicity(hkl);

    // Simplified: intensity proportional to multiplicity
    // Real calculation would include Lorentz-polarization factor,
    // atomic form factors, temperature factor, etc.
    multiplicity as f64
}

/// Calculate multiplicity of an hkl reflection (for cubic systems)
///
/// Number of equivalent planes due to crystal symmetry
fn calculate_multiplicity(hkl: (i32, i32, i32)) -> u32 {
    let (h, k, l) = hkl;
    let h_abs = h.abs();
    let k_abs = k.abs();
    let l_abs = l.abs();

    // Count zeros and distinct non-zero values
    let num_zeros = (if h_abs == 0 { 1 } else { 0 })
        + (if k_abs == 0 { 1 } else { 0 })
        + (if l_abs == 0 { 1 } else { 0 });

    // All three the same
    if h_abs == k_abs && k_abs == l_abs {
        if h_abs == 0 {
            1 // (000) - not actually a reflection
        } else {
            8 // {111} type: all same non-zero
        }
    }
    // Exactly two values the same
    else if (h_abs == k_abs && l_abs != h_abs)
        || (k_abs == l_abs && h_abs != k_abs)
        || (h_abs == l_abs && k_abs != h_abs)
    {
        // Check if the two that are the same are BOTH zero (like 0,0,1) vs both non-zero (like 1,1,0)
        if num_zeros >= 2 {
            6 // {100} type: two zeros, one non-zero like (1,0,0)
        } else if num_zeros == 0 {
            24 // {110} type: two same non-zero, one different non-zero like (1,1,2)
        } else {
            // num_zeros == 1: two non-zero are the same, one is zero like (1,1,0)
            24 // {110} type
        }
    }
    // All different
    else {
        if num_zeros == 0 {
            48 // {123} type: all different, all non-zero
        } else {
            12 // {210} type: all different, one zero
        }
    }
}

/// Generate XRD pattern for a cubic crystal
///
/// # Arguments
/// * `lattice_param` - Lattice parameter a (Å)
/// * `lattice` - Bravais lattice type
/// * `wavelength` - X-ray wavelength (Å)
/// * `max_two_theta` - Maximum 2θ to calculate (degrees)
///
/// # Returns
/// Vector of XRD peaks
pub fn generate_cubic_pattern(
    lattice_param: f64,
    lattice: BravaisLattice,
    wavelength: f64,
    max_two_theta: f64,
) -> Vec<XRDPeak> {
    let mut peaks = Vec::new();

    // Generate hkl indices up to reasonable values
    for h in 0..=6 {
        for k in 0..=h {
            for l in 0..=k {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }

                // Calculate d-spacing
                let d = lattice_param / ((h * h + k * k + l * l) as f64).sqrt();

                // Calculate 2θ
                if let Some(two_theta) = calculate_two_theta(d, wavelength) {
                    if two_theta <= max_two_theta {
                        let intensity = calculate_intensity_simple((h, k, l), lattice);

                        if intensity > 0.0 {
                            peaks.push(XRDPeak {
                                hkl: (h, k, l),
                                two_theta,
                                d_spacing: d,
                                intensity,
                            });
                        }
                    }
                }
            }
        }
    }

    // Normalize intensities to 0-100 scale
    if !peaks.is_empty() {
        let max_intensity = peaks.iter().map(|p| p.intensity).fold(0.0, f64::max);
        for peak in &mut peaks {
            peak.intensity = (peak.intensity / max_intensity) * 100.0;
        }
    }

    // Sort by 2θ
    peaks.sort_by(|a, b| a.two_theta.partial_cmp(&b.two_theta).unwrap());

    peaks
}

/// Calculate crystallite size from peak broadening (Scherrer equation)
///
/// τ = (Kλ)/(β cosθ)
///
/// # Arguments
/// * `wavelength` - X-ray wavelength (Å)
/// * `fwhm_radians` - Full width at half maximum in radians (2θ scale)
/// * `theta_radians` - Bragg angle θ in radians
/// * `k_factor` - Shape factor K (typically 0.9 for spherical crystals)
///
/// # Returns
/// Crystallite size in Å
pub fn scherrer_crystallite_size(
    wavelength: f64,
    fwhm_radians: f64,
    theta_radians: f64,
    k_factor: f64,
) -> f64 {
    (k_factor * wavelength) / (fwhm_radians * theta_radians.cos())
}

/// Convert peak width from degrees to radians
pub fn fwhm_deg_to_rad(fwhm_degrees: f64) -> f64 {
    fwhm_degrees.to_radians()
}

/// Calculate instrumental broadening correction
///
/// β_sample² = β_observed² - β_instrument²
///
/// # Returns
/// Corrected FWHM in same units as input
pub fn correct_instrumental_broadening(
    observed_fwhm: f64,
    instrumental_fwhm: f64,
) -> f64 {
    (observed_fwhm.powi(2) - instrumental_fwhm.powi(2)).sqrt()
}

/// Calculate microstrain from peak broadening
///
/// ε = β/(4 tanθ)
///
/// # Arguments
/// * `fwhm_radians` - FWHM in radians (2θ scale)
/// * `theta_radians` - Bragg angle θ
///
/// # Returns
/// Microstrain (dimensionless)
pub fn calculate_microstrain(fwhm_radians: f64, theta_radians: f64) -> f64 {
    fwhm_radians / (4.0 * theta_radians.tan())
}

/// Williamson-Hall plot: separate size and strain broadening
///
/// β cosθ = Kλ/τ + 4ε sinθ
///
/// Returns (size_contribution, strain_contribution)
///
/// # Arguments
/// * `fwhm_radians` - Peak FWHM in radians
/// * `theta_radians` - Bragg angle θ
/// * `wavelength` - X-ray wavelength
/// * `k_factor` - Scherrer constant (typically 0.9)
///
/// # Returns
/// (y-axis value for W-H plot, x-axis value for W-H plot)
pub fn williamson_hall_values(
    fwhm_radians: f64,
    theta_radians: f64,
) -> (f64, f64) {
    let y = fwhm_radians * theta_radians.cos(); // β cosθ
    let x = 4.0 * theta_radians.sin(); // 4 sinθ
    (y, x)
}

/// Calculate lattice parameter from d-spacing and Miller indices (cubic)
///
/// a = d√(h² + k² + l²)
pub fn lattice_parameter_from_d_cubic(
    d_spacing: f64,
    h: i32,
    k: i32,
    l: i32,
) -> f64 {
    d_spacing * ((h * h + k * k + l * l) as f64).sqrt()
}

/// Refine lattice parameter from multiple peaks (cubic)
///
/// Returns weighted average lattice parameter
pub fn refine_lattice_parameter_cubic(peaks: &[(f64, i32, i32, i32)]) -> f64 {
    let mut sum_a = 0.0;
    let mut sum_weight = 0.0;

    for &(d, h, k, l) in peaks {
        let a = lattice_parameter_from_d_cubic(d, h, k, l);
        let weight = ((h * h + k * k + l * l) as f64).sqrt(); // Higher indices = better precision
        sum_a += a * weight;
        sum_weight += weight;
    }

    sum_a / sum_weight
}

/// Calculate texture coefficient (for preferred orientation analysis)
///
/// TC(hkl) = [I(hkl)/I₀(hkl)] / [(1/N)Σ(I(hkl)/I₀(hkl))]
///
/// # Arguments
/// * `measured_intensity` - Measured peak intensity
/// * `reference_intensity` - Reference intensity (powder pattern)
/// * `avg_ratio` - Average of all I/I₀ ratios in the pattern
///
/// # Returns
/// Texture coefficient (1 = random, >1 = preferred orientation)
pub fn texture_coefficient(
    measured_intensity: f64,
    reference_intensity: f64,
    avg_ratio: f64,
) -> f64 {
    (measured_intensity / reference_intensity) / avg_ratio
}

/// Calculate degree of crystallinity (for semi-crystalline materials)
///
/// X_c = I_crystalline / (I_crystalline + I_amorphous)
///
/// # Arguments
/// * `crystalline_area` - Integrated area of crystalline peaks
/// * `amorphous_area` - Integrated area of amorphous halo
pub fn degree_of_crystallinity(crystalline_area: f64, amorphous_area: f64) -> f64 {
    crystalline_area / (crystalline_area + amorphous_area)
}

/// Calculate Lorentz-polarization factor
///
/// LP = (1 + cos²(2θ)) / (sin²θ cosθ)
///
/// Used to correct intensities in XRD patterns
pub fn lorentz_polarization_factor(two_theta_rad: f64) -> f64 {
    let theta = two_theta_rad / 2.0;
    let cos_2theta = (2.0 * theta).cos();
    let sin_theta = theta.sin();
    let cos_theta = theta.cos();

    (1.0 + cos_2theta * cos_2theta) / (sin_theta * sin_theta * cos_theta)
}

/// Calculate temperature factor (Debye-Waller factor)
///
/// T = exp(-2B(sinθ/λ)²)
///
/// # Arguments
/// * `b_factor` - Atomic displacement parameter (Ų)
/// * `theta_radians` - Bragg angle θ
/// * `wavelength` - X-ray wavelength (Å)
pub fn temperature_factor(b_factor: f64, theta_radians: f64, wavelength: f64) -> f64 {
    let sin_theta_over_lambda = theta_radians.sin() / wavelength;
    (-2.0 * b_factor * sin_theta_over_lambda.powi(2)).exp()
}

/// Calculate absorption correction (for flat plate geometry)
///
/// A = 1/(2μt)·[1 - exp(-2μt/sinθ)]
///
/// # Arguments
/// * `linear_absorption` - Linear absorption coefficient μ (1/cm)
/// * `sample_thickness` - Sample thickness t (cm)
/// * `theta_radians` - Bragg angle θ
pub fn absorption_correction(
    linear_absorption: f64,
    sample_thickness: f64,
    theta_radians: f64,
) -> f64 {
    let mu_t = linear_absorption * sample_thickness;
    let exponent = -2.0 * mu_t / theta_radians.sin();
    (1.0 - exponent.exp()) / (2.0 * mu_t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_theta_calculation() {
        // Silicon (111): d = 3.135 Å with Cu Kα
        let two_theta = calculate_two_theta(3.135, wavelengths::CU_KA).unwrap();
        assert!((two_theta - 28.44).abs() < 0.1); // Should be ~28.44°
    }

    #[test]
    fn test_d_spacing_from_two_theta() {
        // Reverse calculation
        let d = calculate_d_spacing(28.44, wavelengths::CU_KA);
        assert!((d - 3.135).abs() < 0.01);
    }

    #[test]
    fn test_scherrer_equation() {
        // Typical values: λ = 1.54 Å, FWHM = 0.5°, θ = 15°
        let wavelength = 1.54;
        let fwhm_rad = 0.5_f64.to_radians();
        let theta_rad = 15.0_f64.to_radians();

        let size = scherrer_crystallite_size(wavelength, fwhm_rad, theta_rad, 0.9);
        assert!(size > 100.0 && size < 200.0); // Should be ~150 Å
    }

    #[test]
    fn test_lattice_parameter_cubic() {
        // Silicon: FCC, a = 5.431 Å
        // (111) reflection: d = a/√3
        let d = 5.431 / 3.0_f64.sqrt();
        let a = lattice_parameter_from_d_cubic(d, 1, 1, 1);
        assert!((a - 5.431).abs() < 0.001);
    }

    #[test]
    fn test_generate_cubic_pattern() {
        // Silicon FCC with Cu Kα
        let pattern = generate_cubic_pattern(
            5.431,
            BravaisLattice::FaceCenteredCubic,
            wavelengths::CU_KA,
            90.0,
        );

        // Should have several peaks
        assert!(pattern.len() > 5);

        // First peak should be (111) around 28°
        assert!((pattern[0].two_theta - 28.44).abs() < 1.0);
        assert_eq!(pattern[0].hkl, (1, 1, 1));
    }

    #[test]
    fn test_multiplicity() {
        assert_eq!(calculate_multiplicity((1, 0, 0)), 6); // {100}
        assert_eq!(calculate_multiplicity((1, 1, 0)), 24); // {110}
        assert_eq!(calculate_multiplicity((1, 1, 1)), 8); // {111}
        assert_eq!(calculate_multiplicity((2, 1, 0)), 12); // {210} - has zero
    }

    #[test]
    fn test_instrumental_correction() {
        let observed = 0.5; // degrees
        let instrumental = 0.1; // degrees
        let corrected = correct_instrumental_broadening(observed, instrumental);

        assert!((corrected - 0.4899).abs() < 0.001); // √(0.25 - 0.01) ≈ 0.4899
    }

    #[test]
    fn test_lorentz_polarization() {
        // At 2θ = 30°
        let lp = lorentz_polarization_factor(30.0_f64.to_radians());
        assert!(lp > 25.0 && lp < 30.0); // ~27 for 2θ=30°
    }

    #[test]
    fn test_degree_of_crystallinity() {
        let xc = degree_of_crystallinity(80.0, 20.0);
        assert!((xc - 0.8).abs() < 0.001); // 80% crystalline
    }
}
