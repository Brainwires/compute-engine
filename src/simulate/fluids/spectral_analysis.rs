//! Spectral Analysis for Fluid Dynamics
//!
//! Provides tools for analyzing turbulent flows in spectral space:
//! - Kinetic energy spectrum E(k) extraction
//! - Enstrophy spectrum
//! - Dissipation rate calculation
//! - Kolmogorov scale estimation
//!
//! Uses FFT from the signal processing module.

use ndarray::{Array2, Array3};
use rustfft::{FftPlanner, num_complex::Complex};
use std::f64::consts::PI;

/// Result of energy spectrum calculation
#[derive(Debug, Clone)]
pub struct EnergySpectrumResult {
    /// Wavenumber bins k
    pub wavenumbers: Vec<f64>,

    /// Energy per wavenumber E(k)
    pub energy: Vec<f64>,

    /// Total kinetic energy
    pub total_energy: f64,

    /// Grid spacing used
    pub dx: f64,
}

/// Calculate kinetic energy spectrum E(k) from 2D velocity field
///
/// E(k) = (1/2) × Σ_{k ≤ |k'| < k+1} |û(k')|²
///
/// # Arguments
/// * `u` - x-velocity field
/// * `v` - y-velocity field
/// * `dx` - Grid spacing (assumed equal in x and y)
///
/// # Returns
/// Energy spectrum binned by wavenumber magnitude
pub fn energy_spectrum_2d(u: &Array2<f64>, v: &Array2<f64>, dx: f64) -> EnergySpectrumResult {
    let (nx, ny) = u.dim();
    assert_eq!(u.dim(), v.dim(), "Velocity fields must have same dimensions");

    // Flatten and convert to complex for FFT
    let mut u_complex: Vec<Complex<f64>> = u.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut v_complex: Vec<Complex<f64>> = v.iter().map(|&x| Complex::new(x, 0.0)).collect();

    // Compute 2D FFT (as sequence of 1D FFTs)
    let mut planner = FftPlanner::new();

    // FFT rows
    let fft_row = planner.plan_fft_forward(ny);
    for i in 0..nx {
        let start = i * ny;
        let end = start + ny;
        fft_row.process(&mut u_complex[start..end]);
        fft_row.process(&mut v_complex[start..end]);
    }

    // FFT columns (need to extract, process, put back)
    let fft_col = planner.plan_fft_forward(nx);
    let mut col_u = vec![Complex::new(0.0, 0.0); nx];
    let mut col_v = vec![Complex::new(0.0, 0.0); nx];

    for j in 0..ny {
        for i in 0..nx {
            col_u[i] = u_complex[i * ny + j];
            col_v[i] = v_complex[i * ny + j];
        }
        fft_col.process(&mut col_u);
        fft_col.process(&mut col_v);
        for i in 0..nx {
            u_complex[i * ny + j] = col_u[i];
            v_complex[i * ny + j] = col_v[i];
        }
    }

    // Normalize
    let norm = 1.0 / (nx * ny) as f64;
    for c in &mut u_complex {
        *c *= norm;
    }
    for c in &mut v_complex {
        *c *= norm;
    }

    // Calculate wavenumber for each mode and bin energy
    let k_max = (nx.min(ny) / 2) as f64;
    let dk = 2.0 * PI / (dx * nx as f64);
    let n_bins = (k_max as usize).max(1);
    let mut energy_bins = vec![0.0; n_bins];
    let mut counts = vec![0usize; n_bins];

    for i in 0..nx {
        for j in 0..ny {
            // Wavenumber indices (handle negative frequencies)
            let kx_idx = if i <= nx / 2 { i as f64 } else { (i as i64 - nx as i64) as f64 };
            let ky_idx = if j <= ny / 2 { j as f64 } else { (j as i64 - ny as i64) as f64 };

            let k_mag = (kx_idx * kx_idx + ky_idx * ky_idx).sqrt();
            let bin = (k_mag as usize).min(n_bins - 1);

            let idx = i * ny + j;
            let u_hat = u_complex[idx];
            let v_hat = v_complex[idx];

            // Kinetic energy: E = 0.5 * (|û|² + |v̂|²)
            let energy = 0.5 * (u_hat.norm_sqr() + v_hat.norm_sqr());

            energy_bins[bin] += energy;
            counts[bin] += 1;
        }
    }

    // Calculate wavenumbers for output
    let wavenumbers: Vec<f64> = (0..n_bins).map(|k| k as f64 * dk).collect();
    let total_energy: f64 = energy_bins.iter().sum();

    EnergySpectrumResult {
        wavenumbers,
        energy: energy_bins,
        total_energy,
        dx,
    }
}

/// Calculate kinetic energy spectrum E(k) from 3D velocity field
pub fn energy_spectrum_3d(
    u: &Array3<f64>,
    v: &Array3<f64>,
    w: &Array3<f64>,
    dx: f64,
) -> EnergySpectrumResult {
    let (nx, ny, nz) = u.dim();
    assert_eq!(u.dim(), v.dim(), "Velocity fields must have same dimensions");
    assert_eq!(u.dim(), w.dim(), "Velocity fields must have same dimensions");

    // For 3D, we'll use a simpler approach: compute spectra in slices and average
    // Full 3D FFT would require more memory management

    let n_slices = nz;
    let mut accumulated_spectrum: Option<EnergySpectrumResult> = None;

    for k in 0..n_slices {
        // Extract 2D slice
        let u_slice: Array2<f64> = u.slice(ndarray::s![.., .., k]).to_owned();
        let v_slice: Array2<f64> = v.slice(ndarray::s![.., .., k]).to_owned();

        let slice_spectrum = energy_spectrum_2d(&u_slice, &v_slice, dx);

        match &mut accumulated_spectrum {
            None => accumulated_spectrum = Some(slice_spectrum),
            Some(acc) => {
                for (i, e) in slice_spectrum.energy.iter().enumerate() {
                    if i < acc.energy.len() {
                        acc.energy[i] += e;
                    }
                }
                acc.total_energy += slice_spectrum.total_energy;
            }
        }
    }

    // Add contribution from w component (compute 2D spectrum of w in each x-y slice)
    for k in 0..n_slices {
        let w_slice: Array2<f64> = w.slice(ndarray::s![.., .., k]).to_owned();
        let zero_slice = Array2::<f64>::zeros((nx, ny));
        let w_spectrum = energy_spectrum_2d(&w_slice, &zero_slice, dx);

        if let Some(acc) = &mut accumulated_spectrum {
            for (i, e) in w_spectrum.energy.iter().enumerate() {
                if i < acc.energy.len() {
                    acc.energy[i] += e;
                }
            }
            acc.total_energy += w_spectrum.total_energy;
        }
    }

    // Average over slices
    if let Some(ref mut result) = accumulated_spectrum {
        for e in &mut result.energy {
            *e /= n_slices as f64;
        }
        result.total_energy /= n_slices as f64;
    }

    accumulated_spectrum.unwrap_or(EnergySpectrumResult {
        wavenumbers: vec![],
        energy: vec![],
        total_energy: 0.0,
        dx,
    })
}

/// Kolmogorov dissipation scale
///
/// η = (ν³/ε)^(1/4)
///
/// # Arguments
/// * `kinematic_viscosity` - ν = μ/ρ (m²/s)
/// * `dissipation_rate` - ε (m²/s³)
pub fn kolmogorov_scale(kinematic_viscosity: f64, dissipation_rate: f64) -> f64 {
    if kinematic_viscosity <= 0.0 || dissipation_rate <= 0.0 {
        return f64::NAN;
    }
    (kinematic_viscosity.powi(3) / dissipation_rate).powf(0.25)
}

/// Estimate dissipation rate from velocity field
///
/// ε = 2ν × ∫ S_ij S_ij dV
///
/// where S_ij is the strain rate tensor.
pub fn estimate_dissipation_rate_2d(
    u: &Array2<f64>,
    v: &Array2<f64>,
    dx: f64,
    kinematic_viscosity: f64,
) -> f64 {
    let (nx, ny) = u.dim();
    let mut dissipation = 0.0;

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            // Velocity gradients
            let du_dx = (u[[i + 1, j]] - u[[i - 1, j]]) / (2.0 * dx);
            let du_dy = (u[[i, j + 1]] - u[[i, j - 1]]) / (2.0 * dx);
            let dv_dx = (v[[i + 1, j]] - v[[i - 1, j]]) / (2.0 * dx);
            let dv_dy = (v[[i, j + 1]] - v[[i, j - 1]]) / (2.0 * dx);

            // Strain rate tensor components (symmetric part of gradient)
            let s11 = du_dx;
            let s22 = dv_dy;
            let s12 = 0.5 * (du_dy + dv_dx);

            // S_ij S_ij = S11² + S22² + 2*S12²
            let s_sq = s11 * s11 + s22 * s22 + 2.0 * s12 * s12;

            dissipation += s_sq;
        }
    }

    // ε = 2ν × <S_ij S_ij>
    2.0 * kinematic_viscosity * dissipation * dx * dx / ((nx - 2) * (ny - 2)) as f64
}

/// Calculate enstrophy (integral of squared vorticity)
///
/// Ω = ∫ ω² dV
pub fn enstrophy_2d(u: &Array2<f64>, v: &Array2<f64>, dx: f64) -> f64 {
    let (nx, ny) = u.dim();
    let mut enstrophy = 0.0;

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            // Vorticity ω = ∂v/∂x - ∂u/∂y
            let dv_dx = (v[[i + 1, j]] - v[[i - 1, j]]) / (2.0 * dx);
            let du_dy = (u[[i, j + 1]] - u[[i, j - 1]]) / (2.0 * dx);
            let omega = dv_dx - du_dy;

            enstrophy += omega * omega;
        }
    }

    enstrophy * dx * dx
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_energy_spectrum_uniform() {
        // Uniform velocity field should have energy only at k=0
        let n = 64;
        let u = Array2::<f64>::from_elem((n, n), 1.0);
        let v = Array2::<f64>::from_elem((n, n), 0.0);

        let result = energy_spectrum_2d(&u, &v, 0.1);

        // Most energy should be at low k
        assert!(result.energy[0] > result.energy.iter().skip(1).sum::<f64>());
    }

    #[test]
    fn test_kolmogorov_scale() {
        // Water: ν ≈ 1e-6 m²/s, moderate turbulence: ε ≈ 0.01 m²/s³
        let nu = 1e-6;
        let epsilon = 0.01;

        let eta = kolmogorov_scale(nu, epsilon);

        // η should be in the range of tens of micrometers
        assert!(eta > 1e-5);
        assert!(eta < 1e-3);
    }

    #[test]
    fn test_enstrophy_irrotational() {
        // Uniform flow has zero vorticity
        let n = 32;
        let u = Array2::<f64>::from_elem((n, n), 1.0);
        let v = Array2::<f64>::from_elem((n, n), 0.0);

        let omega = enstrophy_2d(&u, &v, 0.1);

        assert!(omega.abs() < 1e-10);
    }

    #[test]
    fn test_enstrophy_vortex() {
        // Create a simple vortex: u = -y, v = x (solid body rotation)
        let n = 64;
        let dx = 2.0 / n as f64;
        let mut u = Array2::<f64>::zeros((n, n));
        let mut v = Array2::<f64>::zeros((n, n));

        for i in 0..n {
            for j in 0..n {
                let x = (i as f64 - n as f64 / 2.0) * dx;
                let y = (j as f64 - n as f64 / 2.0) * dx;
                u[[i, j]] = -y;
                v[[i, j]] = x;
            }
        }

        let omega = enstrophy_2d(&u, &v, dx);

        // Solid body rotation has uniform vorticity ω = 2
        // So enstrophy should be approximately 4 * area
        let expected_enstrophy = 4.0 * (2.0 * 2.0); // 4 * (2x2 domain)
        let relative_error = (omega - expected_enstrophy).abs() / expected_enstrophy;

        assert!(
            relative_error < 0.1,
            "Enstrophy error: {:.1}%",
            relative_error * 100.0
        );
    }
}
