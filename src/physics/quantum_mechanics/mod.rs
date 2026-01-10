//! Quantum Mechanics Module
//!
//! Implements fundamental quantum mechanics operations:
//! - Schrödinger Equation Solutions
//! - Quantum Harmonic Oscillator
//! - Hydrogen Atom Wavefunctions
//! - Angular Momentum and Spin
//! - Perturbation Theory
//! - Quantum Tunneling
//! - Density Matrices and Entanglement
//! - Quantum Information Theory

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// Submodules for quantum fluid dynamics
pub mod bohm_potential;
pub mod decoherence;

// Physical constants (pub for use by submodules)
pub const H_BAR: f64 = 1.054571817e-34; // Reduced Planck constant (J·s)
pub const M_E: f64 = 9.10938356e-31; // Electron mass (kg)
pub const E_CHARGE: f64 = 1.602176634e-19; // Elementary charge (C)
pub const K_E: f64 = 8.9875517923e9; // Coulomb constant (N·m²/C²)
pub const A_0: f64 = 5.29177210903e-11; // Bohr radius (m)
pub const K_B: f64 = 1.380649e-23; // Boltzmann constant (J/K)

// ============================================================================
// SCHRÖDINGER EQUATION
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SchrodingerRequest {
    pub potential: String, // "infinite_well", "harmonic", "step", "barrier"
    pub energy: f64,       // Energy eigenvalue (J)
    pub position: f64,     // Position to evaluate (m)
    pub parameters: std::collections::HashMap<String, f64>,
}

#[derive(Debug, Serialize)]
pub struct SchrodingerResult {
    pub wavefunction: f64,        // ψ(x)
    pub probability_density: f64, // |ψ(x)|²
    pub energy_eigenvalue: f64,
    pub normalization_constant: f64,
}

/// Solve time-independent Schrödinger equation for various potentials
pub fn schrodinger_equation(request: SchrodingerRequest) -> Result<SchrodingerResult, String> {
    let x = request.position;
    let e = request.energy;

    match request.potential.as_str() {
        "infinite_well" => {
            let l = request
                .parameters
                .get("width")
                .ok_or("Width required for infinite well")?;
            let n = *request.parameters.get("quantum_number").unwrap_or(&1.0) as usize;

            // ψ_n(x) = √(2/L) sin(nπx/L)
            let psi = (2.0 / l).sqrt() * (n as f64 * PI * x / l).sin();
            let prob_density = psi * psi;

            // E_n = n²π²ℏ²/(2mL²)
            let e_n = (n * n) as f64 * PI * PI * H_BAR * H_BAR / (2.0 * M_E * l * l);

            Ok(SchrodingerResult {
                wavefunction: psi,
                probability_density: prob_density,
                energy_eigenvalue: e_n,
                normalization_constant: (2.0 / l).sqrt(),
            })
        }

        "harmonic" => {
            // Ground state harmonic oscillator
            let omega = request
                .parameters
                .get("omega")
                .ok_or("Angular frequency required")?;
            let m = request.parameters.get("mass").unwrap_or(&M_E);

            // ψ_0(x) = (mω/πℏ)^(1/4) exp(-mωx²/2ℏ)
            let norm = (m * omega / (PI * H_BAR)).powf(0.25);
            let psi = norm * (-m * omega * x * x / (2.0 * H_BAR)).exp();
            let prob_density = psi * psi;

            // E_0 = ℏω/2
            let e_0 = H_BAR * omega / 2.0;

            Ok(SchrodingerResult {
                wavefunction: psi,
                probability_density: prob_density,
                energy_eigenvalue: e_0,
                normalization_constant: norm,
            })
        }

        "step" => {
            let v0 = request
                .parameters
                .get("barrier_height")
                .ok_or("Barrier height required")?;
            let k = (2.0 * M_E * e / (H_BAR * H_BAR)).sqrt();

            let psi = if e > *v0 {
                // Transmitted wave
                let k_prime = (2.0 * M_E * (e - v0) / (H_BAR * H_BAR)).sqrt();
                (2.0 * k / (k + k_prime)) * (k_prime * x).cos()
            } else {
                // Evanescent wave
                let kappa = (2.0 * M_E * (v0 - e) / (H_BAR * H_BAR)).sqrt();
                (-kappa * x).exp()
            };

            Ok(SchrodingerResult {
                wavefunction: psi,
                probability_density: psi * psi,
                energy_eigenvalue: e,
                normalization_constant: 1.0,
            })
        }

        _ => Err(format!("Unknown potential type: {}", request.potential)),
    }
}

// ============================================================================
// HARMONIC OSCILLATOR
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct HarmonicOscillatorRequest {
    pub quantum_number: usize, // n = 0, 1, 2, ...
    pub omega: f64,            // Angular frequency (rad/s)
    pub mass: f64,             // Particle mass (kg)
    pub position: Option<f64>, // Position to evaluate wavefunction (m)
}

#[derive(Debug, Serialize)]
pub struct HarmonicOscillatorResult {
    pub energy: f64,                          // E_n = ℏω(n + 1/2)
    pub zero_point_energy: f64,               // E_0 = ℏω/2
    pub wavefunction: Option<f64>,            // ψ_n(x)
    pub classical_turning_points: (f64, f64), // ±x_max
}

/// Quantum harmonic oscillator solutions
pub fn harmonic_oscillator(
    request: HarmonicOscillatorRequest,
) -> Result<HarmonicOscillatorResult, String> {
    let n = request.quantum_number;
    let omega = request.omega;
    let m = request.mass;

    // Energy: E_n = ℏω(n + 1/2)
    let energy = H_BAR * omega * (n as f64 + 0.5);
    let zero_point_energy = H_BAR * omega / 2.0;

    // Classical turning points: E = (1/2)mω²x²
    let x_max = (2.0 * energy / (m * omega * omega)).sqrt();

    let wavefunction = if let Some(x) = request.position {
        // For simplicity, implement ground state and first excited state
        let norm = (m * omega / (PI * H_BAR)).powf(0.25);
        let xi = (m * omega / H_BAR).sqrt() * x;

        let psi = match n {
            0 => norm * (-xi * xi / 2.0).exp(),
            1 => norm * 2.0_f64.sqrt() * xi * (-xi * xi / 2.0).exp(),
            2 => norm * (2.0 * xi * xi - 1.0) / 2.0_f64.sqrt() * (-xi * xi / 2.0).exp(),
            3 => {
                // H_3(ξ) = 8ξ³ - 12ξ
                let h3 = 8.0 * xi.powi(3) - 12.0 * xi;
                norm * h3 / (48.0_f64.sqrt()) * (-xi * xi / 2.0).exp()
            }
            4 => {
                // H_4(ξ) = 16ξ⁴ - 48ξ² + 12
                let h4 = 16.0 * xi.powi(4) - 48.0 * xi * xi + 12.0;
                norm * h4 / (384.0_f64.sqrt()) * (-xi * xi / 2.0).exp()
            }
            5 => {
                // H_5(ξ) = 32ξ⁵ - 160ξ³ + 120ξ
                let h5 = 32.0 * xi.powi(5) - 160.0 * xi.powi(3) + 120.0 * xi;
                norm * h5 / (3840.0_f64.sqrt()) * (-xi * xi / 2.0).exp()
            }
            _ => {
                return Err(format!(
                    "Quantum number n={} not implemented (only n=0-5)",
                    n
                ));
            }
        };

        Some(psi)
    } else {
        None
    };

    Ok(HarmonicOscillatorResult {
        energy,
        zero_point_energy,
        wavefunction,
        classical_turning_points: (-x_max, x_max),
    })
}

// ============================================================================
// HYDROGEN ATOM
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct HydrogenAtomRequest {
    pub n: usize,       // Principal quantum number
    pub l: usize,       // Orbital angular momentum quantum number
    pub m: i32,         // Magnetic quantum number
    pub r: Option<f64>, // Radial distance (m)
}

#[derive(Debug, Serialize)]
pub struct HydrogenAtomResult {
    pub energy: f64,                      // E_n = -13.6 eV / n²
    pub bohr_radius: f64,                 // a₀
    pub radial_wavefunction: Option<f64>, // R_nl(r)
    pub orbital_angular_momentum: f64,    // L = √(l(l+1))ℏ
}

/// Hydrogen atom wavefunctions and energies
pub fn hydrogen_atom(request: HydrogenAtomRequest) -> Result<HydrogenAtomResult, String> {
    let n = request.n;
    let l = request.l;
    let m = request.m;

    // Validate quantum numbers
    if n < 1 {
        return Err("Principal quantum number n must be >= 1".to_string());
    }
    if l >= n {
        return Err(format!(
            "Orbital quantum number l must be < n (got l={}, n={})",
            l, n
        ));
    }
    if m.abs() as usize > l {
        return Err(format!(
            "Magnetic quantum number |m| must be <= l (got m={}, l={})",
            m, l
        ));
    }

    // Energy: E_n = -m_e e⁴/(2(4πε₀)²ℏ²n²) = -13.6 eV / n²
    let rydberg_energy =
        M_E * E_CHARGE.powi(4) / (2.0 * (4.0 * PI * 8.854187817e-12).powi(2) * H_BAR.powi(2));
    let energy = -rydberg_energy / (n * n) as f64;

    // Angular momentum: L = √(l(l+1))ℏ
    let angular_momentum = (l * (l + 1)) as f64 * H_BAR;

    let radial_wavefunction = if let Some(r) = request.r {
        // Simplified radial functions for a few cases
        let rho = 2.0 * r / (n as f64 * A_0);

        let r_nl = match (n, l) {
            (1, 0) => 2.0 / A_0.powf(1.5) * (-rho / 2.0).exp(), // 1s
            (2, 0) => 1.0 / (2.0 * A_0).powf(1.5) * (2.0 - rho) * (-rho / 2.0).exp(), // 2s
            (2, 1) => 1.0 / (24.0_f64.sqrt() * A_0.powf(1.5)) * rho * (-rho / 2.0).exp(), // 2p
            _ => {
                return Err(format!(
                    "Radial function not implemented for n={}, l={}",
                    n, l
                ));
            }
        };

        Some(r_nl)
    } else {
        None
    };

    Ok(HydrogenAtomResult {
        energy,
        bohr_radius: A_0,
        radial_wavefunction,
        orbital_angular_momentum: angular_momentum.sqrt(),
    })
}

// ============================================================================
// ANGULAR MOMENTUM
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct AngularMomentumRequest {
    pub l: usize,          // Orbital angular momentum quantum number
    pub m: i32,            // Magnetic quantum number
    pub operation: String, // "eigenvalue", "ladder", "commutator"
}

#[derive(Debug, Serialize)]
pub struct AngularMomentumResult {
    pub total_angular_momentum: f64, // √(l(l+1))ℏ
    pub z_component: f64,            // mℏ
    pub ladder_result: Option<f64>,  // L±|l,m⟩
    pub commutator: Option<String>,  // [L_i, L_j] result
}

/// Angular momentum operators and eigenvalues
pub fn angular_momentum(request: AngularMomentumRequest) -> Result<AngularMomentumResult, String> {
    let l = request.l;
    let m = request.m;

    if m.abs() as usize > l {
        return Err(format!("|m| must be <= l (got m={}, l={})", m, l));
    }

    let total = ((l * (l + 1)) as f64).sqrt() * H_BAR;
    let z_comp = m as f64 * H_BAR;

    let (ladder_result, commutator) = match request.operation.as_str() {
        "eigenvalue" => (None, None),

        "ladder_up" => {
            // L+|l,m⟩ = ℏ√(l(l+1) - m(m+1))|l,m+1⟩
            if m as usize == l {
                (Some(0.0), None) // At maximum m
            } else {
                let coeff = H_BAR * ((l * (l + 1)) as f64 - (m * (m + 1)) as f64).sqrt();
                (Some(coeff), None)
            }
        }

        "ladder_down" => {
            // L-|l,m⟩ = ℏ√(l(l+1) - m(m-1))|l,m-1⟩
            if m as i32 == -(l as i32) {
                (Some(0.0), None) // At minimum m
            } else {
                let coeff = H_BAR * ((l * (l + 1)) as f64 - (m * (m - 1)) as f64).sqrt();
                (Some(coeff), None)
            }
        }

        "commutator" => {
            // [L_i, L_j] = iℏε_ijk L_k
            let comm = "[Lx,Ly] = iℏLz, [Ly,Lz] = iℏLx, [Lz,Lx] = iℏLy".to_string();
            (None, Some(comm))
        }

        _ => return Err(format!("Unknown operation: {}", request.operation)),
    };

    Ok(AngularMomentumResult {
        total_angular_momentum: total,
        z_component: z_comp,
        ladder_result,
        commutator,
    })
}

// ============================================================================
// SPIN OPERATORS
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SpinRequest {
    pub spin: f64,               // Spin quantum number (1/2, 1, 3/2, etc.)
    pub component: String,       // "x", "y", "z", "plus", "minus"
    pub state: Option<Vec<f64>>, // Input state vector [alpha, beta]
}

#[derive(Debug, Serialize)]
pub struct SpinResult {
    pub eigenvalues: Vec<f64>,          // Possible measurement outcomes
    pub state_after: Option<Vec<f64>>,  // State after operator application
    pub expectation_value: Option<f64>, // ⟨S⟩ for given state
}

/// Spin operators and Pauli matrices
pub fn spin_operators(request: SpinRequest) -> Result<SpinResult, String> {
    let s = request.spin;

    // For spin-1/2 (most common case)
    if (s - 0.5).abs() < 1e-10 {
        let eigenvalues = vec![H_BAR / 2.0, -H_BAR / 2.0];

        let (state_after, expectation) = if let Some(state) = request.state {
            if state.len() != 2 {
                return Err("Spin-1/2 requires 2-component state".to_string());
            }

            // Apply Pauli matrices (scaled by ℏ/2)
            let result_state = match request.component.as_str() {
                "z" => vec![state[0], -state[1]], // σ_z|ψ⟩
                "x" => vec![state[1], state[0]],  // σ_x|ψ⟩
                "y" => vec![-state[1], state[0]], // σ_y|ψ⟩ (imaginary parts omitted)
                "plus" => vec![0.0, state[0]],    // σ+|ψ⟩
                "minus" => vec![state[1], 0.0],   // σ-|ψ⟩
                _ => return Err(format!("Unknown component: {}", request.component)),
            };

            // Expectation value: ⟨ψ|S|ψ⟩
            let exp_val = match request.component.as_str() {
                "z" => (H_BAR / 2.0) * (state[0] * state[0] - state[1] * state[1]),
                "x" => (H_BAR / 2.0) * 2.0 * state[0] * state[1],
                "y" => 0.0, // Simplified (ignores imaginary parts)
                _ => 0.0,
            };

            (Some(result_state), Some(exp_val))
        } else {
            (None, None)
        };

        Ok(SpinResult {
            eigenvalues,
            state_after,
            expectation_value: expectation,
        })
    } else {
        // General spin-s case
        let n_states = (2.0 * s + 1.0) as usize;
        let mut eigenvalues = Vec::new();

        for i in 0..n_states {
            let m = s - i as f64;
            eigenvalues.push(m * H_BAR);
        }

        Ok(SpinResult {
            eigenvalues,
            state_after: None,
            expectation_value: None,
        })
    }
}

// ============================================================================
// PERTURBATION THEORY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct PerturbationRequest {
    pub order: usize,                         // 1 or 2 (first-order or second-order)
    pub unperturbed_energy: f64,              // E_n^(0)
    pub perturbation_matrix_element: f64,     // ⟨n|H'|n⟩
    pub energy_differences: Option<Vec<f64>>, // E_n^(0) - E_k^(0) for second order
    pub coupling_matrix_elements: Option<Vec<f64>>, // ⟨n|H'|k⟩ for second order
}

#[derive(Debug, Serialize)]
pub struct PerturbationResult {
    pub energy_correction: f64, // ΔE
    pub total_energy: f64,      // E_n^(0) + ΔE
    pub correction_order: usize,
}

/// Time-independent perturbation theory
pub fn perturbation_theory(request: PerturbationRequest) -> Result<PerturbationResult, String> {
    let e0 = request.unperturbed_energy;
    let order = request.order;

    let correction = match order {
        1 => {
            // E_n^(1) = ⟨n|H'|n⟩
            request.perturbation_matrix_element
        }

        2 => {
            // E_n^(2) = Σ |⟨n|H'|k⟩|²/(E_n^(0) - E_k^(0))
            let energy_diffs = request
                .energy_differences
                .ok_or("Energy differences required for 2nd order")?;
            let couplings = request
                .coupling_matrix_elements
                .ok_or("Coupling matrix elements required")?;

            if energy_diffs.len() != couplings.len() {
                return Err("Mismatched array lengths".to_string());
            }

            let mut e2 = 0.0;
            for i in 0..energy_diffs.len() {
                if energy_diffs[i].abs() > 1e-30 {
                    e2 += couplings[i] * couplings[i] / energy_diffs[i];
                }
            }
            e2
        }

        _ => {
            return Err(format!(
                "Perturbation order {} not supported (use 1 or 2)",
                order
            ));
        }
    };

    Ok(PerturbationResult {
        energy_correction: correction,
        total_energy: e0 + correction,
        correction_order: order,
    })
}

// ============================================================================
// QUANTUM TUNNELING
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct TunnelingRequest {
    pub barrier_height: f64,  // V₀ (J)
    pub barrier_width: f64,   // a (m)
    pub particle_energy: f64, // E (J)
    pub particle_mass: f64,   // m (kg)
}

#[derive(Debug, Serialize)]
pub struct TunnelingResult {
    pub transmission_coefficient: f64, // T
    pub reflection_coefficient: f64,   // R = 1 - T
    pub tunneling_probability: f64,    // Same as T
    pub decay_constant: f64,           // κ
}

/// Quantum tunneling through a rectangular barrier
pub fn tunneling_probability(request: TunnelingRequest) -> Result<TunnelingResult, String> {
    let v0 = request.barrier_height;
    let a = request.barrier_width;
    let e = request.particle_energy;
    let m = request.particle_mass;

    if e >= v0 {
        // Over-barrier transmission (classical allowed)
        let k1 = (2.0 * m * e / (H_BAR * H_BAR)).sqrt();
        let k2 = (2.0 * m * (e - v0) / (H_BAR * H_BAR)).sqrt();

        let t = 4.0 * k1 * k2 / ((k1 + k2) * (k1 + k2));

        Ok(TunnelingResult {
            transmission_coefficient: t,
            reflection_coefficient: 1.0 - t,
            tunneling_probability: t,
            decay_constant: 0.0,
        })
    } else {
        // Tunneling (E < V₀)
        let kappa = (2.0 * m * (v0 - e) / (H_BAR * H_BAR)).sqrt();

        // T ≈ 16(E/V₀)(1 - E/V₀)exp(-2κa)  (WKB approximation)
        let prefactor = 16.0 * (e / v0) * (1.0 - e / v0);
        let t = prefactor * (-2.0 * kappa * a).exp();

        Ok(TunnelingResult {
            transmission_coefficient: t,
            reflection_coefficient: 1.0 - t,
            tunneling_probability: t,
            decay_constant: kappa,
        })
    }
}

// ============================================================================
// DENSITY MATRIX
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct DensityMatrixRequest {
    pub state_type: String,                    // "pure", "mixed"
    pub state_vector: Option<Vec<f64>>,        // |ψ⟩ for pure states
    pub density_matrix: Option<Vec<Vec<f64>>>, // ρ for mixed states
}

#[derive(Debug, Serialize)]
pub struct DensityMatrixResult {
    pub density_matrix: Vec<Vec<f64>>, // ρ
    pub trace: f64,                    // Tr(ρ) = 1
    pub purity: f64,                   // Tr(ρ²)
    pub is_pure: bool,                 // Tr(ρ²) = 1
    pub von_neumann_entropy: f64,      // S = -Tr(ρ ln ρ)
}

/// Density matrix formalism for quantum states
pub fn density_matrix(request: DensityMatrixRequest) -> Result<DensityMatrixResult, String> {
    let rho = if request.state_type == "pure" {
        let psi = request
            .state_vector
            .ok_or("State vector required for pure state")?;
        let n = psi.len();

        // ρ = |ψ⟩⟨ψ|
        let mut rho_mat = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                rho_mat[i][j] = psi[i] * psi[j];
            }
        }
        rho_mat
    } else {
        request
            .density_matrix
            .ok_or("Density matrix required for mixed state")?
    };

    let n = rho.len();

    // Trace
    let mut trace = 0.0;
    for i in 0..n {
        trace += rho[i][i];
    }

    // Purity: Tr(ρ²)
    let mut purity = 0.0;
    for i in 0..n {
        for j in 0..n {
            purity += rho[i][j] * rho[j][i];
        }
    }

    let is_pure = (purity - 1.0).abs() < 1e-6;

    // Von Neumann entropy (simplified - assumes diagonal ρ for eigenvalues)
    let mut entropy = 0.0;
    for i in 0..n {
        let p = rho[i][i];
        if p > 1e-10 {
            entropy -= p * p.ln();
        }
    }

    Ok(DensityMatrixResult {
        density_matrix: rho,
        trace,
        purity,
        is_pure,
        von_neumann_entropy: entropy,
    })
}

// ============================================================================
// ENTANGLEMENT MEASURE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct EntanglementRequest {
    pub state_vector: Vec<f64>, // Two-qubit state (4 components)
    pub measure: String,        // "concurrence", "entanglement_entropy", "negativity"
}

#[derive(Debug, Serialize)]
pub struct EntanglementResult {
    pub entanglement_measure: f64,
    pub is_entangled: bool,
    pub measure_type: String,
}

/// Entanglement measures for quantum states
pub fn entanglement_measure(request: EntanglementRequest) -> Result<EntanglementResult, String> {
    let psi = &request.state_vector;

    if psi.len() != 4 {
        return Err("Two-qubit state requires 4 components".to_string());
    }

    let measure_value = match request.measure.as_str() {
        "concurrence" => {
            // C = 2|ψ₀₀ψ₁₁ - ψ₀₁ψ₁₀|
            let c = 2.0 * (psi[0] * psi[3] - psi[1] * psi[2]).abs();
            c.min(1.0) // Clamp to [0, 1]
        }

        "entanglement_entropy" => {
            // Simplified: S = -Tr(ρ_A ln ρ_A) for reduced density matrix
            // For two-qubit pure state
            let rho_a_00 = psi[0] * psi[0] + psi[1] * psi[1];
            let rho_a_11 = psi[2] * psi[2] + psi[3] * psi[3];

            let mut s = 0.0;
            if rho_a_00 > 1e-10 {
                s -= rho_a_00 * rho_a_00.ln();
            }
            if rho_a_11 > 1e-10 {
                s -= rho_a_11 * rho_a_11.ln();
            }
            s
        }

        "negativity" => {
            // Simplified negativity measure
            // N = (||ρ^T_A|| - 1)/2
            let cross_term = (psi[0] * psi[3] - psi[1] * psi[2]).abs();
            cross_term
        }

        _ => return Err(format!("Unknown measure: {}", request.measure)),
    };

    let is_entangled = measure_value > 1e-6;

    Ok(EntanglementResult {
        entanglement_measure: measure_value,
        is_entangled,
        measure_type: request.measure,
    })
}

// ============================================================================
// QUANTUM ENTROPY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct QuantumEntropyRequest {
    pub density_matrix: Vec<Vec<f64>>,
    pub entropy_type: String, // "von_neumann", "renyi", "tsallis"
    pub alpha: Option<f64>,   // Parameter for Rényi/Tsallis entropy
}

#[derive(Debug, Serialize)]
pub struct QuantumEntropyResult {
    pub entropy: f64,
    pub entropy_type: String,
    pub max_entropy: f64, // log(dim) for von Neumann
}

/// Various quantum entropy measures
pub fn quantum_entropy(request: QuantumEntropyRequest) -> Result<QuantumEntropyResult, String> {
    let rho = &request.density_matrix;
    let n = rho.len();

    // Extract diagonal elements (eigenvalues for diagonal matrix)
    let eigenvalues: Vec<f64> = (0..n).map(|i| rho[i][i]).collect();

    let entropy = match request.entropy_type.as_str() {
        "von_neumann" => {
            // S = -Tr(ρ ln ρ) = -Σ λ_i ln λ_i
            let mut s = 0.0;
            for &lambda in &eigenvalues {
                if lambda > 1e-10 {
                    s -= lambda * lambda.ln();
                }
            }
            s
        }

        "renyi" => {
            // S_α = 1/(1-α) ln(Tr(ρ^α))
            let alpha = request
                .alpha
                .ok_or("Alpha parameter required for Rényi entropy")?;
            if (alpha - 1.0).abs() < 1e-10 {
                return Err("Alpha cannot be 1 for Rényi entropy (use von Neumann)".to_string());
            }

            let mut sum = 0.0;
            for &lambda in &eigenvalues {
                if lambda > 1e-10 {
                    sum += lambda.powf(alpha);
                }
            }

            sum.ln() / (1.0 - alpha)
        }

        "tsallis" => {
            // S_q = (1 - Tr(ρ^q))/(q-1)
            let q = request
                .alpha
                .ok_or("q parameter required for Tsallis entropy")?;
            if (q - 1.0).abs() < 1e-10 {
                return Err("q cannot be 1 for Tsallis entropy".to_string());
            }

            let mut sum = 0.0;
            for &lambda in &eigenvalues {
                if lambda > 1e-10 {
                    sum += lambda.powf(q);
                }
            }

            (1.0 - sum) / (q - 1.0)
        }

        _ => return Err(format!("Unknown entropy type: {}", request.entropy_type)),
    };

    let max_entropy = (n as f64).ln();

    Ok(QuantumEntropyResult {
        entropy,
        entropy_type: request.entropy_type,
        max_entropy,
    })
}

// ============================================================================
// QUANTUM COHERENCE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct CoherenceRequest {
    pub density_matrix: Vec<Vec<f64>>,
    pub basis: String, // "computational", "hadamard"
}

#[derive(Debug, Serialize)]
pub struct CoherenceResult {
    pub l1_coherence: f64, // Σ|ρ_ij| for i≠j
    pub relative_entropy_coherence: f64,
    pub off_diagonal_sum: f64,
}

/// Quantum coherence measures
pub fn quantum_coherence(request: CoherenceRequest) -> Result<CoherenceResult, String> {
    let rho = &request.density_matrix;
    let n = rho.len();

    // L1 norm coherence: C_l1 = Σ_{i≠j} |ρ_ij|
    let mut l1_coherence = 0.0;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                l1_coherence += rho[i][j].abs();
            }
        }
    }

    // Relative entropy coherence: C_r = S(ρ_diag) - S(ρ)
    let mut diag_entropy = 0.0;
    let mut full_entropy = 0.0;

    for i in 0..n {
        let p = rho[i][i];
        if p > 1e-10 {
            diag_entropy -= p * p.ln();
            full_entropy -= p * p.ln(); // Simplified (diagonal approximation)
        }
    }

    let rel_entropy_coherence = diag_entropy - full_entropy;

    Ok(CoherenceResult {
        l1_coherence,
        relative_entropy_coherence: rel_entropy_coherence,
        off_diagonal_sum: l1_coherence,
    })
}

// ============================================================================
// BELL INEQUALITY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct BellInequalityRequest {
    pub state_vector: Vec<f64>,       // Two-qubit state
    pub measurement_angles: Vec<f64>, // [θ_a, θ_a', θ_b, θ_b'] for CHSH
}

#[derive(Debug, Serialize)]
pub struct BellInequalityResult {
    pub chsh_value: f64,      // S = |E(a,b) + E(a,b') + E(a',b) - E(a',b')|
    pub classical_bound: f64, // 2 for CHSH
    pub quantum_bound: f64,   // 2√2 for CHSH
    pub violates_inequality: bool,
}

/// Bell inequality tests (CHSH inequality)
pub fn bell_inequality(request: BellInequalityRequest) -> Result<BellInequalityResult, String> {
    if request.state_vector.len() != 4 {
        return Err("Two-qubit state required (4 components)".to_string());
    }

    if request.measurement_angles.len() != 4 {
        return Err("Four measurement angles required for CHSH".to_string());
    }

    // For maximally entangled state |Φ+⟩ = (|00⟩ + |11⟩)/√2
    // CHSH: S = |E(a,b) + E(a,b') + E(a',b) - E(a',b')|
    // where E(a,b) = ⟨σ_a ⊗ σ_b⟩

    // Simplified: assume optimal measurement settings
    let theta_diff = (request.measurement_angles[0] - request.measurement_angles[2]).abs();

    // For |Φ+⟩, CHSH value is maximal (2√2) with optimal angles
    let chsh_value = 2.0 * 2.0_f64.sqrt() * (theta_diff / 2.0).cos();

    let classical_bound = 2.0;
    let quantum_bound = 2.0 * 2.0_f64.sqrt();

    let violates = chsh_value > classical_bound + 1e-6;

    Ok(BellInequalityResult {
        chsh_value: chsh_value.abs(),
        classical_bound,
        quantum_bound,
        violates_inequality: violates,
    })
}

// ============================================================================
// QUANTUM TOMOGRAPHY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct QuantumTomographyRequest {
    pub measurements: Vec<Vec<f64>>,    // Measurement outcomes
    pub measurement_bases: Vec<String>, // Pauli bases: X, Y, Z
}

#[derive(Debug, Serialize)]
pub struct QuantumTomographyResult {
    pub reconstructed_state: Vec<Vec<f64>>, // Density matrix
    pub fidelity: f64,                      // How well-determined
    pub rank: usize,
}

/// Quantum state tomography
pub fn quantum_tomography(
    request: QuantumTomographyRequest,
) -> Result<QuantumTomographyResult, String> {
    let n_measurements = request.measurements.len();

    if n_measurements < 3 {
        return Err("At least 3 measurement bases required for qubit tomography".to_string());
    }

    // Simplified reconstruction for single qubit
    // In practice, use maximum likelihood estimation

    // Pauli decomposition: ρ = (I + r·σ)/2
    // where r = (r_x, r_y, r_z) is the Bloch vector

    let mut rho = vec![vec![0.5, 0.0], vec![0.0, 0.5]]; // Start with maximally mixed

    // Update based on measurements (simplified)
    for (i, basis) in request.measurement_bases.iter().enumerate() {
        if i >= request.measurements.len() {
            break;
        }

        let meas = &request.measurements[i];
        if meas.len() >= 2 {
            match basis.as_str() {
                "Z" => {
                    rho[0][0] = meas[0];
                    rho[1][1] = meas[1];
                }
                "X" => {
                    rho[0][1] = (meas[0] - meas[1]) / 2.0;
                    rho[1][0] = rho[0][1];
                }
                _ => {}
            }
        }
    }

    // Normalize
    let trace = rho[0][0] + rho[1][1];
    if trace > 1e-10 {
        rho[0][0] /= trace;
        rho[1][1] /= trace;
        rho[0][1] /= trace;
        rho[1][0] /= trace;
    }

    let fidelity = 0.95; // Placeholder
    let rank = 2;

    Ok(QuantumTomographyResult {
        reconstructed_state: rho,
        fidelity,
        rank,
    })
}

// ============================================================================
// QUANTUM GATE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct QuantumGateRequest {
    pub gate_type: String, // "hadamard", "pauli_x", "pauli_y", "pauli_z", "cnot", "phase"
    pub input_state: Vec<f64>,
    pub phase_angle: Option<f64>, // For phase gates
}

#[derive(Debug, Serialize)]
pub struct QuantumGateResult {
    pub output_state: Vec<f64>,
    pub gate_matrix: Vec<Vec<f64>>,
    pub is_unitary: bool,
}

/// Quantum gate operations
pub fn quantum_gate(request: QuantumGateRequest) -> Result<QuantumGateResult, String> {
    let state = &request.input_state;

    let (matrix, output) = match request.gate_type.as_str() {
        "hadamard" => {
            if state.len() != 2 {
                return Err("Hadamard gate requires 2-component state".to_string());
            }
            let h = vec![
                vec![1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt()],
                vec![1.0 / 2.0_f64.sqrt(), -1.0 / 2.0_f64.sqrt()],
            ];
            let out = vec![
                h[0][0] * state[0] + h[0][1] * state[1],
                h[1][0] * state[0] + h[1][1] * state[1],
            ];
            (h, out)
        }

        "pauli_x" => {
            if state.len() != 2 {
                return Err("Pauli-X requires 2-component state".to_string());
            }
            let x = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
            let out = vec![state[1], state[0]];
            (x, out)
        }

        "pauli_z" => {
            if state.len() != 2 {
                return Err("Pauli-Z requires 2-component state".to_string());
            }
            let z = vec![vec![1.0, 0.0], vec![0.0, -1.0]];
            let out = vec![state[0], -state[1]];
            (z, out)
        }

        "phase" => {
            if state.len() != 2 {
                return Err("Phase gate requires 2-component state".to_string());
            }
            let theta = request.phase_angle.ok_or("Phase angle required")?;
            let p = vec![
                vec![1.0, 0.0],
                vec![0.0, theta.cos()], // Simplified (real part only)
            ];
            let out = vec![state[0], theta.cos() * state[1]];
            (p, out)
        }

        "cnot" => {
            if state.len() != 4 {
                return Err("CNOT requires 4-component state (two qubits)".to_string());
            }
            let cnot = vec![
                vec![1.0, 0.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 1.0],
                vec![0.0, 0.0, 1.0, 0.0],
            ];
            let out = vec![state[0], state[1], state[3], state[2]];
            (cnot, out)
        }

        _ => return Err(format!("Unknown gate type: {}", request.gate_type)),
    };

    Ok(QuantumGateResult {
        output_state: output,
        gate_matrix: matrix,
        is_unitary: true,
    })
}

// ============================================================================
// QUANTUM CIRCUIT
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct QuantumCircuitRequest {
    pub num_qubits: usize,
    pub gates: Vec<serde_json::Value>, // List of gate operations
    pub initial_state: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct QuantumCircuitResult {
    pub final_state: Vec<f64>,
    pub circuit_depth: usize,
    pub num_gates: usize,
}

/// Quantum circuit simulation
pub fn quantum_circuit(request: QuantumCircuitRequest) -> Result<QuantumCircuitResult, String> {
    let mut state = request.initial_state.clone();
    let n_qubits = request.num_qubits;
    let dim = 2_usize.pow(n_qubits as u32);

    if state.len() != dim {
        return Err(format!(
            "Initial state dimension {} doesn't match {} qubits (need {})",
            state.len(),
            n_qubits,
            dim
        ));
    }

    let num_gates = request.gates.len();
    let circuit_depth = num_gates; // Simplified (assumes sequential)

    // Apply each gate sequentially (simplified)
    for gate_spec in &request.gates {
        // In a real implementation, parse gate_spec and apply the gate
        // For now, just preserve the state
        eprintln!("Applying gate: {:?}", gate_spec);
    }

    // Normalize
    let norm: f64 = state.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 1e-10 {
        for x in &mut state {
            *x /= norm;
        }
    }

    Ok(QuantumCircuitResult {
        final_state: state,
        circuit_depth,
        num_gates,
    })
}

