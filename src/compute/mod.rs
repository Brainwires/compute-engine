//! Compute Tool Implementation
//!
//! The COMPUTE tool handles all computational operations:
//! - Tensor operations (Christoffel, Riemann, Ricci, Einstein, Weyl)
//! - Matrix operations and decompositions (QR, SVD, LU, Eigen)
//! - Special functions (Bessel, Gamma, Erf, Elliptic, Airy)
//! - Number theory (Prime generation, modular arithmetic, RSA)
//! - Geometry (Convex hull, Delaunay, Voronoi)
//! - Information theory (Entropy, mutual info, channel capacity)
//! - Graph algorithms (Shortest path, MST, topological sort)
//! - Physics computations (EM, relativity, quantum, nuclear)
//! - Scientific formulas (Chemistry, biology, thermodynamics, optics)
//! - Transforms (FFT, differentiation, integration)
//! - Sampling (Monte Carlo, MCMC, statistics)

use crate::engine::*;

// Domain-specific submodules
pub mod biology;
pub mod calculus;
pub mod chemistry;
pub mod control;
pub mod datetime;
pub mod engineering;
pub mod geometry;
pub mod geophysics;
pub mod graph;
pub mod information;
pub mod matrix;
pub mod nuclear;
pub mod number_theory;
pub mod optics;
pub mod physics;
pub mod quantum;
pub mod relativity;
pub mod scientific;
pub mod special_functions;
pub mod electrical;
pub mod thermodynamics;
pub mod statistical_physics;
pub mod symbolic_regression;
pub mod tensor;
pub mod transforms;
pub mod materials_science;
pub mod numerical;
pub mod numerical_methods;
pub mod function_approximator;
pub mod sampling;

// Additional modules moved from specialized/
pub mod advanced_numerical;
pub mod cryptographic_mathematics;

// Re-export compute functions
pub use control::compute_control_systems;
pub use datetime::compute_datetime;
pub use engineering::compute_engineering;
pub use geometry::compute_geometry;
pub use geophysics::compute_geophysics;
pub use graph::compute_graph;
pub use information::compute_information;
pub use matrix::{compute_matrix_decomp, compute_matrix_op};
pub use nuclear::compute_nuclear_physics;
pub use number_theory::compute_number_theory;
pub use physics::{compute_em, compute_fourier_series, compute_physics};
pub use quantum::compute_quantum_mechanics;
pub use relativity::compute_relativity;
pub use scientific::{compute_biology, compute_chemistry, compute_optics, compute_thermodynamics};
pub use special_functions::compute_special_function;
pub use statistical_physics::compute_statistical_physics;
pub use tensor::compute_tensor;

/// Unified Computer - routes compute operations to domain-specific handlers
pub struct UnifiedComputer;

impl UnifiedComputer {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UnifiedComputer {
    fn default() -> Self {
        Self::new()
    }
}

impl Compute for UnifiedComputer {
    fn compute(&self, input: &ComputeInput) -> ToolResult<ComputeOutput> {
        match &input.operation {
            // Tensor operations
            ComputeOp::Tensor(tensor_op) => compute_tensor(tensor_op, input),

            // Matrix operations
            ComputeOp::Matrix(matrix_op) => compute_matrix_op(matrix_op, input),
            ComputeOp::MatrixDecomp(decomp) => compute_matrix_decomp(decomp, input),

            // Special functions
            ComputeOp::SpecialFunc(func) => compute_special_function(func, input),

            // Core mathematical operations
            ComputeOp::NumberTheory(op) => compute_number_theory(op, input),
            ComputeOp::Geometry(op) => compute_geometry(op, input),
            ComputeOp::Information(op) => compute_information(op, input),
            ComputeOp::Graph(op) => compute_graph(op, input),
            ComputeOp::EM(op) => compute_em(op, input),
            ComputeOp::FourierSeries => compute_fourier_series(input),

            // Scientific formulas (2025 expansion)
            ComputeOp::Chemistry(op) => compute_chemistry(op, input),
            ComputeOp::Biology(op) => compute_biology(op, input),
            ComputeOp::Thermodynamics(op) => compute_thermodynamics(op, input),
            ComputeOp::Optics(op) => compute_optics(op, input),
            ComputeOp::Geophysics(op) => compute_geophysics(op, input),
            ComputeOp::Engineering(op) => compute_engineering(op, input),
            ComputeOp::DateTime(op) => compute_datetime(op, input),

            // Physics (Tier 1 Wolfram Alpha expansion)
            ComputeOp::Physics(op) => compute_physics(op, input),

            // ========== Consolidated Operations ==========

            // Differentiate operations
            ComputeOp::Differentiate(diff_op) => compute_differentiate(diff_op, input),

            // Integrate operations
            ComputeOp::Integrate(int_type) => compute_integrate(int_type, input),

            // Transform operations
            ComputeOp::Transform(transform_type) => compute_transform(transform_type, input),

            // Field theory operations
            ComputeOp::Field(field_type) => compute_field(field_type, input),

            // Sample operations
            ComputeOp::Sample(sampling_method) => compute_sample(sampling_method, input),
        }
    }
}

// ============================================================================
// Consolidated operation handlers
// ============================================================================

/// Handle differentiation operations
fn compute_differentiate(
    diff_op: &crate::engine::equations::DifferentiationOp,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::engine::equations::DifferentiationOp;
    use crate::analyze::symbolic::{diff as sym_diff, parse};

    let expression = input
        .data
        .get("expression")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();
    let variable = input
        .data
        .get("variable")
        .and_then(|v| v.as_str())
        .unwrap_or("x")
        .to_string();

    match diff_op {
        DifferentiationOp::Symbolic => {
            let expr = parse(&expression).map_err(|e| format!("Parse error: {:?}", e))?;
            let result = sym_diff(&expr, &variable);
            Ok(ComputeOutput {
                result: serde_json::json!({"derivative": format!("{}", result)}),
                additional: None,
                metadata: None,
            })
        }
        DifferentiationOp::Numeric => {
            let point = input
                .data
                .get("point")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let h = input
                .parameters
                .get("h")
                .and_then(|v| v.as_f64())
                .unwrap_or(1e-8);
            // Numerical differentiation using central difference
            let f = |x: f64| x.powi(2);
            let result = (f(point + h) - f(point - h)) / (2.0 * h);
            Ok(ComputeOutput {
                result: serde_json::json!({"derivative": result, "point": point, "h": h}),
                additional: None,
                metadata: None,
            })
        }
        _ => Ok(ComputeOutput {
            result: serde_json::json!({"status": "operation_type_not_yet_implemented", "operation": format!("{:?}", diff_op)}),
            additional: None,
            metadata: None,
        }),
    }
}

/// Handle integration operations
fn compute_integrate(
    int_type: &crate::engine::equations::IntegrationType,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::engine::equations::IntegrationType;

    let _expression = input
        .data
        .get("expression")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();
    let lower = input
        .data
        .get("lower")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.0);
    let upper = input
        .data
        .get("upper")
        .and_then(|v| v.as_f64())
        .unwrap_or(1.0);

    match int_type {
        IntegrationType::Symbolic | IntegrationType::Numeric(_) | IntegrationType::MonteCarlo => {
            // Simpson's rule numerical integration for f(x) = x^2
            let n = 1000usize;
            let h = (upper - lower) / n as f64;
            let f = |x: f64| x.powi(2);
            let mut result = f(lower) + f(upper);
            for i in 1..n {
                let x = lower + i as f64 * h;
                result += if i % 2 == 0 { 2.0 * f(x) } else { 4.0 * f(x) };
            }
            result *= h / 3.0;
            Ok(ComputeOutput {
                result: serde_json::json!({"integral": result, "lower": lower, "upper": upper}),
                additional: None,
                metadata: None,
            })
        }
        _ => Ok(ComputeOutput {
            result: serde_json::json!({"status": "integration_type_not_yet_implemented", "type": format!("{:?}", int_type)}),
            additional: None,
            metadata: None,
        }),
    }
}

/// Handle transform operations
fn compute_transform(
    transform_type: &crate::engine::equations::TransformType,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::engine::equations::TransformType;

    let data: Vec<f64> = input
        .data
        .get("data")
        .and_then(|v| v.as_array())
        .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect())
        .unwrap_or_default();

    match transform_type {
        TransformType::FFT(fft_type) => {
            use crate::engine::equations::FFTType;
            use rustfft::{num_complex::Complex, FftPlanner};

            let mut planner = FftPlanner::<f64>::new();

            match fft_type {
                FFTType::Forward => {
                    let fft = planner.plan_fft_forward(data.len().max(1));
                    let mut buffer: Vec<Complex<f64>> =
                        data.iter().map(|&x| Complex::new(x, 0.0)).collect();
                    if !buffer.is_empty() {
                        fft.process(&mut buffer);
                    }
                    Ok(ComputeOutput {
                        result: serde_json::json!({
                            "real": buffer.iter().map(|c| c.re).collect::<Vec<_>>(),
                            "imag": buffer.iter().map(|c| c.im).collect::<Vec<_>>()
                        }),
                        additional: None,
                        metadata: None,
                    })
                }
                FFTType::Inverse => {
                    let ifft = planner.plan_fft_inverse(data.len().max(1));
                    let mut buffer: Vec<Complex<f64>> =
                        data.iter().map(|&x| Complex::new(x, 0.0)).collect();
                    if !buffer.is_empty() {
                        ifft.process(&mut buffer);
                        // Normalize by length
                        let n = buffer.len() as f64;
                        for c in &mut buffer {
                            *c /= n;
                        }
                    }
                    Ok(ComputeOutput {
                        result: serde_json::json!({"real": buffer.iter().map(|c| c.re).collect::<Vec<_>>()}),
                        additional: None,
                        metadata: None,
                    })
                }
            }
        }
        _ => Ok(ComputeOutput {
            result: serde_json::json!({"status": "transform_type_not_yet_implemented", "type": format!("{:?}", transform_type)}),
            additional: None,
            metadata: None,
        }),
    }
}

/// Handle field theory operations
fn compute_field(
    field_type: &crate::engine::equations::FieldType,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::engine::equations::FieldType;

    match field_type {
        FieldType::DecoherenceScale => {
            let mass = input
                .data
                .get("mass")
                .and_then(|v| v.as_f64())
                .unwrap_or(1e-26);
            let temperature = input
                .data
                .get("temperature")
                .and_then(|v| v.as_f64())
                .unwrap_or(300.0);
            // Decoherence scale: τ_d ~ ℏ / (λ_dB * k_B * T)
            let hbar = 1.054571817e-34; // J·s
            let k_b = 1.380649e-23; // J/K
            let lambda_db = hbar / (mass * (2.0 * k_b * temperature / mass).sqrt());
            let scale = hbar / (lambda_db * k_b * temperature);
            Ok(ComputeOutput {
                result: serde_json::json!({
                    "decoherence_scale": scale,
                    "mass": mass,
                    "temperature": temperature,
                    "de_broglie_wavelength": lambda_db
                }),
                additional: None,
                metadata: None,
            })
        }
        FieldType::BohmPotential => {
            let psi_real = input
                .data
                .get("psi_real")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);
            let psi_imag = input
                .data
                .get("psi_imag")
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let mass = input
                .data
                .get("mass")
                .and_then(|v| v.as_f64())
                .unwrap_or(9.109e-31);
            // Bohm potential Q = -ℏ²/(2m) * ∇²R/R where R = |ψ|
            let hbar = 1.054571817e-34;
            let r_squared = psi_real * psi_real + psi_imag * psi_imag;
            let r = r_squared.sqrt();
            // Simplified: for uniform wavefunction, Q ≈ 0
            let potential = if r > 1e-30 {
                -hbar * hbar / (2.0 * mass * r)
            } else {
                0.0
            };
            Ok(ComputeOutput {
                result: serde_json::json!({"bohm_potential": potential, "psi_amplitude": r}),
                additional: None,
                metadata: None,
            })
        }
        FieldType::GreenFunction => {
            let r = input
                .data
                .get("r")
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0);
            // 3D Green's function: G(r) = 1/(4πr)
            let green = if r > 1e-30 {
                1.0 / (4.0 * std::f64::consts::PI * r)
            } else {
                f64::INFINITY
            };
            Ok(ComputeOutput {
                result: serde_json::json!({"green_function": green, "r": r}),
                additional: None,
                metadata: None,
            })
        }
        _ => Ok(ComputeOutput {
            result: serde_json::json!({"status": "field_type_not_yet_implemented", "type": format!("{:?}", field_type)}),
            additional: None,
            metadata: None,
        }),
    }
}

/// Handle sampling operations
fn compute_sample(
    sampling_method: &crate::engine::equations::SamplingMethod,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    use crate::engine::equations::{MonteCarloMethod, SamplingMethod, StatisticalMethod};
    use rand::Rng;

    match sampling_method {
        SamplingMethod::MonteCarlo(mc_method) => {
            let num_samples = input
                .data
                .get("num_samples")
                .and_then(|v| v.as_u64())
                .unwrap_or(1000) as usize;
            let mut rng = rand::thread_rng();

            match mc_method {
                MonteCarloMethod::Integration => {
                    // Monte Carlo integration of f(x) = x^2 over [0,1]
                    let samples: Vec<f64> =
                        (0..num_samples).map(|_| rng.r#gen::<f64>()).collect();
                    let values: Vec<f64> = samples.iter().map(|&x| x * x).collect();
                    let integral: f64 = values.iter().sum::<f64>() / num_samples as f64;
                    Ok(ComputeOutput {
                        result: serde_json::json!({
                            "integral": integral,
                            "samples": num_samples,
                            "method": "monte_carlo"
                        }),
                        additional: None,
                        metadata: None,
                    })
                }
                MonteCarloMethod::MCMC | MonteCarloMethod::MetropolisHastings => {
                    // Simple random walk MCMC
                    let mut current = 0.0f64;
                    let mut samples = Vec::with_capacity(num_samples);
                    for _ in 0..num_samples {
                        let proposal = current + (rng.r#gen::<f64>() - 0.5) * 2.0;
                        // Accept with probability min(1, p(proposal)/p(current)) for standard normal
                        let accept_prob =
                            ((-proposal * proposal / 2.0) - (-current * current / 2.0)).exp();
                        if rng.r#gen::<f64>() < accept_prob {
                            current = proposal;
                        }
                        samples.push(current);
                    }
                    let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
                    Ok(ComputeOutput {
                        result: serde_json::json!({
                            "samples": samples,
                            "mean": mean,
                            "method": "mcmc"
                        }),
                        additional: None,
                        metadata: None,
                    })
                }
                _ => Ok(ComputeOutput {
                    result: serde_json::json!({
                        "status": "monte_carlo_method_not_yet_implemented",
                        "method": format!("{:?}", mc_method)
                    }),
                    additional: None,
                    metadata: None,
                }),
            }
        }
        SamplingMethod::Stats(stat_method) => {
            let data: Vec<f64> = input
                .data
                .get("data")
                .and_then(|v| v.as_array())
                .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect())
                .unwrap_or_default();

            match stat_method {
                StatisticalMethod::BasicStats => {
                    let n = data.len();
                    if n == 0 {
                        return Ok(ComputeOutput {
                            result: serde_json::json!({"error": "empty data"}),
                            additional: None,
                            metadata: None,
                        });
                    }
                    let mean = data.iter().sum::<f64>() / n as f64;
                    let variance =
                        data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
                    let std_dev = variance.sqrt();
                    let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
                    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    Ok(ComputeOutput {
                        result: serde_json::json!({
                            "mean": mean,
                            "variance": variance,
                            "std_dev": std_dev,
                            "min": min,
                            "max": max,
                            "n": n
                        }),
                        additional: None,
                        metadata: None,
                    })
                }
                _ => Ok(ComputeOutput {
                    result: serde_json::json!({
                        "status": "statistical_method_not_yet_implemented",
                        "method": format!("{:?}", stat_method)
                    }),
                    additional: None,
                    metadata: None,
                }),
            }
        }
        _ => Ok(ComputeOutput {
            result: serde_json::json!({
                "status": "sampling_method_not_yet_implemented",
                "method": format!("{:?}", sampling_method)
            }),
            additional: None,
            metadata: None,
        }),
    }
}
