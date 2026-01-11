//! List of all available operations across all modules
//!
//! This provides a comprehensive list of operations for discovery/documentation purposes.

use std::collections::HashMap;

/// Returns a list of all available tools (8 tools)
pub fn list_tools() -> Vec<&'static str> {
    vec![
        "solve", "compute", "analyze", "simulate", "ml", "chaos", "units", "validate",
    ]
}

/// List all available operations across all modules
pub fn list_all_operations() -> HashMap<String, Vec<String>> {
    let mut ops = HashMap::new();

    // PRIMARY TOOLS
    ops.insert(
        "solve".to_string(),
        vec![
            "root_finding".to_string(),
            "linear_system".to_string(),
            "polynomial".to_string(),
            "differential_ode".to_string(),
            "differential_pde".to_string(),
            "differential_bvp".to_string(),
            "algebraic_symbolic".to_string(),
            "algebraic_numeric".to_string(),
            "einstein_vacuum".to_string(),
            "einstein_with_source".to_string(),
            "einstein_axisymmetric".to_string(),
            "navier_stokes".to_string(),
            "fluid_euler".to_string(),
            "fluid_potential".to_string(),
            "heat_equation".to_string(),
            "wave_equation".to_string(),
            "schrodinger".to_string(),
            "maxwell".to_string(),
            "optimize_minimize".to_string(),
            "optimize_maximize".to_string(),
            "optimize_constrained".to_string(),
            "optimize_least_squares".to_string(),
            "game_theory_nash".to_string(),
            "game_theory_minimax".to_string(),
            "linear_programming_simplex".to_string(),
            "linear_programming_dual".to_string(),
        ],
    );

    ops.insert(
        "compute".to_string(),
        vec![
            // Matrix operations
            "matrix_determinant".to_string(),
            "matrix_inverse".to_string(),
            "matrix_eigenvalues".to_string(),
            "matrix_eigenvectors".to_string(),
            "matrix_svd".to_string(),
            "matrix_qr".to_string(),
            "matrix_lu".to_string(),
            "matrix_cholesky".to_string(),
            "matrix_norm".to_string(),
            "matrix_rank".to_string(),
            "matrix_trace".to_string(),
            "matrix_exp".to_string(),
            // Calculus
            "differentiate_symbolic".to_string(),
            "differentiate_numeric".to_string(),
            "differentiate_partial".to_string(),
            "integrate_definite".to_string(),
            "integrate_indefinite".to_string(),
            "integrate_numeric".to_string(),
            "integrate_improper".to_string(),
            "integrate_multiple".to_string(),
            // Transforms
            "transform_fourier".to_string(),
            "transform_laplace".to_string(),
            "transform_z".to_string(),
            "transform_wavelet".to_string(),
            "transform_hilbert".to_string(),
            // Field theory
            "field_electric".to_string(),
            "field_magnetic".to_string(),
            "field_gravitational".to_string(),
            "field_green_function".to_string(),
            "field_decoherence_scale".to_string(),
            "field_bohm_potential".to_string(),
            // Sampling
            "sample_monte_carlo".to_string(),
            "sample_mcmc".to_string(),
            "sample_importance".to_string(),
            "sample_rejection".to_string(),
            // Special functions
            "special_bessel".to_string(),
            "special_gamma".to_string(),
            "special_erf".to_string(),
            "special_elliptic".to_string(),
            "special_airy".to_string(),
            "special_legendre".to_string(),
            // Tensor operations
            "tensor_christoffel".to_string(),
            "tensor_riemann".to_string(),
            "tensor_ricci".to_string(),
            "tensor_weyl".to_string(),
            "tensor_einstein".to_string(),
        ],
    );

    ops.insert(
        "analyze".to_string(),
        vec![
            "simplify".to_string(),
            "expand".to_string(),
            "factor".to_string(),
            "series_taylor".to_string(),
            "series_laurent".to_string(),
            "series_fourier".to_string(),
            "limit".to_string(),
            "asymptotic".to_string(),
            "roots".to_string(),
            "extrema".to_string(),
            "stability_eigenvalue".to_string(),
            "stability_lyapunov".to_string(),
            "stability_routh_hurwitz".to_string(),
            "validate_expression".to_string(),
        ],
    );

    ops.insert(
        "simulate".to_string(),
        vec![
            // Time evolution
            "time_evolution_euler".to_string(),
            "time_evolution_runge_kutta".to_string(),
            "time_evolution_adaptive".to_string(),
            "time_evolution_symplectic".to_string(),
            // Stochastic
            "stochastic_brownian".to_string(),
            "stochastic_geometric_brownian".to_string(),
            "stochastic_ornstein_uhlenbeck".to_string(),
            "stochastic_jump_diffusion".to_string(),
            // Finance
            "finance_black_scholes".to_string(),
            "finance_heston".to_string(),
            "finance_sabr".to_string(),
            // Fluid dynamics
            "fluid_navier_stokes_2d".to_string(),
            "fluid_lattice_boltzmann".to_string(),
            "fluid_quantum_navier_stokes".to_string(),
            // N-body
            "nbody_gravitational".to_string(),
            "nbody_molecular".to_string(),
            // Cellular automata
            "cellular_automaton".to_string(),
        ],
    );

    // NEW TOOLS (from 8-tool expansion)
    ops.insert(
        "ml".to_string(),
        vec![
            "clustering_kmeans".to_string(),
            "clustering_dbscan".to_string(),
            "clustering_hierarchical".to_string(),
            "clustering_gaussian_mixture".to_string(),
            "regression_linear".to_string(),
            "regression_logistic".to_string(),
            "regression_ridge".to_string(),
            "regression_lasso".to_string(),
            "neural_network_create".to_string(),
            "neural_network_train".to_string(),
            "neural_network_predict".to_string(),
            "dimensionality_reduction_pca".to_string(),
            "dimensionality_reduction_tsne".to_string(),
            "classification_svm".to_string(),
            "classification_decision_tree".to_string(),
        ],
    );

    ops.insert(
        "chaos".to_string(),
        vec![
            "fractal_mandelbrot".to_string(),
            "fractal_julia".to_string(),
            "fractal_burning_ship".to_string(),
            "fractal_koch".to_string(),
            "fractal_sierpinski".to_string(),
            "attractor_lorenz".to_string(),
            "attractor_rossler".to_string(),
            "attractor_henon".to_string(),
            "attractor_logistic".to_string(),
            "lyapunov_exponent".to_string(),
            "lyapunov_spectrum".to_string(),
            "bifurcation_logistic".to_string(),
            "bifurcation_period_doubling".to_string(),
            "dimension_box_counting".to_string(),
            "dimension_correlation".to_string(),
            "dimension_kaplan_yorke".to_string(),
        ],
    );

    ops.insert(
        "units".to_string(),
        vec![
            "convert".to_string(),
            "analyze".to_string(),
            "check_compatibility".to_string(),
            "get_base".to_string(),
            "derive".to_string(),
            "parse".to_string(),
            "simplify".to_string(),
        ],
    );

    ops.insert(
        "validate".to_string(),
        vec![
            "equation".to_string(),
            "dimensions".to_string(),
            "conservation".to_string(),
            "symmetry".to_string(),
            "physics".to_string(),
            "bounds".to_string(),
            "singularities".to_string(),
        ],
    );

    // DOMAIN MODULES (for reference/discovery)
    ops.insert(
        "tensor_calculus".to_string(),
        vec![
            "solve_vacuum_einstein".to_string(),
            "solve_einstein_with_source".to_string(),
            "compute_christoffel_symbols".to_string(),
            "compute_riemann_tensor".to_string(),
            "compute_ricci_tensor".to_string(),
            "compute_ricci_scalar".to_string(),
            "compute_einstein_tensor".to_string(),
            "compute_weyl_tensor".to_string(),
            "tensor_product".to_string(),
            "tensor_contraction".to_string(),
            "covariant_derivative".to_string(),
            "lie_derivative".to_string(),
            "parallel_transport".to_string(),
            "geodesic_equation".to_string(),
        ],
    );

    ops.insert(
        "advanced_calculus".to_string(),
        vec![
            "calculus_of_variations".to_string(),
            "green_function".to_string(),
            "partial_fraction".to_string(),
            "series_expansion".to_string(),
            "limit_calculation".to_string(),
            "laplacian".to_string(),
            "gradient".to_string(),
            "divergence".to_string(),
            "curl".to_string(),
            "line_integral".to_string(),
            "surface_integral".to_string(),
            "volume_integral".to_string(),
            "greens_theorem".to_string(),
            "stokes_theorem".to_string(),
            "divergence_theorem".to_string(),
            "fourier_transform".to_string(),
            "inverse_fourier_transform".to_string(),
            "laplace_transform".to_string(),
            "inverse_laplace_transform".to_string(),
            "residue_theorem".to_string(),
            "laurent_series".to_string(),
            "conformal_mapping".to_string(),
        ],
    );

    ops.insert(
        "signal_processing".to_string(),
        vec![
            "fft".to_string(),
            "ifft".to_string(),
            "filter".to_string(),
            "spectral_analysis".to_string(),
            "autocorrelation".to_string(),
            "cross_correlation".to_string(),
            "power_spectrum".to_string(),
            "wavelet_transform".to_string(),
            "window_functions".to_string(),
        ],
    );

    ops.insert(
        "stochastic_processes".to_string(),
        vec![
            "brownian_motion".to_string(),
            "geometric_brownian_motion".to_string(),
            "ornstein_uhlenbeck".to_string(),
            "poisson_process".to_string(),
            "levy_process".to_string(),
            "jump_diffusion".to_string(),
            "heston_model".to_string(),
            "sabr_model".to_string(),
        ],
    );

    ops.insert(
        "cryptographic_mathematics".to_string(),
        vec![
            "generate_prime".to_string(),
            "is_prime".to_string(),
            "mod_exp".to_string(),
            "mod_inv".to_string(),
            "gcd".to_string(),
            "lcm".to_string(),
            "euler_totient".to_string(),
            "discrete_log".to_string(),
        ],
    );

    ops.insert(
        "linear_algebra".to_string(),
        vec![
            "matrix_norm".to_string(),
            "matrix_power".to_string(),
            "matrix_exp".to_string(),
            "cholesky".to_string(),
            "qr_decomposition".to_string(),
            "svd".to_string(),
            "eigendecomposition".to_string(),
            "pseudoinverse".to_string(),
        ],
    );

    ops.insert(
        "statistics".to_string(),
        vec![
            "descriptive_statistics".to_string(),
            "hypothesis_test".to_string(),
            "anova".to_string(),
            "regression_analysis".to_string(),
            "correlation".to_string(),
            "monte_carlo_integration".to_string(),
            "mcmc_sampling".to_string(),
        ],
    );

    ops.insert(
        "graph_theory".to_string(),
        vec![
            "shortest_path".to_string(),
            "minimum_spanning_tree".to_string(),
            "connected_components".to_string(),
            "topological_sort".to_string(),
            "graph_properties".to_string(),
        ],
    );

    ops.insert(
        "information_theory".to_string(),
        vec![
            "shannon_entropy".to_string(),
            "mutual_information".to_string(),
            "channel_capacity".to_string(),
            "huffman_coding".to_string(),
            "conditional_entropy".to_string(),
            "relative_entropy".to_string(),
        ],
    );

    ops.insert(
        "electromagnetism".to_string(),
        vec![
            "maxwell_equations".to_string(),
            "em_wave".to_string(),
            "antenna".to_string(),
            "transmission_line".to_string(),
            "waveguide".to_string(),
            "scattering".to_string(),
            "poynting_vector".to_string(),
        ],
    );

    ops.insert(
        "special_functions".to_string(),
        vec![
            "bessel".to_string(),
            "gamma".to_string(),
            "error_function".to_string(),
            "elliptic".to_string(),
            "orthogonal_polynomial".to_string(),
            "airy".to_string(),
        ],
    );

    ops.insert(
        "numerical_methods".to_string(),
        vec![
            "solve_ode".to_string(),
            "find_root".to_string(),
            "integrate".to_string(),
            "interpolate".to_string(),
            "solve_linear_system".to_string(),
            "differentiate".to_string(),
            "solve_pde".to_string(),
        ],
    );

    ops.insert(
        "chemistry".to_string(),
        vec![
            "balance_equation".to_string(),
            "thermodynamics".to_string(),
            "electrochemistry".to_string(),
            "kinetics".to_string(),
            "gas_law".to_string(),
            "acid_base".to_string(),
            "molar_mass".to_string(),
        ],
    );

    ops
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_list_tools() {
        let tools = list_tools();
        assert_eq!(tools.len(), 8);
        assert!(tools.contains(&"solve"));
        assert!(tools.contains(&"compute"));
        assert!(tools.contains(&"analyze"));
        assert!(tools.contains(&"simulate"));
        assert!(tools.contains(&"ml"));
        assert!(tools.contains(&"chaos"));
        assert!(tools.contains(&"units"));
        assert!(tools.contains(&"validate"));
    }

    #[test]
    fn test_list_all_operations() {
        let ops = list_all_operations();
        assert!(ops.contains_key("solve"));
        assert!(ops.contains_key("compute"));
        assert!(ops.contains_key("analyze"));
        assert!(ops.contains_key("simulate"));
        assert!(ops.contains_key("ml"));
        assert!(ops.contains_key("chaos"));
        assert!(ops.contains_key("units"));
        assert!(ops.contains_key("validate"));
    }
}
