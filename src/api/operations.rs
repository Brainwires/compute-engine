//! List of all available operations across all modules

use std::collections::HashMap;

/// List all available operations across all modules
pub fn list_all_operations() -> HashMap<String, Vec<String>> {
    let mut ops = HashMap::new();

    // TENSOR CALCULUS - 14 operations
    ops.insert("tensor_calculus".to_string(), vec![
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
    ]);

    // ADVANCED CALCULUS - 22 operations
    ops.insert("advanced_calculus".to_string(), vec![
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
        "vector_field_analysis".to_string(),
        "scalar_field_analysis".to_string(),
        "tensor_field_analysis".to_string(),
        "differential_forms".to_string(),
        "exterior_derivative".to_string(),
    ]);

    // FLUID DYNAMICS - 22 operations
    ops.insert("fluid_dynamics".to_string(), vec![
        "cavity_flow".to_string(),
        "navier_stokes".to_string(),
        "channel_flow".to_string(),
        "lid_driven_cavity".to_string(),
        "all_no_slip".to_string(),
        "apply_velocity_boundary".to_string(),
        "apply_pressure_boundary".to_string(),
        "uniform".to_string(),
        "analyze".to_string(),
        "solve".to_string(),
        "new".to_string(),
        "cfl_condition".to_string(),
        "diffusion_condition".to_string(),
        "stable_timestep".to_string(),
        "is_boundary".to_string(),
        "is_interior".to_string(),
        "is_u_boundary".to_string(),
        "is_v_boundary".to_string(),
        "x_coord".to_string(),
        "y_coord".to_string(),
        "x_coordinates".to_string(),
        "y_coordinates".to_string(),
    ]);

    // SIGNAL PROCESSING - 9 operations
    ops.insert("signal_processing".to_string(), vec![
        "autocorrelation".to_string(),
        "power_spectrum".to_string(),
        "coherence".to_string(),
        "cepstrum".to_string(),
        "fft".to_string(),
        "ifft".to_string(),
        "filter".to_string(),
        "spectral_analysis".to_string(),
        "detect_peaks".to_string(),
        "compute_fourier_series".to_string(),
        "wavelet_transform".to_string(),
        "windowing_functions".to_string(),
        "cross_correlation".to_string(),
    ]);

    // STOCHASTIC PROCESSES - 14 operations
    ops.insert("stochastic_processes".to_string(), vec![
        "brownian_motion".to_string(),
        "geometric_brownian_motion".to_string(),
        "ornstein_uhlenbeck".to_string(),
        "poisson_process".to_string(),
        "levy_process".to_string(),
        "jump_diffusion".to_string(),
        "fractional_brownian_motion".to_string(),
        "mean_reverting_process".to_string(),
        "variance_gamma_process".to_string(),
        "heston_model".to_string(),
        "sabr_model".to_string(),
        "stochastic_volatility".to_string(),
        "path_generation".to_string(),
        "moment_calculation".to_string(),
    ]);

    // CRYPTOGRAPHIC MATHEMATICS - 10 operations
    ops.insert("cryptographic_mathematics".to_string(), vec![
        "generate_prime".to_string(),
        "is_prime".to_string(),
        "mod_exp".to_string(),
        "mod_inv".to_string(),
        "gcd".to_string(),
        "lcm".to_string(),
        "euler_totient".to_string(),
        "carmichael_lambda".to_string(),
        "discrete_log".to_string(),
        "elliptic_curve_point_add".to_string(),
    ]);

    // SYMBOLIC REGRESSION - 8 operations
    ops.insert("symbolic_regression".to_string(), vec![
        "fit_polynomial".to_string(),
        "fit_exponential".to_string(),
        "fit_logarithmic".to_string(),
        "fit_power_law".to_string(),
        "fit_rational".to_string(),
        "symbolic_differentiation".to_string(),
        "symbolic_integration".to_string(),
        "expression_simplification".to_string(),
    ]);

    // DIMENSIONAL ANALYSIS - 9 operations
    ops.insert("dimensional_analysis".to_string(), vec![
        "check_dimensional_consistency".to_string(),
        "derive_units".to_string(),
        "convert_units".to_string(),
        "buckingham_pi".to_string(),
        "dimensionless_groups".to_string(),
        "unit_analysis".to_string(),
        "validate_equation_dimensions".to_string(),
        "infer_dimensions".to_string(),
        "scale_analysis".to_string(),
    ]);

    // EQUATION VALIDATION - 8 operations
    ops.insert("equation_validation".to_string(), vec![
        "validate_equation".to_string(),
        "parse_equation".to_string(),
        "extract_variables".to_string(),
        "check_mathematical_correctness".to_string(),
        "check_dimensional_consistency".to_string(),
        "check_physics_compliance".to_string(),
        "check_conservation_laws".to_string(),
        "check_symmetries".to_string(),
    ]);

    // COMPUTATIONAL GEOMETRY - 5 operations
    ops.insert("computational_geometry".to_string(), vec![
        "convex_hull".to_string(),
        "delaunay_triangulation".to_string(),
        "voronoi_diagram".to_string(),
        "polygon_area".to_string(),
        "point_in_polygon".to_string(),
    ]);

    // LINEAR ALGEBRA - 6 operations
    ops.insert("linear_algebra".to_string(), vec![
        "matrix_norm".to_string(),
        "matrix_power".to_string(),
        "matrix_exp".to_string(),
        "cholesky".to_string(),
        "compute_qr".to_string(),
        "compute_svd".to_string(),
        "compute_eigendecomposition".to_string(),
        "compute_pca".to_string(),
        "compute_matrix_rank".to_string(),
        "compute_pseudoinverse".to_string(),
    ]);

    // STATISTICS - 6 operations
    ops.insert("statistics".to_string(), vec![
        "hypothesis_test".to_string(),
        "anova".to_string(),
        "regression_analysis".to_string(),
        "time_series".to_string(),
        "statistics".to_string(),
        "monte_carlo_integration".to_string(),
        "mcmc_sampling".to_string(),
        "correlation".to_string(),
        "kl_divergence".to_string(),
        "mutual_information".to_string(),
    ]);

    // OPTIMIZATION - 3 operations
    ops.insert("optimization".to_string(), vec![
        "curve_fitting".to_string(),
        "gradient_descent".to_string(),
        "nelder_mead".to_string(),
    ]);

    // GRAPH THEORY - 5 operations
    ops.insert("graph_theory".to_string(), vec![
        "shortest_path".to_string(),
        "minimum_spanning_tree".to_string(),
        "connected_components".to_string(),
        "graph_properties".to_string(),
        "topological_sort".to_string(),
    ]);

    // INFORMATION THEORY - 7 operations
    ops.insert("information_theory".to_string(), vec![
        "shannon_entropy".to_string(),
        "mutual_information".to_string(),
        "channel_capacity".to_string(),
        "huffman_coding".to_string(),
        "kolmogorov_complexity".to_string(),
        "conditional_entropy".to_string(),
        "relative_entropy".to_string(),
    ]);

    // ELECTROMAGNETISM - 8 operations
    ops.insert("electromagnetism".to_string(), vec![
        "maxwell_equations".to_string(),
        "em_wave".to_string(),
        "antenna".to_string(),
        "transmission_line".to_string(),
        "waveguide".to_string(),
        "scattering".to_string(),
        "poynting_vector".to_string(),
        "skin_effect".to_string(),
    ]);

    // SPECIAL FUNCTIONS - 6 operations
    ops.insert("special_functions".to_string(), vec![
        "bessel".to_string(),
        "gamma".to_string(),
        "error_function".to_string(),
        "elliptic".to_string(),
        "orthogonal_polynomial".to_string(),
        "airy".to_string(),
    ]);

    // NUMERICAL METHODS - 8 operations
    ops.insert("numerical_methods".to_string(), vec![
        "solve_ode".to_string(),
        "find_root".to_string(),
        "integrate".to_string(),
        "interpolate".to_string(),
        "solve_linear_system".to_string(),
        "differentiate".to_string(),
        "solve_pde".to_string(),
        "matrix_decomposition".to_string(),
    ]);

    // CHEMISTRY - 7 operations
    ops.insert("chemistry".to_string(), vec![
        "balance_equation".to_string(),
        "thermodynamics".to_string(),
        "electrochemistry".to_string(),
        "kinetics".to_string(),
        "gas_law".to_string(),
        "acid_base".to_string(),
        "molar_mass".to_string(),
    ]);

    ops
}
