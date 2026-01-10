//! Comprehensive CLI help system that programmatically discovers all operations
//!
//! This module provides detailed help for 4 primary tools and 6 legacy tools.
//!
//! ## 4-Tool Architecture (Primary API)
//!
//! The engine has been consolidated to 4 primary tools:
//! - **Solve**: Equations, systems, root finding, optimization
//! - **Compute**: Calculus, transforms, field theory, sampling, matrix ops
//! - **Analyze**: Series, limits, stability, simplification
//! - **Simulate**: Time evolution, stochastic, fluid dynamics
//!
//! ## Legacy Tools (Backward Compatibility)
//!
//! Legacy tool names are still supported for backward compatibility:
//! - `differentiate` → routes to `compute` with `operation: {differentiate: ...}`
//! - `integrate` → routes to `compute` with `operation: {integrate: ...}`
//! - `transform` → routes to `compute` with `operation: {transform: ...}`
//! - `fieldtheory` → routes to `compute` with `operation: {field: ...}`
//! - `sample` → routes to `compute` with `operation: {sample: ...}`
//! - `optimize` → routes to `solve` with `equation_type: {optimize: ...}`

use crate::engine::equations::*;
use std::fmt;

/// Tool names (primary 4 + legacy 6)
#[derive(Debug, Clone, Copy)]
pub enum ToolName {
    // Primary 4 tools
    Solve,
    Compute,
    Analyze,
    Simulate,
    // Legacy tools (backward compatibility)
    Differentiate,
    Integrate,
    Transform,
    FieldTheory,
    Sample,
    Optimize,
}

impl fmt::Display for ToolName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ToolName::Solve => write!(f, "solve"),
            ToolName::Differentiate => write!(f, "differentiate"),
            ToolName::Integrate => write!(f, "integrate"),
            ToolName::Analyze => write!(f, "analyze"),
            ToolName::Simulate => write!(f, "simulate"),
            ToolName::Compute => write!(f, "compute"),
            ToolName::Transform => write!(f, "transform"),
            ToolName::FieldTheory => write!(f, "fieldtheory"),
            ToolName::Sample => write!(f, "sample"),
            ToolName::Optimize => write!(f, "optimize"),
        }
    }
}

impl ToolName {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "solve" => Some(ToolName::Solve),
            "differentiate" | "diff" => Some(ToolName::Differentiate),
            "integrate" | "int" => Some(ToolName::Integrate),
            "analyze" => Some(ToolName::Analyze),
            "simulate" | "sim" => Some(ToolName::Simulate),
            "compute" | "calc" | "calculate" => Some(ToolName::Compute),
            "transform" => Some(ToolName::Transform),
            "fieldtheory" | "field" => Some(ToolName::FieldTheory),
            "sample" => Some(ToolName::Sample),
            "optimize" | "opt" => Some(ToolName::Optimize),
            _ => None,
        }
    }

    /// Returns primary 4 tools only
    pub fn primary() -> Vec<ToolName> {
        vec![
            ToolName::Solve,
            ToolName::Compute,
            ToolName::Analyze,
            ToolName::Simulate,
        ]
    }

    /// Returns legacy tools (backward compatibility)
    pub fn legacy() -> Vec<ToolName> {
        vec![
            ToolName::Differentiate,
            ToolName::Integrate,
            ToolName::Transform,
            ToolName::FieldTheory,
            ToolName::Sample,
            ToolName::Optimize,
        ]
    }

    /// Returns all tools (primary + legacy)
    pub fn all() -> Vec<ToolName> {
        let mut all = Self::primary();
        all.extend(Self::legacy());
        all
    }

    pub fn description(&self) -> &str {
        match self {
            // Primary 4 tools
            ToolName::Solve => "Solve equations, systems, optimization, and mathematical problems",
            ToolName::Compute => "Compute: calculus, transforms, fields, sampling, matrix ops",
            ToolName::Analyze => "Analyze expressions (validate, simplify, parse, stability)",
            ToolName::Simulate => "Simulate time evolution, stochastic, and fluid dynamics",
            // Legacy tools (backward compatibility)
            ToolName::Differentiate => "[Legacy → compute] Derivatives and gradients",
            ToolName::Integrate => "[Legacy → compute] Integrals (definite, indefinite, etc.)",
            ToolName::Transform => "[Legacy → compute] Fourier, Laplace, wavelet transforms",
            ToolName::FieldTheory => "[Legacy → compute] Physics fields (EM, quantum, gravity)",
            ToolName::Sample => "[Legacy → compute] Statistical sampling and analysis",
            ToolName::Optimize => "[Legacy → solve] Curve fitting and optimization",
        }
    }

    pub fn is_legacy(&self) -> bool {
        matches!(
            self,
            ToolName::Differentiate
                | ToolName::Integrate
                | ToolName::Transform
                | ToolName::FieldTheory
                | ToolName::Sample
                | ToolName::Optimize
        )
    }
}

/// Print the main help screen (not currently used - kept for future reference)
#[allow(dead_code)]
pub fn print_main_help() {
    println!("Computational Engine - Consolidated 4-Tool Architecture");
    println!();
    println!("USAGE:");
    println!("  brainwires-compute-engine list-ops [TOOL]");
    println!("  brainwires-compute-engine <SUBCOMMAND>");
    println!();
    println!("DESCRIPTION:");
    println!("  A unified computational engine providing 194+ mathematical operations");
    println!("  accessible through a clean 4-tool API with backward compatibility.");
    println!();
    println!("PRIMARY TOOLS (Recommended):");
    for tool in ToolName::primary() {
        println!("  {:15} - {}", tool.to_string(), tool.description());
    }
    println!();
    println!("LEGACY TOOLS (Backward Compatibility):");
    for tool in ToolName::legacy() {
        println!("  {:15} - {}", tool.to_string(), tool.description());
    }
    println!();
    println!("USAGE EXAMPLES:");
    println!("  brainwires-compute-engine list-ops             # Show ALL operations by module");
    println!("  brainwires-compute-engine list-ops compute     # Show compute tool operations");
    println!("  brainwires-compute-engine list-ops solve       # Show solve tool operations");
    println!("  brainwires-compute-engine help                 # Show built-in help");
    println!("  brainwires-compute-engine info                 # Show version and modules");
    println!();
    println!("For detailed help on a specific tool, use:");
    println!("  brainwires-compute-engine list-ops <TOOL>");
}

/// Print detailed help for a specific tool
pub fn print_tool_help(tool: ToolName) {
    match tool {
        ToolName::Solve => print_solve_help(),
        ToolName::Differentiate => print_differentiate_help(),
        ToolName::Integrate => print_integrate_help(),
        ToolName::Analyze => print_analyze_help(),
        ToolName::Simulate => print_simulate_help(),
        ToolName::Compute => print_compute_help(),
        ToolName::Transform => print_transform_help(),
        ToolName::FieldTheory => print_fieldtheory_help(),
        ToolName::Sample => print_sample_help(),
        ToolName::Optimize => print_optimize_help(),
    }
}

fn print_solve_help() {
    println!("TOOL: solve");
    println!("DESCRIPTION: Solve equations, systems, and mathematical problems");
    println!();
    println!("EQUATION TYPES:");
    println!();

    println!("  Einstein Equations:");
    println!("    • vacuum              - Vacuum Einstein equations (Rμν = 0)");
    println!("    • with_source         - Einstein equations with energy-momentum tensor");
    println!("    • schwarzschild       - Schwarzschild black hole solution");
    println!("    • kerr_newman         - Kerr-Newman rotating charged black hole");
    println!("    • friedmann_robertson_walker - FRW cosmological metric");
    println!();

    println!("  Fluid Equations:");
    println!("    • navier_stokes       - Navier-Stokes equations for viscous flow");
    println!("    • cavity_flow         - Lid-driven cavity flow");
    println!("    • channel_flow        - Flow in a channel");
    println!("    • euler               - Inviscid Euler equations");
    println!("    • bernoulli           - Bernoulli equation for fluid flow");
    println!();

    println!("  Differential Equations:");
    println!("    • ode                 - Ordinary differential equations");
    println!("    • pde                 - Partial differential equations");
    println!("    • boundary_value      - Boundary value problems");
    println!("    • initial_value       - Initial value problems");
    println!();

    println!("  Electromagnetic Equations:");
    println!("    • maxwell             - Maxwell's equations");
    println!("    • wave                - Wave equations");
    println!("    • transmission_line   - Transmission line equations");
    println!("    • waveguide           - Waveguide modes");
    println!("    • helmholtz           - Helmholtz equation");
    println!();

    println!("  Chemical Equations:");
    println!("    • balance             - Chemical equation balancing");
    println!("    • thermodynamic       - Thermodynamic equilibrium");
    println!("    • kinetics            - Reaction kinetics");
    println!("    • gas_law             - Ideal and real gas laws");
    println!("    • acid_base           - Acid-base equilibria");
    println!("    • electrochemistry    - Electrochemical reactions");
    println!();

    println!("  Other Types:");
    println!("    • linear_system       - Systems of linear equations (Ax = b)");
    println!("    • root_finding        - Finding roots of equations");
    println!("    • number_theory       - Number theory problems (primality, factorization)");
    println!("    • differential_geometry - Geodesics, parallel transport, minimal surfaces");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "solve",
    "input": {{
      "equation_type": "linear_system",
      "equations": ["2*x + 3*y = 8", "x - y = 1"],
      "variables": ["x", "y"]
    }}
  }}"#);
}

fn print_differentiate_help() {
    println!("TOOL: differentiate");
    println!("DESCRIPTION: Compute derivatives, gradients, and differential operators");
    println!();
    println!("OPERATION TYPES:");
    println!();

    println!("  Vector Calculus:");
    println!("    • gradient            - ∇f (gradient of scalar field)");
    println!("    • divergence          - ∇·F (divergence of vector field)");
    println!("    • curl                - ∇×F (curl of vector field)");
    println!("    • laplacian           - ∇²f (Laplacian operator)");
    println!("    • directional         - Directional derivative");
    println!();

    println!("  Tensor Calculus:");
    println!("    • covariant           - Covariant derivative");
    println!("    • lie                 - Lie derivative");
    println!("    • exterior_derivative - Exterior derivative of differential forms");
    println!();

    println!("  Other Operations:");
    println!("    • variational         - Variational derivatives (Euler-Lagrange)");
    println!("    • differential_forms  - Operations on differential forms");
    println!("    • numeric             - Numeric differentiation");
    println!("    • symbolic            - Symbolic differentiation");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "differentiate",
    "input": {{
      "operation": {{"vector_calc": "gradient"}},
      "expression": "x^2 + y^2 + z^2",
      "variables": ["x", "y", "z"]
    }}
  }}"#);
}

fn print_integrate_help() {
    println!("TOOL: integrate");
    println!("DESCRIPTION: Compute integrals (definite, indefinite, line, surface, volume)");
    println!();
    println!("INTEGRATION TYPES:");
    println!();

    println!("  Geometric Integrals:");
    println!("    • line                - Line integrals ∫C F·dr");
    println!("    • surface             - Surface integrals ∬S F·dS");
    println!("    • volume              - Volume integrals ∭V f dV");
    println!("    • contour             - Contour integrals in complex plane");
    println!();

    println!("  Integral Theorems:");
    println!("    • greens              - Green's theorem");
    println!("    • stokes              - Stokes' theorem");
    println!("    • divergence          - Divergence theorem (Gauss)");
    println!("    • cauchy_integral     - Cauchy integral formula");
    println!();

    println!("  Complex Analysis:");
    println!("    • residue             - Residue theorem");
    println!("    • cauchy              - Cauchy integral");
    println!("    • contour             - Contour integration");
    println!();

    println!("  Numeric Integration:");
    println!("    • trapezoidal         - Trapezoidal rule");
    println!("    • simpson             - Simpson's rule");
    println!("    • gauss_quadrature    - Gaussian quadrature");
    println!("    • adaptive            - Adaptive integration");
    println!();

    println!("  Other:");
    println!("    • symbolic            - Symbolic integration");
    println!("    • monte_carlo         - Monte Carlo integration");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "integrate",
    "input": {{
      "integration_type": {{"numeric": "gauss_quadrature"}},
      "expression": "sin(x)",
      "variables": ["x"],
      "limits": [[0.0, 3.14159]]
    }}
  }}"#);
}

fn print_analyze_help() {
    println!("TOOL: analyze");
    println!("DESCRIPTION: Analyze expressions (validate, simplify, parse, check)");
    println!();
    println!("OPERATIONS:");
    println!();

    println!("  Expression Operations:");
    println!("    • simplify            - Algebraic simplification");
    println!("    • parse               - Parse mathematical expressions");
    println!("    • extract_variables   - Extract variables from expressions");
    println!();

    println!("  Validation:");
    println!("    • validate            - General validation");
    println!("    • check_correctness   - Check mathematical correctness");
    println!("    • check_dimensions    - Dimensional analysis");
    println!("    • check_physics       - Physics consistency checks");
    println!("    • check_conservation  - Conservation law verification");
    println!("    • check_symmetries    - Symmetry analysis");
    println!();

    println!("  Series and Limits:");
    println!("    • partial_fraction    - Partial fraction decomposition");
    println!("    • series_expansion    - Taylor/Maclaurin series");
    println!("    • laurent_series      - Laurent series expansion");
    println!("    • limit               - Limit computation");
    println!();

    println!("  Field Analysis:");
    println!("    • field_analysis      - Analyze vector/scalar/tensor fields");
    println!();

    println!("  Dimensional Analysis:");
    println!("    • dimensional_check   - Check dimensional consistency");
    println!("    • validate_dimensions - Validate dimensions");
    println!("    • infer_dimensions    - Infer unknown dimensions");
    println!("    • scale_analysis      - Scaling analysis");
    println!();

    println!("  Other:");
    println!("    • units_derive        - Derive units");
    println!("    • units_analyze       - Analyze units");
    println!("    • graph_components    - Graph component analysis");
    println!("    • graph_properties    - Graph property analysis");
    println!("    • is_prime            - Primality testing");
    println!("    • fluid_analysis      - Fluid dynamics analysis");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "analyze",
    "input": {{
      "operation": "simplify",
      "expression": "(x^2 - 1)/(x - 1)"
    }}
  }}"#);
}

fn print_simulate_help() {
    println!("TOOL: simulate");
    println!("DESCRIPTION: Simulate time evolution and stochastic processes");
    println!();
    println!("SIMULATION MODELS:");
    println!();

    println!("  Stochastic Processes:");
    println!("    • brownian_motion     - Standard Brownian motion");
    println!("    • geometric_brownian  - Geometric Brownian motion (stock prices)");
    println!("    • ornstein_uhlenbeck  - Mean-reverting process");
    println!("    • poisson             - Poisson process");
    println!("    • levy                - Lévy process");
    println!("    • jump_diffusion      - Jump-diffusion process");
    println!("    • fractional_brownian - Fractional Brownian motion");
    println!("    • mean_reverting      - Mean-reverting processes");
    println!("    • variance_gamma      - Variance-gamma process");
    println!();

    println!("  Finance Models:");
    println!("    • heston              - Heston stochastic volatility");
    println!("    • sabr                - SABR volatility model");
    println!("    • stochastic_volatility - General stochastic volatility");
    println!("    • black_scholes       - Black-Scholes model");
    println!();

    println!("  Fluid Dynamics:");
    println!("    • navier_stokes       - Navier-Stokes simulation");
    println!("    • euler               - Euler equations");
    println!("    • lattice_boltzmann   - Lattice Boltzmann method");
    println!();

    println!("  Time Evolution:");
    println!("    • euler               - Euler method");
    println!("    • runge_kutta4        - 4th-order Runge-Kutta");
    println!("    • adaptive_step       - Adaptive step size");
    println!("    • implicit_euler      - Implicit Euler (stiff equations)");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "simulate",
    "input": {{
      "model": {{"stochastic": "brownian_motion"}},
      "equations": ["dx = sigma * dW"],
      "variables": ["x"],
      "parameters": {{"sigma": 1.0}},
      "initial_conditions": {{"x": 0.0}},
      "range": [0.0, 1.0],
      "steps": 1000,
      "num_paths": 100
    }}
  }}"#);
}

fn print_compute_help() {
    println!("TOOL: compute");
    println!("DESCRIPTION: Compute tensor operations, matrix decompositions, special functions");
    println!();
    println!("OPERATION CATEGORIES (100+ operations):");
    println!();

    println!("  Tensor Operations:");
    println!("    • christoffel         - Christoffel symbols");
    println!("    • riemann             - Riemann curvature tensor");
    println!("    • ricci               - Ricci tensor");
    println!("    • ricci_scalar        - Ricci scalar");
    println!("    • einstein            - Einstein tensor");
    println!("    • weyl                - Weyl tensor");
    println!("    • product             - Tensor products");
    println!("    • contraction         - Tensor contraction");
    println!("    • parallel_transport  - Parallel transport");
    println!();

    println!("  Matrix Operations:");
    println!("    • norm                - Matrix norms");
    println!("    • power               - Matrix powers");
    println!("    • exp                 - Matrix exponential");
    println!("    • rank                - Matrix rank");
    println!("    • pseudoinverse       - Moore-Penrose pseudoinverse");
    println!("    • determinant         - Determinant");
    println!("    • trace               - Trace");
    println!("    • inverse             - Matrix inverse");
    println!();

    println!("  Matrix Decompositions:");
    println!("    • qr                  - QR decomposition");
    println!("    • svd                 - Singular value decomposition");
    println!("    • eigen               - Eigenvalue decomposition");
    println!("    • pca                 - Principal component analysis");
    println!("    • cholesky            - Cholesky decomposition");
    println!("    • lu                  - LU decomposition");
    println!("    • schur               - Schur decomposition");
    println!();

    println!("  Special Functions:");
    println!("    • bessel              - Bessel functions (J, Y, I, K)");
    println!("    • gamma               - Gamma and related functions");
    println!("    • erf                 - Error function");
    println!("    • elliptic            - Elliptic integrals");
    println!("    • orthogonal_poly     - Orthogonal polynomials");
    println!("    • airy                - Airy functions");
    println!("    • hypergeometric      - Hypergeometric functions");
    println!();

    println!("  Number Theory:");
    println!("    • generate_prime      - Generate prime numbers");
    println!("    • mod_exp             - Modular exponentiation");
    println!("    • mod_inv             - Modular inverse");
    println!("    • gcd                 - Greatest common divisor");
    println!("    • lcm                 - Least common multiple");
    println!("    • euler_totient       - Euler's totient function");
    println!("    • carmichael_lambda   - Carmichael function");
    println!("    • ec_point_add        - Elliptic curve point addition");
    println!("    • rsa_keypair         - RSA key generation");
    println!("    • rsa_encrypt         - RSA encryption");
    println!("    • rsa_decrypt         - RSA decryption");
    println!("    • sha256              - SHA-256 hash");
    println!("    • sha3_256            - SHA3-256 hash");
    println!("    • chinese_remainder   - Chinese remainder theorem");
    println!("    • discrete_log        - Discrete logarithm");
    println!("    • primality_test      - Primality testing");
    println!();

    println!("  Computational Geometry:");
    println!("    • convex_hull         - Convex hull computation");
    println!("    • delaunay            - Delaunay triangulation");
    println!("    • voronoi             - Voronoi diagram");
    println!("    • polygon_area        - Area of polygon");
    println!("    • point_in_polygon    - Point containment test");
    println!();

    println!("  Information Theory:");
    println!("    • entropy             - Shannon entropy");
    println!("    • mutual_info         - Mutual information");
    println!("    • channel_capacity    - Channel capacity");
    println!("    • huffman             - Huffman coding");
    println!("    • kolmogorov          - Kolmogorov complexity");
    println!("    • conditional_entropy - Conditional entropy");
    println!("    • kl_divergence       - Kullback-Leibler divergence");
    println!();

    println!("  Physics (57 operations):");
    println!("    Relativity (10): lorentz_transform, time_dilation, gravitational_lensing...");
    println!("    Statistical Physics (10): partition_function, maxwell_boltzmann...");
    println!("    Quantum Mechanics (15): schrodinger_equation, entanglement_measure...");
    println!("    Control Systems (12): transfer_function, bode_plot, state_space...");
    println!("    Nuclear Physics (8): radioactive_decay, binding_energy, fusion_energy...");
    println!();

    println!("  Chemistry (10 operations):");
    println!("    • molar_mass, ph_calculation, arrhenius, nernst_equation...");
    println!();

    println!("  Biology (6 operations):");
    println!("    • michaelis_menten, hardy_weinberg, allometric_scaling...");
    println!();

    println!("  Thermodynamics (5 operations):");
    println!("    • conduction, convection, radiation, thermal_resistance...");
    println!();

    println!("  Optics (4 operations):");
    println!("    • thin_lens, snells_law, diffraction_grating, fresnel_equations");
    println!();

    println!("  Geophysics (4 operations):");
    println!("    • seismology, atmosphere, radiometric_dating, planetary_science");
    println!();

    println!("  Engineering (11 operations):");
    println!("    • sound_pressure_level, stress, pid_controller...");
    println!();

    println!("  Date/Time (11 operations):");
    println!("    • add_interval, business_days, is_leap_year...");
    println!();

    println!("  Graph Theory (3 operations):");
    println!("    • topological_sort, shortest_path, minimum_spanning_tree");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "compute",
    "input": {{
      "operation": {{"matrix_decomp": "svd"}},
      "data": [[1, 2], [3, 4], [5, 6]]
    }}
  }}"#);
}

fn print_transform_help() {
    println!("TOOL: transform");
    println!("DESCRIPTION: Apply transforms (Fourier, Laplace, wavelet, filters)");
    println!();
    println!("TRANSFORM TYPES:");
    println!();

    println!("  Fourier Transforms:");
    println!("    • forward             - Forward Fourier transform");
    println!("    • inverse             - Inverse Fourier transform");
    println!();

    println!("  Laplace Transforms:");
    println!("    • forward             - Forward Laplace transform");
    println!("    • inverse             - Inverse Laplace transform");
    println!();

    println!("  Wavelet Transforms:");
    println!("    • haar                - Haar wavelet");
    println!("    • daubechies          - Daubechies wavelets");
    println!("    • morlet              - Morlet wavelet");
    println!("    • mexican             - Mexican hat (Ricker) wavelet");
    println!();

    println!("  FFT:");
    println!("    • forward             - Forward FFT");
    println!("    • inverse             - Inverse FFT");
    println!();

    println!("  Filters:");
    println!("    • low_pass            - Low-pass filter");
    println!("    • high_pass           - High-pass filter");
    println!("    • band_pass           - Band-pass filter");
    println!("    • band_stop           - Band-stop filter");
    println!();

    println!("  Window Functions:");
    println!("    • hamming             - Hamming window");
    println!("    • hanning             - Hann window");
    println!("    • blackman            - Blackman window");
    println!("    • kaiser              - Kaiser window");
    println!();

    println!("  Other:");
    println!("    • conformal           - Conformal mappings");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "transform",
    "input": {{
      "transform_type": {{"fft": "forward"}},
      "data": [1.0, 2.0, 3.0, 4.0],
      "sampling_rate": 1000.0
    }}
  }}"#);
}

fn print_fieldtheory_help() {
    println!("TOOL: fieldtheory");
    println!("DESCRIPTION: Work with physics fields (EM, gravity, quantum)");
    println!();
    println!("FIELD TYPES:");
    println!();

    println!("  Electromagnetic Fields:");
    println!("    • antenna             - Antenna radiation patterns");
    println!("    • waveguide           - Waveguide field solutions");
    println!("    • scattering          - EM scattering problems");
    println!();

    println!("  Green's Functions:");
    println!("    • green_function      - Green's function solutions");
    println!();

    println!("  Quantum Field Theory:");
    println!("    • scalar_field        - Scalar field theory");
    println!("    • dirac_field         - Dirac (fermion) fields");
    println!("    • gauge_field         - Gauge field theories");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "fieldtheory",
    "input": {{
      "field_type": {{"em": "antenna"}},
      "configuration": {{
        "frequency": 2.4e9,
        "power": 1.0
      }},
      "points": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    }}
  }}"#);
}

fn print_sample_help() {
    println!("TOOL: sample");
    println!("DESCRIPTION: Sample and analyze statistical data");
    println!();
    println!("SAMPLING METHODS:");
    println!();

    println!("  Stochastic:");
    println!("    • path_generation     - Generate stochastic paths");
    println!("    • moments             - Compute statistical moments");
    println!();

    println!("  Monte Carlo:");
    println!("    • integration         - Monte Carlo integration");
    println!("    • mcmc                - Markov Chain Monte Carlo");
    println!("    • metropolis_hastings - Metropolis-Hastings algorithm");
    println!("    • gibbs               - Gibbs sampling");
    println!();

    println!("  Statistical Methods:");
    println!("    • hypothesis_test     - Statistical hypothesis testing");
    println!("    • anova               - Analysis of variance");
    println!("    • regression          - Regression analysis");
    println!("    • time_series         - Time series analysis");
    println!("    • basic_stats         - Mean, variance, std dev, etc.");
    println!("    • correlation         - Correlation analysis");
    println!();

    println!("  Signal Analysis:");
    println!("    • spectral_analysis   - Spectral analysis");
    println!("    • autocorrelation     - Autocorrelation function");
    println!("    • cross_correlation   - Cross-correlation");
    println!("    • power_spectrum      - Power spectral density");
    println!("    • coherence           - Coherence function");
    println!("    • cepstrum            - Cepstral analysis");
    println!("    • peak_detection      - Peak finding algorithms");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "sample",
    "input": {{
      "method": {{"stats": "basic_stats"}},
      "data": [1.2, 2.3, 3.1, 4.5, 5.0]
    }}
  }}"#);
}

fn print_optimize_help() {
    println!("TOOL: optimize");
    println!("DESCRIPTION: Optimize and fit data to models");
    println!();
    println!("OPTIMIZATION METHODS:");
    println!();

    println!("  Curve Fitting:");
    println!("    • polynomial          - Polynomial fitting");
    println!("    • exponential         - Exponential fitting (y = a*exp(bx))");
    println!("    • logarithmic         - Logarithmic fitting (y = a + b*ln(x))");
    println!("    • power_law           - Power law fitting (y = a*x^b)");
    println!("    • rational            - Rational function fitting");
    println!("    • trigonometric       - Trigonometric fitting");
    println!("    • custom              - Custom function fitting");
    println!();

    println!("  Minimization:");
    println!("    • gradient_descent    - Gradient descent");
    println!("    • nelder_mead         - Nelder-Mead simplex");
    println!("    • conjugate_gradient  - Conjugate gradient method");
    println!("    • bfgs                - BFGS quasi-Newton");
    println!("    • levenberg_marquardt - Levenberg-Marquardt algorithm");
    println!();

    println!("  Interpolation:");
    println!("    • linear              - Linear interpolation");
    println!("    • polynomial          - Polynomial interpolation");
    println!("    • spline              - Spline interpolation");
    println!("    • cubic               - Cubic interpolation");
    println!();

    println!("  Dimensional Analysis:");
    println!("    • buckingham_pi       - Buckingham Pi theorem");
    println!("    • dimensionless_groups - Find dimensionless groups");
    println!("    • similarity_analysis - Similarity analysis");
    println!();

    println!("  Symbolic Regression:");
    println!("    • symbolic_regression - Evolutionary function discovery");
    println!();

    println!("  Automatic Model Selection:");
    println!("    • auto                - Auto-select best fit using AIC/BIC/AICc");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "optimize",
    "input": {{
      "method": {{"fit": "polynomial"}},
      "data": [[1.0, 2.0, 3.0], [2.1, 4.2, 6.1]],
      "parameters": {{"degree": 2}}
    }}
  }}"#);
}

/// Count total operations across all tools
pub fn count_operations() -> usize {
    // This is an approximation - could be made more precise by
    // programmatically counting enum variants
    194
}
