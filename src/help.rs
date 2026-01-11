//! Comprehensive CLI help system
//!
//! ## 8-Tool Architecture
//!
//! The engine provides 8 tools:
//! - **Solve**: Equations, systems, root finding, optimization
//! - **Compute**: Calculus, transforms, field theory, sampling, matrix ops
//! - **Analyze**: Series, limits, stability, simplification
//! - **Simulate**: Time evolution, stochastic, fluid dynamics
//! - **ML**: Machine learning (clustering, neural nets, regression)
//! - **Chaos**: Chaos theory (fractals, attractors, Lyapunov)
//! - **Units**: Dimensional analysis and unit conversion
//! - **Validate**: Equation and physics validation

use std::fmt;

/// Tool names (8 tools)
#[derive(Debug, Clone, Copy)]
pub enum ToolName {
    Solve,
    Compute,
    Analyze,
    Simulate,
    ML,
    Chaos,
    Units,
    Validate,
}

impl fmt::Display for ToolName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ToolName::Solve => write!(f, "solve"),
            ToolName::Compute => write!(f, "compute"),
            ToolName::Analyze => write!(f, "analyze"),
            ToolName::Simulate => write!(f, "simulate"),
            ToolName::ML => write!(f, "ml"),
            ToolName::Chaos => write!(f, "chaos"),
            ToolName::Units => write!(f, "units"),
            ToolName::Validate => write!(f, "validate"),
        }
    }
}

impl ToolName {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "solve" => Some(ToolName::Solve),
            "compute" | "calc" | "calculate" => Some(ToolName::Compute),
            "analyze" => Some(ToolName::Analyze),
            "simulate" | "sim" => Some(ToolName::Simulate),
            "ml" | "machine_learning" => Some(ToolName::ML),
            "chaos" => Some(ToolName::Chaos),
            "units" => Some(ToolName::Units),
            "validate" => Some(ToolName::Validate),
            _ => None,
        }
    }

    /// Returns all 8 tools
    pub fn all() -> Vec<ToolName> {
        vec![
            ToolName::Solve,
            ToolName::Compute,
            ToolName::Analyze,
            ToolName::Simulate,
            ToolName::ML,
            ToolName::Chaos,
            ToolName::Units,
            ToolName::Validate,
        ]
    }

    pub fn description(&self) -> &str {
        match self {
            ToolName::Solve => "Solve equations, systems, optimization, and mathematical problems",
            ToolName::Compute => "Compute: calculus, transforms, fields, sampling, matrix ops",
            ToolName::Analyze => "Analyze expressions (validate, simplify, parse, stability)",
            ToolName::Simulate => "Simulate time evolution, stochastic, and fluid dynamics",
            ToolName::ML => "Machine learning: clustering, neural networks, regression",
            ToolName::Chaos => "Chaos theory: fractals, attractors, Lyapunov exponents",
            ToolName::Units => "Dimensional analysis and unit conversion",
            ToolName::Validate => "Equation and physics validation",
        }
    }
}

/// Print the main help screen
pub fn print_main_help() {
    println!("Computational Engine - 8-Tool Architecture");
    println!();
    println!("USAGE:");
    println!("  brainwires-compute-engine help [TOOL]");
    println!("  brainwires-compute-engine <SUBCOMMAND>");
    println!();
    println!("DESCRIPTION:");
    println!("  A unified computational engine providing 200+ mathematical operations");
    println!("  accessible through a clean 8-tool API.");
    println!();
    println!("TOOLS:");
    for tool in ToolName::all() {
        println!("  {:15} - {}", tool.to_string(), tool.description());
    }
    println!();
    println!("USAGE EXAMPLES:");
    println!("  brainwires-compute-engine help compute   # Show compute tool operations");
    println!("  brainwires-compute-engine help solve     # Show solve tool operations");
    println!("  brainwires-compute-engine info           # Show version and modules");
    println!();
    println!("For detailed help on a specific tool, use:");
    println!("  brainwires-compute-engine help <TOOL>");
}

/// Print detailed help for a specific tool
pub fn print_tool_help(tool: ToolName) {
    match tool {
        ToolName::Solve => print_solve_help(),
        ToolName::Compute => print_compute_help(),
        ToolName::Analyze => print_analyze_help(),
        ToolName::Simulate => print_simulate_help(),
        ToolName::ML => print_ml_help(),
        ToolName::Chaos => print_chaos_help(),
        ToolName::Units => print_units_help(),
        ToolName::Validate => print_validate_help(),
    }
}

fn print_solve_help() {
    println!("TOOL: solve");
    println!("DESCRIPTION: Solve equations, systems, and mathematical problems");
    println!();
    println!("EQUATION TYPES:");
    println!();

    println!("  Einstein Equations:");
    println!("    - vacuum              - Vacuum Einstein equations (Rμν = 0)");
    println!("    - with_source         - Einstein equations with energy-momentum tensor");
    println!("    - schwarzschild       - Schwarzschild black hole solution");
    println!("    - kerr_newman         - Kerr-Newman rotating charged black hole");
    println!();

    println!("  Fluid Equations:");
    println!("    - navier_stokes       - Navier-Stokes equations for viscous flow");
    println!("    - euler               - Inviscid Euler equations");
    println!();

    println!("  Differential Equations:");
    println!("    - ode                 - Ordinary differential equations");
    println!("    - pde                 - Partial differential equations");
    println!("    - boundary_value      - Boundary value problems");
    println!();

    println!("  Other Types:");
    println!("    - linear_system       - Systems of linear equations (Ax = b)");
    println!("    - root_finding        - Finding roots of equations");
    println!("    - optimization        - Optimization problems");
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

fn print_compute_help() {
    println!("TOOL: compute");
    println!("DESCRIPTION: Compute numerical operations, matrix decompositions, special functions");
    println!();
    println!("OPERATION CATEGORIES:");
    println!();

    println!("  Matrix Operations:");
    println!("    - determinant, inverse, eigenvalues, eigenvectors");
    println!("    - svd, qr, lu, cholesky decompositions");
    println!("    - norm, rank, trace, exp");
    println!();

    println!("  Calculus:");
    println!("    - differentiate (symbolic, numeric, partial)");
    println!("    - integrate (definite, indefinite, numeric, multiple)");
    println!();

    println!("  Transforms:");
    println!("    - fourier, laplace, z, wavelet, hilbert");
    println!();

    println!("  Field Theory:");
    println!("    - electric, magnetic, gravitational fields");
    println!("    - green_function, decoherence_scale, bohm_potential");
    println!();

    println!("  Sampling:");
    println!("    - monte_carlo, mcmc, importance, rejection");
    println!();

    println!("  Special Functions:");
    println!("    - bessel, gamma, erf, elliptic, airy, legendre");
    println!();

    println!("  Physics (60+ operations):");
    println!("    - Relativity, Quantum, Statistical physics, Control theory");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "compute",
    "input": {{
      "operation": {{"matrix": "determinant"}},
      "data": {{"matrix": [[1, 2], [3, 4]]}}
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
    println!("    - simplify            - Algebraic simplification");
    println!("    - expand              - Expand expressions");
    println!("    - factor              - Factor expressions");
    println!();

    println!("  Series and Limits:");
    println!("    - series              - Taylor/Laurent series");
    println!("    - limit               - Limit computation");
    println!("    - asymptotic          - Asymptotic analysis");
    println!();

    println!("  Analysis:");
    println!("    - roots               - Find roots");
    println!("    - extrema             - Find extrema");
    println!("    - stability           - Stability analysis");
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

    println!("  Time Evolution:");
    println!("    - euler, runge_kutta4, adaptive_step");
    println!();

    println!("  Stochastic Processes:");
    println!("    - brownian_motion, geometric_brownian, ornstein_uhlenbeck");
    println!("    - jump_diffusion, levy process");
    println!();

    println!("  Finance Models:");
    println!("    - black_scholes, heston, sabr");
    println!();

    println!("  Fluid Dynamics:");
    println!("    - navier_stokes_2d, lattice_boltzmann");
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
      "steps": 1000
    }}
  }}"#);
}

fn print_ml_help() {
    println!("TOOL: ml");
    println!("DESCRIPTION: Machine learning operations");
    println!();
    println!("OPERATIONS:");
    println!();

    println!("  Clustering:");
    println!("    - kmeans, dbscan, hierarchical, gaussian_mixture");
    println!();

    println!("  Regression:");
    println!("    - linear, logistic, ridge, lasso");
    println!();

    println!("  Neural Networks:");
    println!("    - create, train, predict, forward, backward");
    println!();

    println!("  Dimensionality Reduction:");
    println!("    - pca, tsne, umap");
    println!();

    println!("  Classification:");
    println!("    - svm, decision_tree, random_forest");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "ml",
    "input": {{
      "operation": {{"clustering": "kmeans"}},
      "data": [[1.0, 2.0], [1.5, 1.8], [5.0, 8.0], [6.0, 9.0]],
      "parameters": {{"n_clusters": 2}}
    }}
  }}"#);
}

fn print_chaos_help() {
    println!("TOOL: chaos");
    println!("DESCRIPTION: Chaos theory operations");
    println!();
    println!("OPERATIONS:");
    println!();

    println!("  Fractals:");
    println!("    - mandelbrot, julia, burning_ship, koch, sierpinski");
    println!();

    println!("  Attractors:");
    println!("    - lorenz, rossler, henon, logistic");
    println!();

    println!("  Lyapunov:");
    println!("    - exponent, spectrum");
    println!();

    println!("  Bifurcation:");
    println!("    - logistic, period_doubling");
    println!();

    println!("  Dimension:");
    println!("    - box_counting, correlation, kaplan_yorke");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "chaos",
    "input": {{
      "operation": {{"fractal": "mandelbrot"}},
      "parameters": {{"max_iter": 1000}},
      "resolution": [800, 600]
    }}
  }}"#);
}

fn print_units_help() {
    println!("TOOL: units");
    println!("DESCRIPTION: Dimensional analysis and unit conversion");
    println!();
    println!("OPERATIONS:");
    println!();

    println!("  Conversion:");
    println!("    - convert             - Convert between units");
    println!();

    println!("  Analysis:");
    println!("    - analyze             - Dimensional analysis");
    println!("    - check_compatibility - Check if units are compatible");
    println!("    - get_base            - Get SI base units");
    println!("    - derive              - Derive units for quantity");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "units",
    "input": {{
      "operation": "convert",
      "value": 100.0,
      "from_unit": "m",
      "to_unit": "ft"
    }}
  }}"#);
}

fn print_validate_help() {
    println!("TOOL: validate");
    println!("DESCRIPTION: Equation and physics validation");
    println!();
    println!("OPERATIONS:");
    println!();

    println!("  Validation:");
    println!("    - equation            - Validate equation syntax/structure");
    println!("    - dimensions          - Check dimensional consistency");
    println!("    - conservation        - Verify conservation laws");
    println!("    - symmetry            - Analyze symmetries");
    println!("    - physics             - Check physics compliance");
    println!("    - bounds              - Check bounds and ranges");
    println!("    - singularities       - Detect singularities");
    println!();

    println!("JSON EXAMPLE:");
    println!(r#"  {{
    "tool": "validate",
    "input": {{
      "operation": "equation",
      "expression": "E = m*c^2",
      "variable_units": {{
        "E": "J",
        "m": "kg",
        "c": "m/s"
      }}
    }}
  }}"#);
}

/// Count total operations across all tools
pub fn count_operations() -> usize {
    // Approximation - covers all 8 tools
    229
}
