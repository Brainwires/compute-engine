//! FULLY AUTOMATIC help system that derives ALL documentation from enum definitions
//!
//! **THIS IS A SELF-DOCUMENTING SYSTEM:**
//! When you add a new enum variant, it automatically appears in help with ZERO manual updates.
//!
//! Uses strum::EnumIter to programmatically discover all enum variants at compile time.

use crate::engine::equations::*;
use strum::IntoEnumIterator;

/// Print hierarchical help for a tool and optional category
pub fn print_hierarchical_help(tool: &str, category: Option<&str>) {
    match tool.to_lowercase().as_str() {
        "compute" => print_compute_help(category),
        "solve" => print_solve_help(category),
        "differentiate" => print_differentiate_help(category),
        "integrate" => print_integrate_help(category),
        "analyze" => print_analyze_help(category),
        "simulate" => print_simulate_help(category),
        "transform" => print_transform_help(category),
        "fieldtheory" => print_fieldtheory_help(category),
        "sample" => print_sample_help(category),
        "optimize" => print_optimize_help(category),
        _ => {
            eprintln!("Unknown tool: {}", tool);
            eprintln!("Available tools: solve, differentiate, integrate, analyze, simulate,");
            eprintln!("                 compute, transform, fieldtheory, sample, optimize");
        }
    }
}

// ============================================================================
// COMPUTE TOOL - Fully automatic from enums
// ============================================================================

fn print_compute_help(category: Option<&str>) {
    println!("TOOL: compute");
    println!("DESCRIPTION: Compute tensor operations, matrix decompositions, special functions");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "tensor" => print_enum_ops::<TensorOp>("TENSOR OPERATIONS", "compute", "tensor"),
            "matrix" => print_enum_ops::<MatrixOp>("MATRIX OPERATIONS", "compute", "matrix"),
            "decomp" | "decomposition" => print_enum_ops::<MatrixDecomp>("MATRIX DECOMPOSITIONS", "compute", "matrix_decomp"),
            "special" | "functions" => print_enum_ops::<SpecialFunction>("SPECIAL FUNCTIONS", "compute", "special_func"),
            "number" | "numbertheory" => print_enum_ops::<NumberTheoryOp>("NUMBER THEORY & CRYPTOGRAPHY", "compute", "number_theory"),
            "geometry" => print_enum_ops::<GeometryOp>("COMPUTATIONAL GEOMETRY", "compute", "geometry"),
            "information" => print_enum_ops::<InformationOp>("INFORMATION THEORY", "compute", "information"),
            "chemistry" => print_enum_ops::<ChemistryOp>("CHEMISTRY", "compute", "chemistry"),
            "biology" => print_enum_ops::<BiologyOp>("BIOLOGY", "compute", "biology"),
            "thermodynamics" | "thermo" => print_enum_ops::<ThermodynamicsOp>("THERMODYNAMICS", "compute", "thermodynamics"),
            "optics" => print_enum_ops::<OpticsOp>("OPTICS", "compute", "optics"),
            "geophysics" | "geo" => print_enum_ops::<GeophysicsOp>("GEOPHYSICS", "compute", "geophysics"),
            "engineering" => print_enum_ops::<EngineeringOp>("ENGINEERING", "compute", "engineering"),
            "datetime" | "time" => print_enum_ops::<DateTimeOp>("DATE/TIME CALCULATIONS", "compute", "date_time"),
            "graph" => print_enum_ops::<GraphOp>("GRAPH THEORY", "compute", "graph"),
            "physics" => print_physics_subcategories(None),
            "physics_relativity" | "relativity" => print_enum_ops::<RelativityOp>("PHYSICS → RELATIVITY", "compute", "physics.relativity"),
            "physics_statistical" | "statistical" => print_enum_ops::<StatPhysicsOp>("PHYSICS → STATISTICAL PHYSICS", "compute", "physics.statistical_physics"),
            "physics_quantum" | "quantum" => print_enum_ops::<QuantumMechOp>("PHYSICS → QUANTUM MECHANICS", "compute", "physics.quantum_mechanics"),
            "physics_control" | "control" => print_enum_ops::<ControlSystemsOp>("PHYSICS → CONTROL SYSTEMS", "compute", "physics.control_systems"),
            "physics_nuclear" | "nuclear" => print_enum_ops::<NuclearOp>("PHYSICS → NUCLEAR PHYSICS", "compute", "physics.nuclear_physics"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                eprintln!();
                print_compute_categories();
            }
        }
    } else {
        print_compute_categories();
    }
}

fn print_compute_categories() {
    println!("CATEGORIES:");
    println!("  Use: brainwires-compute-engine help compute <CATEGORY>");
    println!();

    // Automatically count operations in each category
    let tensor_count = TensorOp::iter().count();
    let matrix_count = MatrixOp::iter().count();
    let decomp_count = MatrixDecomp::iter().count();
    let special_count = SpecialFunction::iter().count();
    let number_count = NumberTheoryOp::iter().count();
    let geometry_count = GeometryOp::iter().count();
    let info_count = InformationOp::iter().count();
    let chem_count = ChemistryOp::iter().count();
    let bio_count = BiologyOp::iter().count();
    let thermo_count = ThermodynamicsOp::iter().count();
    let optics_count = OpticsOp::iter().count();
    let geo_count = GeophysicsOp::iter().count();
    let eng_count = EngineeringOp::iter().count();
    let dt_count = DateTimeOp::iter().count();
    let graph_count = GraphOp::iter().count();
    let rel_count = RelativityOp::iter().count();
    let stat_count = StatPhysicsOp::iter().count();
    let quant_count = QuantumMechOp::iter().count();
    let control_count = ControlSystemsOp::iter().count();
    let nuclear_count = NuclearOp::iter().count();

    println!("  tensor           - Tensor operations ({} ops)", tensor_count);
    println!("  matrix           - Matrix operations ({} ops)", matrix_count);
    println!("  decomposition    - Matrix decompositions ({} ops)", decomp_count);
    println!("  special          - Special functions ({} ops)", special_count);
    println!("  numbertheory     - Number theory & cryptography ({} ops)", number_count);
    println!("  geometry         - Computational geometry ({} ops)", geometry_count);
    println!("  information      - Information theory ({} ops)", info_count);
    println!("  chemistry        - Chemistry formulas ({} ops)", chem_count);
    println!("  biology          - Biology formulas ({} ops)", bio_count);
    println!("  thermodynamics   - Thermodynamics ({} ops)", thermo_count);
    println!("  optics           - Optics ({} ops)", optics_count);
    println!("  geophysics       - Geophysics ({} ops)", geo_count);
    println!("  engineering      - Engineering ({} ops)", eng_count);
    println!("  datetime         - Date/time calculations ({} ops)", dt_count);
    println!("  graph            - Graph theory ({} ops)", graph_count);
    println!("  physics          - Physics ({} ops total)", rel_count + stat_count + quant_count + control_count + nuclear_count);
    println!();
    println!("Total: {} operations",
        tensor_count + matrix_count + decomp_count + special_count + number_count +
        geometry_count + info_count + chem_count + bio_count + thermo_count +
        optics_count + geo_count + eng_count + dt_count + graph_count +
        rel_count + stat_count + quant_count + control_count + nuclear_count
    );
}

fn print_physics_subcategories(_sub: Option<&str>) {
    println!("COMPUTE → PHYSICS");
    println!();
    println!("SUBCATEGORIES:");
    println!("  Use: brainwires-compute-engine help compute physics_<SUBCATEGORY>");
    println!();

    let rel_count = RelativityOp::iter().count();
    let stat_count = StatPhysicsOp::iter().count();
    let quant_count = QuantumMechOp::iter().count();
    let control_count = ControlSystemsOp::iter().count();
    let nuclear_count = NuclearOp::iter().count();

    println!("  physics_relativity   - Special & General Relativity ({} ops)", rel_count);
    println!("  physics_statistical  - Statistical Physics ({} ops)", stat_count);
    println!("  physics_quantum      - Quantum Mechanics ({} ops)", quant_count);
    println!("  physics_control      - Control Systems ({} ops)", control_count);
    println!("  physics_nuclear      - Nuclear Physics ({} ops)", nuclear_count);
    println!();
    println!("Total: {} operations", rel_count + stat_count + quant_count + control_count + nuclear_count);
}

/// Generic function that prints ALL operations from any enum using strum iteration
fn print_enum_ops<T: IntoEnumIterator + std::fmt::Debug>(title: &str, tool: &str, op_path: &str) {
    println!("COMPUTE → {}", title);
    println!();
    println!("OPERATIONS:");

    let ops: Vec<T> = T::iter().collect();
    for op in &ops {
        // Format the debug string to snake_case
        let debug_str = format!("{:?}", op);
        let snake_case = to_snake_case(&debug_str);
        println!("  • {}", snake_case);
    }

    println!();
    println!("Total: {} operations", ops.len());
    println!();
    println!("JSON EXAMPLE:");
    let first_op = ops.first().map(|o| to_snake_case(&format!("{:?}", o))).unwrap_or_default();
    println!(r#"  {{
    "tool": "{}",
    "input": {{
      "operation": {{"{}": "{}"}},
      "data": {{"example": "data"}}
    }}
  }}"#, tool, op_path, first_op);
}

/// Convert PascalCase to snake_case
fn to_snake_case(s: &str) -> String {
    let mut result = String::new();
    for (i, ch) in s.chars().enumerate() {
        if ch.is_uppercase() {
            if i > 0 {
                result.push('_');
            }
            result.push(ch.to_lowercase().next().unwrap());
        } else {
            result.push(ch);
        }
    }
    result
}

// ============================================================================
// OTHER TOOLS - Fully automatic
// ============================================================================

fn print_solve_help(category: Option<&str>) {
    println!("TOOL: solve");
    println!("DESCRIPTION: Solve equations, systems, and mathematical problems");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "einstein" => print_enum_ops::<EinsteinEquation>("EINSTEIN EQUATIONS", "solve", "einstein"),
            "fluid" => print_enum_ops::<FluidEquation>("FLUID EQUATIONS", "solve", "fluid"),
            "differential" => print_enum_ops::<DifferentialEquation>("DIFFERENTIAL EQUATIONS", "solve", "differential"),
            "em" | "electromagnetic" => print_enum_ops::<EMEquation>("ELECTROMAGNETIC EQUATIONS", "solve", "electromagnetic"),
            "chemical" => print_enum_ops::<ChemicalEquation>("CHEMICAL EQUATIONS", "solve", "chemical"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_solve_categories();
            }
        }
    } else {
        print_solve_categories();
    }
}

fn print_solve_categories() {
    println!("EQUATION TYPES:");
    println!("  Use: brainwires-compute-engine help solve <CATEGORY>");
    println!();

    let einstein_count = EinsteinEquation::iter().count();
    let fluid_count = FluidEquation::iter().count();
    let diff_count = DifferentialEquation::iter().count();
    let em_count = EMEquation::iter().count();
    let chem_count = ChemicalEquation::iter().count();

    println!("  einstein         - Einstein equations ({} types)", einstein_count);
    println!("  fluid            - Fluid equations ({} types)", fluid_count);
    println!("  differential     - Differential equations ({} types)", diff_count);
    println!("  electromagnetic  - Electromagnetic equations ({} types)", em_count);
    println!("  chemical         - Chemical equations ({} types)", chem_count);
    println!("  linear_system    - Systems of linear equations");
    println!("  root_finding     - Root finding algorithms");
    println!();
    println!("Total: {} equation types", einstein_count + fluid_count + diff_count + em_count + chem_count + 2);
}

fn print_differentiate_help(category: Option<&str>) {
    println!("TOOL: differentiate");
    println!("DESCRIPTION: Compute derivatives, gradients, and differential operators");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "vector" | "vectorcalc" => print_enum_ops::<VectorCalcOp>("VECTOR CALCULUS", "differentiate", "vector_calc"),
            "tensor" | "tensorcalc" => print_enum_ops::<TensorDiffOp>("TENSOR CALCULUS", "differentiate", "tensor_calc"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_differentiate_categories();
            }
        }
    } else {
        print_differentiate_categories();
    }
}

fn print_differentiate_categories() {
    println!("OPERATION TYPES:");
    println!("  Use: brainwires-compute-engine help differentiate <CATEGORY>");
    println!();

    let vector_count = VectorCalcOp::iter().count();
    let tensor_count = TensorDiffOp::iter().count();

    println!("  vector           - Vector calculus ({} ops)", vector_count);
    println!("  tensor           - Tensor calculus ({} ops)", tensor_count);
    println!("  variational      - Variational calculus");
    println!("  differential_forms - Differential forms");
    println!("  numeric          - Numeric differentiation");
    println!("  symbolic         - Symbolic differentiation");
    println!();
    println!("Total: {} operation types", vector_count + tensor_count + 4);
}

fn print_integrate_help(category: Option<&str>) {
    println!("TOOL: integrate");
    println!("DESCRIPTION: Compute integrals (definite, indefinite, line, surface, volume)");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "geometric" => print_enum_ops::<GeometricIntegral>("GEOMETRIC INTEGRALS", "integrate", "geometric"),
            "theorem" | "theorems" => print_enum_ops::<IntegralTheorem>("INTEGRAL THEOREMS", "integrate", "theorem"),
            "complex" => print_enum_ops::<ComplexIntegral>("COMPLEX ANALYSIS", "integrate", "complex_analysis"),
            "numeric" => print_enum_ops::<NumericIntegration>("NUMERIC INTEGRATION", "integrate", "numeric"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_integrate_categories();
            }
        }
    } else {
        print_integrate_categories();
    }
}

fn print_integrate_categories() {
    println!("INTEGRATION TYPES:");
    println!("  Use: brainwires-compute-engine help integrate <CATEGORY>");
    println!();

    let geom_count = GeometricIntegral::iter().count();
    let theorem_count = IntegralTheorem::iter().count();
    let complex_count = ComplexIntegral::iter().count();
    let numeric_count = NumericIntegration::iter().count();

    println!("  geometric        - Line/surface/volume integrals ({} types)", geom_count);
    println!("  theorem          - Integral theorems ({} types)", theorem_count);
    println!("  complex          - Complex analysis ({} types)", complex_count);
    println!("  numeric          - Numeric integration ({} types)", numeric_count);
    println!("  symbolic         - Symbolic integration");
    println!("  monte_carlo      - Monte Carlo integration");
    println!();
    println!("Total: {} integration types", geom_count + theorem_count + complex_count + numeric_count + 2);
}

fn print_analyze_help(_category: Option<&str>) {
    println!("TOOL: analyze");
    println!("DESCRIPTION: Analyze expressions (validate, simplify, parse, check)");
    println!();
    print_enum_ops::<AnalysisOp>("ANALYSIS OPERATIONS", "analyze", "operation");
}

fn print_simulate_help(category: Option<&str>) {
    println!("TOOL: simulate");
    println!("DESCRIPTION: Simulate time evolution and stochastic processes");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "stochastic" => print_enum_ops::<StochasticProcess>("STOCHASTIC PROCESSES", "simulate", "stochastic"),
            "finance" => print_enum_ops::<FinanceModel>("FINANCE MODELS", "simulate", "finance"),
            "fluid" => print_enum_ops::<FluidSim>("FLUID DYNAMICS", "simulate", "fluid_dynamics"),
            "time" | "evolution" => print_enum_ops::<TimeEvolutionMethod>("TIME EVOLUTION", "simulate", "time_evolution"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_simulate_categories();
            }
        }
    } else {
        print_simulate_categories();
    }
}

fn print_simulate_categories() {
    println!("SIMULATION MODELS:");
    println!("  Use: brainwires-compute-engine help simulate <CATEGORY>");
    println!();

    let stoch_count = StochasticProcess::iter().count();
    let finance_count = FinanceModel::iter().count();
    let fluid_count = FluidSim::iter().count();
    let time_count = TimeEvolutionMethod::iter().count();

    println!("  stochastic       - Stochastic processes ({} types)", stoch_count);
    println!("  finance          - Finance models ({} types)", finance_count);
    println!("  fluid            - Fluid dynamics ({} types)", fluid_count);
    println!("  evolution        - Time evolution methods ({} types)", time_count);
    println!();
    println!("Total: {} simulation types", stoch_count + finance_count + fluid_count + time_count);
}

fn print_transform_help(category: Option<&str>) {
    println!("TOOL: transform");
    println!("DESCRIPTION: Apply transforms (Fourier, Laplace, wavelet, filters)");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "fourier" => print_enum_ops::<FourierTransform>("FOURIER TRANSFORMS", "transform", "fourier"),
            "laplace" => print_enum_ops::<LaplaceTransform>("LAPLACE TRANSFORMS", "transform", "laplace"),
            "wavelet" => print_enum_ops::<WaveletType>("WAVELET TRANSFORMS", "transform", "wavelet"),
            "fft" => print_enum_ops::<FFTType>("FFT", "transform", "fft"),
            "filter" => print_enum_ops::<FilterType>("FILTERS", "transform", "filter"),
            "window" => print_enum_ops::<WindowType>("WINDOW FUNCTIONS", "transform", "window"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_transform_categories();
            }
        }
    } else {
        print_transform_categories();
    }
}

fn print_transform_categories() {
    println!("TRANSFORM TYPES:");
    println!("  Use: brainwires-compute-engine help transform <CATEGORY>");
    println!();

    let fourier_count = FourierTransform::iter().count();
    let laplace_count = LaplaceTransform::iter().count();
    let wavelet_count = WaveletType::iter().count();
    let fft_count = FFTType::iter().count();
    let filter_count = FilterType::iter().count();
    let window_count = WindowType::iter().count();

    println!("  fourier          - Fourier transforms ({} types)", fourier_count);
    println!("  laplace          - Laplace transforms ({} types)", laplace_count);
    println!("  wavelet          - Wavelet transforms ({} types)", wavelet_count);
    println!("  fft              - Fast Fourier Transform ({} types)", fft_count);
    println!("  filter           - Digital filters ({} types)", filter_count);
    println!("  window           - Window functions ({} types)", window_count);
    println!("  conformal        - Conformal mappings");
    println!();
    println!("Total: {} transform types", fourier_count + laplace_count + wavelet_count + fft_count + filter_count + window_count + 1);
}

fn print_fieldtheory_help(category: Option<&str>) {
    println!("TOOL: fieldtheory");
    println!("DESCRIPTION: Work with physics fields (EM, gravity, quantum)");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "em" | "electromagnetic" => print_enum_ops::<EMField>("ELECTROMAGNETIC FIELDS", "fieldtheory", "em"),
            "quantum" => print_enum_ops::<QuantumFieldType>("QUANTUM FIELD THEORY", "fieldtheory", "quantum_field"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_fieldtheory_categories();
            }
        }
    } else {
        print_fieldtheory_categories();
    }
}

fn print_fieldtheory_categories() {
    println!("FIELD TYPES:");
    println!("  Use: brainwires-compute-engine help fieldtheory <CATEGORY>");
    println!();

    let em_count = EMField::iter().count();
    let quantum_count = QuantumFieldType::iter().count();

    println!("  em               - Electromagnetic fields ({} types)", em_count);
    println!("  quantum          - Quantum field theory ({} types)", quantum_count);
    println!("  green_function   - Green's functions");
    println!();
    println!("Total: {} field types", em_count + quantum_count + 1);
}

fn print_sample_help(category: Option<&str>) {
    println!("TOOL: sample");
    println!("DESCRIPTION: Sample and analyze statistical data");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "montecarlo" | "monte_carlo" => print_enum_ops::<MonteCarloMethod>("MONTE CARLO METHODS", "sample", "monte_carlo"),
            "stats" | "statistical" => print_enum_ops::<StatisticalMethod>("STATISTICAL METHODS", "sample", "stats"),
            "signal" => print_enum_ops::<SignalMethod>("SIGNAL ANALYSIS", "sample", "signal_analysis"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_sample_categories();
            }
        }
    } else {
        print_sample_categories();
    }
}

fn print_sample_categories() {
    println!("SAMPLING METHODS:");
    println!("  Use: brainwires-compute-engine help sample <CATEGORY>");
    println!();

    let mc_count = MonteCarloMethod::iter().count();
    let stats_count = StatisticalMethod::iter().count();
    let signal_count = SignalMethod::iter().count();

    println!("  montecarlo       - Monte Carlo methods ({} types)", mc_count);
    println!("  statistical      - Statistical methods ({} types)", stats_count);
    println!("  signal           - Signal analysis ({} types)", signal_count);
    println!("  path_generation  - Path generation");
    println!("  moments          - Moment calculation");
    println!();
    println!("Total: {} sampling methods", mc_count + stats_count + signal_count + 2);
}

fn print_optimize_help(category: Option<&str>) {
    println!("TOOL: optimize");
    println!("DESCRIPTION: Optimize and fit data to models");
    println!();

    if let Some(cat) = category {
        match cat.to_lowercase().as_str() {
            "fit" | "fitting" => print_enum_ops::<FitMethod>("CURVE FITTING METHODS", "optimize", "fit"),
            "minimize" | "minimization" => print_enum_ops::<MinimizationMethod>("MINIMIZATION METHODS", "optimize", "minimize"),
            "interpolation" => print_enum_ops::<InterpolationMethod>("INTERPOLATION METHODS", "optimize", "interpolation"),
            "dimensional" => print_enum_ops::<DimAnalysisMethod>("DIMENSIONAL ANALYSIS", "optimize", "dimensional_analysis"),
            _ => {
                eprintln!("Unknown category: {}", cat);
                print_optimize_categories();
            }
        }
    } else {
        print_optimize_categories();
    }
}

fn print_optimize_categories() {
    println!("OPTIMIZATION METHODS:");
    println!("  Use: brainwires-compute-engine help optimize <CATEGORY>");
    println!();

    let fit_count = FitMethod::iter().count();
    let min_count = MinimizationMethod::iter().count();
    let interp_count = InterpolationMethod::iter().count();
    let dim_count = DimAnalysisMethod::iter().count();

    println!("  fit              - Curve fitting ({} types)", fit_count);
    println!("  minimize         - Minimization ({} types)", min_count);
    println!("  interpolation    - Interpolation ({} types)", interp_count);
    println!("  dimensional      - Dimensional analysis ({} types)", dim_count);
    println!("  symbolic_regression - Symbolic regression");
    println!("  auto             - Automatic model selection");
    println!();
    println!("Total: {} optimization methods", fit_count + min_count + interp_count + dim_count + 2);
}
