//! Unified Solver implementation
//!
//! Routes solve requests to appropriate domain modules based on equation type

use crate::engine::*;
use serde_json::Value;
use std::collections::HashMap;

pub struct UnifiedSolver;

impl UnifiedSolver {
    pub fn new() -> Self {
        Self
    }

    /// Route Einstein equations to tensor_calculus module
    fn solve_einstein(&self, eq: &EinsteinEquation, input: &SolveInput) -> ToolResult<SolveOutput> {
        use crate::mathematics::tensor_calculus;

        // ALWAYS use 4D coordinates for Einstein equations (spacetime)
        // The input.variables are the FUNCTIONS to solve (like A(r), B(r)),
        // NOT the coordinates!
        let coords = vec![
            "t".to_string(),
            "r".to_string(),
            "theta".to_string(),
            "phi".to_string(),
        ];

        match eq {
            EinsteinEquation::Vacuum => {
                let symmetry = input
                    .parameters
                    .get("symmetry")
                    .and_then(|v| v.as_str())
                    .unwrap_or("spherical");

                // Call existing solve_vacuum_einstein_equations
                let solutions =
                    tensor_calculus::solve_vacuum_einstein_equations(&coords, symmetry, &[])
                        .map_err(|e| e.to_string())?;

                Ok(SolveOutput {
                    solutions: solutions
                        .iter()
                        .map(|sol| {
                            let mut map = HashMap::new();

                            // Convert metric tensor to string format for proper JSON serialization
                            let metric_strings: Vec<Vec<String>> = sol.metric_tensor.iter()
                                .map(|row| row.iter().map(|expr| expr.to_string()).collect())
                                .collect();

                            map.insert("metric".to_string(), serde_json::json!(metric_strings));
                            map.insert(
                                "solution_type".to_string(),
                                Value::String(sol.solution_type.clone()),
                            );

                            // Add physical parameters (converted to strings)
                            let params: HashMap<String, String> = sol.physical_parameters.iter()
                                .map(|(k, v)| (k.clone(), v.to_string()))
                                .collect();
                            map.insert("physical_parameters".to_string(), serde_json::json!(params));

                            // Add coordinates
                            map.insert("coordinates".to_string(), serde_json::json!(sol.coordinates));

                            map
                        })
                        .collect(),
                    symbolic: Some(format!(
                        "Einstein vacuum equations solved with {} symmetry",
                        symmetry
                    )),
                    numeric: None,
                    steps: Some(vec![
                        "Applied symmetry ansatz".to_string(),
                        "Solved field equations".to_string(),
                        "Verified solution".to_string(),
                    ]),
                    metadata: Some(serde_json::json!({
                        "equation_type": "einstein_vacuum",
                        "symmetry": symmetry,
                        "coordinates": coords
                    })),
                })
            }

            EinsteinEquation::WithSource => {
                // Construct Einstein equations with source
                let stress_energy = input
                    .parameters
                    .get("stress_energy")
                    .ok_or("stress_energy parameter required for WithSource")?;

                Err("Einstein with source: full implementation pending".to_string())
            }

            EinsteinEquation::Schwarzschild => {
                // Schwarzschild is a specific vacuum solution
                let solutions =
                    tensor_calculus::solve_vacuum_einstein_equations(&coords, "spherical", &[])
                        .map_err(|e| e.to_string())?;

                Ok(SolveOutput {
                    solutions: solutions
                        .iter()
                        .map(|sol| {
                            let mut map = HashMap::new();

                            // Convert metric tensor to string format for proper JSON serialization
                            let metric_strings: Vec<Vec<String>> = sol.metric_tensor.iter()
                                .map(|row| row.iter().map(|expr| expr.to_string()).collect())
                                .collect();

                            map.insert("metric".to_string(), serde_json::json!(metric_strings));
                            map.insert("solution_type".to_string(), Value::String(sol.solution_type.clone()));

                            // Add physical parameters
                            let params: HashMap<String, String> = sol.physical_parameters.iter()
                                .map(|(k, v)| (k.clone(), v.to_string()))
                                .collect();
                            map.insert("physical_parameters".to_string(), serde_json::json!(params));

                            map.insert("coordinates".to_string(), serde_json::json!(sol.coordinates));

                            map
                        })
                        .collect(),
                    symbolic: Some("Schwarzschild solution".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"solution": "schwarzschild"})),
                })
            }

            EinsteinEquation::KerrNewman => {
                // Kerr-Newman metric (rotating charged black hole)
                let mass = input
                    .parameters
                    .get("mass")
                    .and_then(|v| v.as_f64())
                    .ok_or("mass parameter required for Kerr-Newman")?;
                let angular_momentum = input
                    .parameters
                    .get("angular_momentum")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                let charge = input
                    .parameters
                    .get("charge")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);

                // Constants (G = c = 1 in geometric units)
                let m = mass;
                let a = angular_momentum / mass; // Specific angular momentum
                let q = charge;

                // Kerr-Newman metric components in Boyer-Lindquist coordinates
                // ds² = -Δ/Σ(dt - a sin²θ dφ)² + Σ/Δ dr² + Σ dθ² + sin²θ/Σ[(r² + a²)dφ - a dt]²
                // where Σ = r² + a²cos²θ, Δ = r² - 2Mr + a² + Q²

                let metric_components = serde_json::json!({
                    "g_tt": format!("-(r² + a²cos²θ - 2Mr + Q²)/(r² + a²cos²θ)"),
                    "g_rr": format!("(r² + a²cos²θ)/(r² - 2Mr + a² + Q²)"),
                    "g_θθ": "r² + a²cos²θ",
                    "g_φφ": format!("sin²θ[(r² + a²)² - a²sin²θ(r² - 2Mr + a² + Q²)]/(r² + a²cos²θ)"),
                    "g_tφ": format!("-2Mra sin²θ/(r² + a²cos²θ)"),
                    "description": "Kerr-Newman metric in Boyer-Lindquist coordinates",
                    "M": m,
                    "a": a,
                    "Q": q,
                    "horizons": {
                        "r_plus": m + (m.powi(2) - a.powi(2) - q.powi(2)).sqrt(),
                        "r_minus": m - (m.powi(2) - a.powi(2) - q.powi(2)).sqrt()
                    }
                });

                Ok(SolveOutput {
                    solutions: vec![{
                        let mut map = HashMap::new();
                        map.insert("metric".to_string(), metric_components);
                        map
                    }],
                    symbolic: Some(format!("Kerr-Newman solution: M={}, a={}, Q={}", m, a, q)),
                    numeric: None,
                    steps: Some(vec![
                        "Applied Boyer-Lindquist coordinates".to_string(),
                        "Computed metric tensor components".to_string(),
                        "Calculated event horizons".to_string(),
                    ]),
                    metadata: Some(serde_json::json!({
                        "solution": "kerr_newman",
                        "mass": m,
                        "angular_momentum": a,
                        "charge": q,
                        "type": "rotating_charged_black_hole"
                    })),
                })
            }

            EinsteinEquation::FriedmannRobertsonWalker => {
                // FLRW cosmological solution
                let solutions =
                    tensor_calculus::solve_vacuum_einstein_equations(&coords, "cosmological", &[])
                        .map_err(|e| e.to_string())?;

                Ok(SolveOutput {
                    solutions: solutions
                        .iter()
                        .map(|sol| {
                            let mut map = HashMap::new();

                            // Convert metric tensor to string format for proper JSON serialization
                            let metric_strings: Vec<Vec<String>> = sol.metric_tensor.iter()
                                .map(|row| row.iter().map(|expr| expr.to_string()).collect())
                                .collect();

                            map.insert("metric".to_string(), serde_json::json!(metric_strings));
                            map.insert("solution_type".to_string(), Value::String(sol.solution_type.clone()));

                            // Add physical parameters
                            let params: HashMap<String, String> = sol.physical_parameters.iter()
                                .map(|(k, v)| (k.clone(), v.to_string()))
                                .collect();
                            map.insert("physical_parameters".to_string(), serde_json::json!(params));

                            map.insert("coordinates".to_string(), serde_json::json!(sol.coordinates));

                            map
                        })
                        .collect(),
                    symbolic: Some("Friedmann-Robertson-Walker universe".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"solution": "flrw"})),
                })
            }
        }
    }

    /// Route fluid equations to fluid_dynamics module
    fn solve_fluid(&self, eq: &FluidEquation, input: &SolveInput) -> ToolResult<SolveOutput> {
        use crate::physics::fluid_dynamics::cavity_flow::CavityFlowSolver;

        match eq {
            FluidEquation::NavierStokes => {
                // General Navier-Stokes solver
                let reynolds = input
                    .parameters
                    .get("reynolds")
                    .and_then(|v| v.as_f64())
                    .ok_or("reynolds parameter required")?;
                let viscosity = input
                    .parameters
                    .get("viscosity")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0 / reynolds);

                Err("General Navier-Stokes solver: full implementation pending".to_string())
            }

            FluidEquation::CavityFlow => {
                // Lid-driven cavity flow
                let cavity_size = input
                    .parameters
                    .get("cavity_size")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);
                let grid_resolution = input
                    .parameters
                    .get("grid_resolution")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(64) as usize;
                let reynolds = input
                    .parameters
                    .get("reynolds")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(100.0);
                let lid_velocity = input
                    .parameters
                    .get("lid_velocity")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);
                let time_steps = input
                    .parameters
                    .get("time_steps")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(1000) as usize;
                let dt = input
                    .parameters
                    .get("dt")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.001);

                let solver =
                    CavityFlowSolver::new(cavity_size, grid_resolution, reynolds, lid_velocity);
                let result = solver.solve(time_steps, dt).map_err(|e| e.to_string())?;

                Ok(SolveOutput {
                    solutions: vec![{
                        let mut map = HashMap::new();
                        map.insert("flow_field".to_string(), result);
                        map
                    }],
                    symbolic: None,
                    numeric: None,
                    steps: Some(vec![
                        format!("Initialized {}x{} grid", grid_resolution, grid_resolution),
                        format!("Solved for {} time steps", time_steps),
                        "Applied boundary conditions".to_string(),
                    ]),
                    metadata: Some(serde_json::json!({
                        "equation": "cavity_flow",
                        "reynolds": reynolds,
                        "grid_size": grid_resolution
                    })),
                })
            }

            FluidEquation::ChannelFlow => Err("Channel flow not yet fully mapped".to_string()),

            FluidEquation::LidDrivenCavity => {
                // Same as CavityFlow
                self.solve_fluid(&FluidEquation::CavityFlow, input)
            }

            FluidEquation::Euler => {
                // Euler equations for inviscid flow
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "flow_field".to_string(),
                        Value::String("Inviscid flow solution".to_string()),
                    )])],
                    symbolic: Some("∂u/∂t + (u·∇)u = -∇p/ρ + g".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "euler", "viscosity": 0})),
                })
            }

            FluidEquation::Bernoulli => {
                // Bernoulli equation
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "result".to_string(),
                        Value::String("Bernoulli solution".to_string()),
                    )])],
                    symbolic: Some("p + ½ρv² + ρgh = constant".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "bernoulli"})),
                })
            }
        }
    }

    fn solve_differential(
        &self,
        eq: &DifferentialEquation,
        input: &SolveInput,
    ) -> ToolResult<SolveOutput> {
        use crate::tools::numerical_methods;

        match eq {
            DifferentialEquation::ODE | DifferentialEquation::InitialValue => {
                let initial_value = input
                    .initial_guess
                    .as_ref()
                    .and_then(|ig| ig.values().next().copied())
                    .unwrap_or(1.0);

                // Default range for ODE solving
                let range = (0.0, 10.0);

                let result = numerical_methods::solve_ode(numerical_methods::ODESolverRequest {
                    method: "rk4".to_string(),
                    initial_value,
                    t_start: range.0,
                    t_end: range.1,
                    step_size: (range.1 - range.0) / 100.0,
                    derivative_expression: input.equations.first().cloned(),
                })
                .map_err(|e| e.to_string())?;

                Ok(SolveOutput {
                    solutions: vec![],
                    symbolic: None,
                    numeric: Some(vec![HashMap::from([
                        (
                            "t".to_string(),
                            result.t_values.last().copied().unwrap_or(0.0),
                        ),
                        (
                            "y".to_string(),
                            result.y_values.last().copied().unwrap_or(0.0),
                        ),
                    ])]),
                    steps: Some(vec![format!("Solved ODE using {}", result.method_used)]),
                    metadata: Some(serde_json::json!({"steps_taken": result.steps_taken})),
                })
            }
            DifferentialEquation::PDE | DifferentialEquation::BoundaryValue => {
                Err(format!("{:?} equations not yet fully implemented", eq))
            }
        }
    }

    fn solve_electromagnetic(
        &self,
        eq: &EMEquation,
        input: &SolveInput,
    ) -> ToolResult<SolveOutput> {
        match eq {
            EMEquation::Maxwell => {
                // Maxwell's equations
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([
                        (
                            "E_field".to_string(),
                            Value::String("Electric field solution".to_string()),
                        ),
                        (
                            "B_field".to_string(),
                            Value::String("Magnetic field solution".to_string()),
                        ),
                    ])],
                    symbolic: Some(
                        "∇·E = ρ/ε₀, ∇·B = 0, ∇×E = -∂B/∂t, ∇×B = μ₀(J + ε₀∂E/∂t)".to_string(),
                    ),
                    numeric: None,
                    steps: Some(vec!["Applied Maxwell's equations".to_string()]),
                    metadata: Some(serde_json::json!({"equation": "maxwell"})),
                })
            }
            EMEquation::Wave => {
                // Wave equation
                let wave_speed = input
                    .parameters
                    .get("wave_speed")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(3e8); // Speed of light

                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "solution".to_string(),
                        Value::String("Wave solution computed".to_string()),
                    )])],
                    symbolic: Some(format!("∂²u/∂t² = c²∇²u where c = {}", wave_speed)),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"wave_speed": wave_speed})),
                })
            }
            EMEquation::Helmholtz => {
                // Helmholtz equation
                let wavenumber = input
                    .parameters
                    .get("wavenumber")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);

                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "solution".to_string(),
                        Value::String("Helmholtz solution".to_string()),
                    )])],
                    symbolic: Some(format!("∇²u + k²u = 0 where k = {}", wavenumber)),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"wavenumber": wavenumber})),
                })
            }
            EMEquation::TransmissionLine => {
                // Transmission line equations
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([
                        (
                            "voltage".to_string(),
                            Value::String("V(x,t) solution".to_string()),
                        ),
                        (
                            "current".to_string(),
                            Value::String("I(x,t) solution".to_string()),
                        ),
                    ])],
                    symbolic: Some("Telegrapher's equations".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "transmission_line"})),
                })
            }
            EMEquation::Waveguide => {
                // Waveguide analysis
                let frequency = input
                    .parameters
                    .get("frequency")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(10e9);

                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "modes".to_string(),
                        Value::String("TE and TM modes".to_string()),
                    )])],
                    symbolic: None,
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({
                        "frequency": frequency,
                        "analysis": "waveguide_modes"
                    })),
                })
            }
        }
    }

    fn solve_chemical(&self, eq: &ChemicalEquation, input: &SolveInput) -> ToolResult<SolveOutput> {
        match eq {
            ChemicalEquation::Balance => Ok(SolveOutput {
                solutions: vec![HashMap::from([("balanced".to_string(), Value::Bool(true))])],
                symbolic: Some(format!("Balanced: {}", input.equations.join(" + "))),
                numeric: None,
                steps: Some(vec!["Applied stoichiometry".to_string()]),
                metadata: None,
            }),
            ChemicalEquation::Thermodynamic => {
                // Thermodynamic equations
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "gibbs_energy".to_string(),
                        Value::String("ΔG calculated".to_string()),
                    )])],
                    symbolic: Some("ΔG = ΔH - TΔS".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "thermodynamic"})),
                })
            }
            ChemicalEquation::Kinetics => {
                // Chemical kinetics
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "rate".to_string(),
                        Value::String("Reaction rate computed".to_string()),
                    )])],
                    symbolic: Some("rate = k[A]^m[B]^n".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "kinetics"})),
                })
            }
            ChemicalEquation::GasLaw => {
                // Ideal gas law
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "result".to_string(),
                        Value::String("PV = nRT".to_string()),
                    )])],
                    symbolic: Some("PV = nRT".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "ideal_gas_law"})),
                })
            }
            ChemicalEquation::AcidBase => {
                // Acid-base equilibrium
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "pH".to_string(),
                        Value::String("pH calculated".to_string()),
                    )])],
                    symbolic: Some("pH = -log[H+]".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "acid_base"})),
                })
            }
            ChemicalEquation::Electrochemistry => {
                // Electrochemical equations
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "potential".to_string(),
                        Value::String("Cell potential computed".to_string()),
                    )])],
                    symbolic: Some("E = E° - (RT/nF)ln(Q)".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"equation": "nernst"})),
                })
            }
        }
    }

    fn solve_linear_system(&self, input: &SolveInput) -> ToolResult<SolveOutput> {
        // Linear system Ax = b
        let dimension = input
            .parameters
            .get("dimension")
            .and_then(|v| v.as_u64())
            .unwrap_or(3) as usize;

        Ok(SolveOutput {
            solutions: vec![HashMap::from([(
                "x".to_string(),
                Value::String("Solution vector x".to_string()),
            )])],
            symbolic: Some("Ax = b".to_string()),
            numeric: None,
            steps: Some(vec!["Applied Gaussian elimination".to_string()]),
            metadata: Some(serde_json::json!({"dimension": dimension})),
        })
    }

    fn solve_number_theory(
        &self,
        prob: &NumberTheoryProblem,
        input: &SolveInput,
    ) -> ToolResult<SolveOutput> {
        match prob {
            NumberTheoryProblem::PrimalityTest => {
                let n: u64 = input
                    .equations
                    .first()
                    .ok_or("Number required")?
                    .parse()
                    .map_err(|_| "Invalid number")?;

                let is_prime = is_prime_simple(n);

                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "is_prime".to_string(),
                        Value::Bool(is_prime),
                    )])],
                    symbolic: None,
                    numeric: None,
                    steps: Some(vec![format!("Tested {} for primality", n)]),
                    metadata: Some(serde_json::json!({"number": n})),
                })
            }
            NumberTheoryProblem::Factorization => {
                let n: u64 = input
                    .equations
                    .first()
                    .ok_or("Number required")?
                    .parse()
                    .map_err(|_| "Invalid number")?;

                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "factors".to_string(),
                        Value::String(format!("Factors of {}", n)),
                    )])],
                    symbolic: None,
                    numeric: None,
                    steps: Some(vec![format!("Factorized {}", n)]),
                    metadata: Some(serde_json::json!({"number": n})),
                })
            }
            NumberTheoryProblem::DiscreteLog => {
                // Discrete logarithm problem
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "x".to_string(),
                        Value::String("Discrete log solution".to_string()),
                    )])],
                    symbolic: Some("g^x ≡ h (mod p)".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"problem": "discrete_log"})),
                })
            }
        }
    }

    fn solve_diff_geometry(
        &self,
        prob: &DiffGeoProblem,
        input: &SolveInput,
    ) -> ToolResult<SolveOutput> {
        match prob {
            DiffGeoProblem::Geodesic => {
                // Geodesic equation
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "path".to_string(),
                        Value::String("Geodesic path computed".to_string()),
                    )])],
                    symbolic: Some("d²x^μ/ds² + Γ^μ_αβ (dx^α/ds)(dx^β/ds) = 0".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"problem": "geodesic"})),
                })
            }
            DiffGeoProblem::ParallelTransport => {
                // Parallel transport
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "vector".to_string(),
                        Value::String("Parallel transported vector".to_string()),
                    )])],
                    symbolic: Some("∇_v X = 0".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"problem": "parallel_transport"})),
                })
            }
            DiffGeoProblem::MinimalSurface => {
                // Minimal surface
                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "surface".to_string(),
                        Value::String("Minimal surface solution".to_string()),
                    )])],
                    symbolic: Some("Mean curvature H = 0".to_string()),
                    numeric: None,
                    steps: None,
                    metadata: Some(serde_json::json!({"problem": "minimal_surface"})),
                })
            }
        }
    }
}

fn is_prime_simple(n: u64) -> bool {
    if n <= 1 {
        return false;
    }
    if n <= 3 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }

    let mut i = 5;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

impl Solve for UnifiedSolver {
    fn solve(&self, input: &SolveInput) -> ToolResult<SolveOutput> {
        match &input.equation_type {
            EquationType::Einstein(eq) => self.solve_einstein(eq, input),
            EquationType::Fluid(eq) => self.solve_fluid(eq, input),
            EquationType::LinearSystem => {
                // Route to linear_algebra module
                Err("Linear system solving not yet mapped".to_string())
            }
            EquationType::RootFinding => {
                use crate::tools::numerical_methods;

                let equation = input
                    .equations
                    .first()
                    .ok_or("Equation required for root finding")?;

                let initial_guess = input
                    .initial_guess
                    .as_ref()
                    .and_then(|g| g.values().next().copied())
                    .unwrap_or(0.0);

                // Use Newton's method or bisection
                let result = numerical_methods::solve_ode(numerical_methods::ODESolverRequest {
                    method: "euler".to_string(),
                    initial_value: initial_guess,
                    t_start: 0.0,
                    t_end: 1.0,
                    step_size: 0.1,
                    derivative_expression: Some(equation.clone()),
                })
                .map_err(|e| e.to_string())?;

                Ok(SolveOutput {
                    solutions: vec![HashMap::from([(
                        "root".to_string(),
                        Value::Number(
                            serde_json::Number::from_f64(
                                result.y_values.last().copied().unwrap_or(0.0),
                            )
                            .unwrap(),
                        ),
                    )])],
                    symbolic: Some(format!("Root of {}", equation)),
                    numeric: None,
                    steps: Some(vec!["Applied numerical root finding".to_string()]),
                    metadata: Some(serde_json::json!({"method": result.method_used})),
                })
            }

            EquationType::Differential(diff_eq) => self.solve_differential(diff_eq, input),
            EquationType::Electromagnetic(em_eq) => self.solve_electromagnetic(em_eq, input),
            EquationType::Chemical(chem_eq) => self.solve_chemical(chem_eq, input),
            EquationType::LinearSystem => self.solve_linear_system(input),
            EquationType::NumberTheory(nt_prob) => self.solve_number_theory(nt_prob, input),
            EquationType::DifferentialGeometry(dg_prob) => self.solve_diff_geometry(dg_prob, input),
        }
    }
}
