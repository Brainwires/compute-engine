//! Comprehensive fluid dynamics integration tests
//!
//! Tests for all SIMULATE tool fluid dynamics operations including:
//! - Navier-Stokes solvers (2D, 3D)
//! - Quantum Navier-Stokes (1D, 2D)
//! - Lattice Boltzmann method
//! - Analytical solutions
//! - Flow analysis utilities

use computational_engine::create_default_dispatcher;
use computational_engine::engine::equations::*;
use computational_engine::engine::*;
use computational_engine::simulate::fluids::*;
use std::collections::HashMap;

// Helper function to create SimulateInput for fluid tests
fn create_fluid_input(
    model: FluidSim,
    range: Option<[f64; 2]>,
    steps: Option<usize>,
    parameters: HashMap<String, f64>,
) -> SimulateInput {
    SimulateInput {
        model: SimulationModel::FluidDynamics(model),
        equations: Vec::new(),
        variables: vec!["velocity_x".to_string(), "velocity_y".to_string()],
        parameters,
        initial_conditions: None,
        range,
        steps,
        method: None,
        num_paths: None,
    }
}

// ============================================================================
// LATTICE BOLTZMANN METHOD TESTS
// ============================================================================

#[test]
fn test_lattice_boltzmann_basic() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 20.0);
    params.insert("ny".to_string(), 20.0);
    params.insert("tau".to_string(), 0.6);
    params.insert("lid_velocity".to_string(), 0.1);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::LatticeBotzmann,
        None,
        Some(100),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Lattice Boltzmann: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("velocity_x"), "Should have velocity_x");
        assert!(output.results.contains_key("velocity_y"), "Should have velocity_y");
        assert!(output.results.contains_key("density"), "Should have density");

        let vel_x = &output.results["velocity_x"];
        assert_eq!(vel_x.len(), 400, "Should have 20x20 = 400 grid points");

        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "lattice_boltzmann_d2q9");
    }
}

#[test]
fn test_lattice_boltzmann_lid_driven_cavity() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 32.0);
    params.insert("ny".to_string(), 32.0);
    params.insert("tau".to_string(), 0.7);
    params.insert("lid_velocity".to_string(), 0.05);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::LatticeBotzmann,
        None,
        Some(500),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let vel_x = &output.results["velocity_x"];
        let vel_y = &output.results["velocity_y"];
        let density = &output.results["density"];

        // All values should be finite
        assert!(vel_x.iter().all(|&v| v.is_finite()), "All velocities should be finite");
        assert!(vel_y.iter().all(|&v| v.is_finite()));
        assert!(density.iter().all(|&d| d.is_finite() && d > 0.0), "Density should be positive");
    }
}

#[test]
fn test_lattice_boltzmann_different_grid_sizes() {
    let dispatcher = create_default_dispatcher();

    // Test with larger grid
    let mut params = HashMap::new();
    params.insert("nx".to_string(), 50.0);
    params.insert("ny".to_string(), 30.0);
    params.insert("tau".to_string(), 0.6);
    params.insert("lid_velocity".to_string(), 0.1);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::LatticeBotzmann,
        None,
        Some(50),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let vel_x = &output.results["velocity_x"];
        assert_eq!(vel_x.len(), 50 * 30, "Should have 50x30 = 1500 grid points");
    }
}

#[test]
fn test_lattice_boltzmann_metadata() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 16.0);
    params.insert("ny".to_string(), 16.0);
    params.insert("tau".to_string(), 0.55);
    params.insert("lid_velocity".to_string(), 0.08);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::LatticeBotzmann,
        None,
        Some(10),
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok());

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.expect("Should have metadata");
        assert_eq!(metadata["method"], "lattice_boltzmann_d2q9");
        assert_eq!(metadata["grid_size"][0], 16);
        assert_eq!(metadata["grid_size"][1], 16);
        assert_eq!(metadata["steps"], 10);
        assert_eq!(metadata["tau"], 0.55);
        assert_eq!(metadata["lid_velocity"], 0.08);
    }
}

// ============================================================================
// QUANTUM NAVIER-STOKES 1D TESTS
// ============================================================================

#[test]
fn test_qns_1d_basic() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 100.0);
    params.insert("domain_length".to_string(), 1.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("enable_quantum".to_string(), 1.0);
    params.insert("cfl".to_string(), 0.3);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes1D,
        Some([0.0, 0.1]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate QNS 1D: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("density"), "Should have density");
        assert!(output.results.contains_key("velocity"), "Should have velocity");
        assert!(output.results.contains_key("pressure"), "Should have pressure");
        assert!(output.results.contains_key("x"), "Should have x coordinates");

        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "quantum_navier_stokes_1d");
        assert_eq!(metadata["enable_quantum"], true);
    }
}

#[test]
fn test_qns_1d_sod_shock_tube() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 200.0);
    params.insert("domain_length".to_string(), 1.0);
    params.insert("viscosity".to_string(), 0.001);
    params.insert("enable_quantum".to_string(), 1.0);
    params.insert("cfl".to_string(), 0.5);
    params.insert("init_type".to_string(), 1.0); // Sod shock tube

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes1D,
        Some([0.0, 0.2]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Sod shock tube: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let density = &output.results["density"];
        let velocity = &output.results["velocity"];
        let pressure = &output.results["pressure"];

        // Check that we have the right number of points
        assert_eq!(density.len(), 200);
        assert_eq!(velocity.len(), 200);
        assert_eq!(pressure.len(), 200);

        // All values should be finite
        assert!(density.iter().all(|&d| d.is_finite() && d > 0.0));
        assert!(velocity.iter().all(|&v| v.is_finite()));
        assert!(pressure.iter().all(|&p| p.is_finite() && p > 0.0));
    }
}

#[test]
fn test_qns_1d_without_quantum() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 100.0);
    params.insert("domain_length".to_string(), 1.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("enable_quantum".to_string(), 0.0); // Disable quantum effects
    params.insert("cfl".to_string(), 0.5);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes1D,
        Some([0.0, 0.05]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should work without quantum effects");

    if let Ok(ToolResponse::Simulate(output)) = result {
        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["enable_quantum"], false);
    }
}

// ============================================================================
// QUANTUM NAVIER-STOKES 2D TESTS
// ============================================================================

#[test]
fn test_qns_2d_basic() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 32.0);
    params.insert("ny".to_string(), 32.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("enable_quantum".to_string(), 1.0);
    params.insert("cfl".to_string(), 0.3);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes2D,
        Some([0.0, 0.1]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate QNS 2D: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("density"), "Should have density");
        assert!(output.results.contains_key("velocity_x"), "Should have velocity_x");
        assert!(output.results.contains_key("velocity_y"), "Should have velocity_y");
        assert!(output.results.contains_key("pressure"), "Should have pressure");

        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "quantum_navier_stokes_2d");
    }
}

#[test]
fn test_qns_2d_taylor_green() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 64.0);
    params.insert("ny".to_string(), 64.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("enable_quantum".to_string(), 1.0);
    params.insert("cfl".to_string(), 0.3);
    params.insert("init_type".to_string(), 1.0); // Taylor-Green vortex
    params.insert("u0".to_string(), 1.0);
    params.insert("wavenumber".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes2D,
        Some([0.0, 0.5]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate Taylor-Green vortex: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        // Taylor-Green should have non-zero moments
        let moments = output.moments.expect("Should have moments");
        assert!(moments.contains_key("kinetic_energy"));
        assert!(moments.contains_key("enstrophy"));
        assert!(moments.contains_key("max_vorticity"));
        assert!(moments.contains_key("max_divergence"));

        // Kinetic energy should be positive for non-trivial flow
        let ke = moments["kinetic_energy"];
        assert!(ke >= 0.0, "Kinetic energy should be non-negative");
    }
}

#[test]
fn test_qns_2d_decaying_turbulence() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 32.0);
    params.insert("ny".to_string(), 32.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("enable_quantum".to_string(), 0.0);
    params.insert("cfl".to_string(), 0.3);
    params.insert("init_type".to_string(), 2.0); // Decaying turbulence
    params.insert("energy".to_string(), 0.5);
    params.insert("k_peak".to_string(), 4.0);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes2D,
        Some([0.0, 0.2]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate decaying turbulence: {:?}", result);
}

// ============================================================================
// NAVIER-STOKES 3D TESTS
// ============================================================================

#[test]
fn test_ns_3d_basic() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 16.0);
    params.insert("ny".to_string(), 16.0);
    params.insert("nz".to_string(), 16.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("cfl".to_string(), 0.5);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::NavierStokes3D,
        Some([0.0, 0.1]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate NS 3D: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        assert!(output.results.contains_key("velocity_x"), "Should have velocity_x");
        assert!(output.results.contains_key("velocity_y"), "Should have velocity_y");
        assert!(output.results.contains_key("velocity_z"), "Should have velocity_z");
        assert!(output.results.contains_key("pressure"), "Should have pressure");

        let metadata = output.metadata.unwrap();
        assert_eq!(metadata["method"], "navier_stokes_3d");
        assert_eq!(metadata["grid_size"][0], 16);
        assert_eq!(metadata["grid_size"][1], 16);
        assert_eq!(metadata["grid_size"][2], 16);
    }
}

#[test]
fn test_ns_3d_taylor_green() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 16.0);
    params.insert("ny".to_string(), 16.0);
    params.insert("nz".to_string(), 16.0);
    params.insert("viscosity".to_string(), 0.01);
    params.insert("cfl".to_string(), 0.5);
    params.insert("init_type".to_string(), 1.0); // Taylor-Green
    params.insert("u0".to_string(), 1.0);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::NavierStokes3D,
        Some([0.0, 0.2]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate 3D Taylor-Green: {:?}", result);

    if let Ok(ToolResponse::Simulate(output)) = result {
        let moments = output.moments.expect("Should have moments");
        assert!(moments.contains_key("kinetic_energy"));
        assert!(moments.contains_key("enstrophy"));
        assert!(moments.contains_key("reynolds_number"));
    }
}

#[test]
fn test_ns_3d_uniform_flow() {
    let dispatcher = create_default_dispatcher();

    let mut params = HashMap::new();
    params.insert("nx".to_string(), 8.0);
    params.insert("ny".to_string(), 8.0);
    params.insert("nz".to_string(), 8.0);
    params.insert("viscosity".to_string(), 0.1);
    params.insert("cfl".to_string(), 0.5);
    params.insert("init_type".to_string(), 0.0); // Uniform
    params.insert("u0".to_string(), 1.0);
    params.insert("v0".to_string(), 0.0);
    params.insert("w0".to_string(), 0.0);

    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::NavierStokes3D,
        Some([0.0, 0.1]),
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_ok(), "Should simulate uniform flow: {:?}", result);
}

// ============================================================================
// NAVIER-STOKES SOLVER (DIRECT) TESTS
// ============================================================================

#[test]
fn test_navier_stokes_solver_basic() {
    use serde_json::json;

    let config = json!({
        "domain": {
            "width": 1.0,
            "height": 1.0,
            "grid_points_x": 11,
            "grid_points_y": 11
        },
        "fluid_properties": {
            "viscosity": 0.01,
            "density": 1.0
        },
        "boundary_conditions": {
            "top": "no_slip",
            "bottom": "no_slip",
            "left": "no_slip",
            "right": "no_slip"
        },
        "initial_conditions": {
            "velocity_x": 0.0,
            "velocity_y": 0.0,
            "pressure": 0.0
        },
        "simulation_parameters": {
            "time_step": 0.001,
            "total_time": 0.01
        }
    });

    let solver = NavierStokesSolver::new(config);
    let result = solver.solve();
    assert!(result.is_ok(), "Should solve Navier-Stokes: {:?}", result);

    let output = result.unwrap();
    assert!(output["solution"]["velocity_field"]["u"].is_array());
    assert!(output["solution"]["velocity_field"]["v"].is_array());
    assert!(output["solution"]["pressure_field"].is_array());
    assert!(output["solution"]["derived_quantities"]["vorticity"].is_array());
    assert!(output["solution"]["derived_quantities"]["divergence"].is_array());
}

#[test]
fn test_navier_stokes_solver_channel_flow() {
    use serde_json::json;

    let config = json!({
        "domain": {
            "width": 2.0,
            "height": 1.0,
            "grid_points_x": 21,
            "grid_points_y": 11
        },
        "fluid_properties": {
            "viscosity": 0.1,
            "density": 1.0
        },
        "boundary_conditions": {
            "top": "no_slip",
            "bottom": "no_slip",
            "left": "velocity_inlet",
            "right": "pressure_outlet"
        },
        "initial_conditions": {
            "velocity_x": 0.0,
            "velocity_y": 0.0,
            "pressure": 0.0
        },
        "simulation_parameters": {
            "time_step": 0.001,
            "total_time": 0.02
        }
    });

    let solver = NavierStokesSolver::new(config);
    let result = solver.solve();
    assert!(result.is_ok(), "Should solve channel flow: {:?}", result);

    let output = result.unwrap();
    let sim_params = &output["simulation_parameters"];
    assert_eq!(sim_params["grid"]["nx"], 21);
    assert_eq!(sim_params["grid"]["ny"], 11);
}

#[test]
fn test_navier_stokes_solver_free_slip() {
    use serde_json::json;

    let config = json!({
        "domain": {
            "width": 1.0,
            "height": 1.0,
            "grid_points_x": 11,
            "grid_points_y": 11
        },
        "fluid_properties": {
            "viscosity": 0.01,
            "density": 1.0
        },
        "boundary_conditions": {
            "top": "free_slip",
            "bottom": "free_slip",
            "left": "no_slip",
            "right": "no_slip"
        },
        "initial_conditions": {
            "velocity_x": 0.1,
            "velocity_y": 0.0,
            "pressure": 0.0
        },
        "simulation_parameters": {
            "time_step": 0.001,
            "total_time": 0.01
        }
    });

    let solver = NavierStokesSolver::new(config);
    let result = solver.solve();
    assert!(result.is_ok(), "Should handle free slip: {:?}", result);
}

// ============================================================================
// GRID TESTS
// ============================================================================

#[test]
fn test_grid_creation() {
    let grid = Grid2D::uniform(11, 11, 1.0, 1.0);
    assert_eq!(grid.nx, 11);
    assert_eq!(grid.ny, 11);
    assert!((grid.dx - 0.1).abs() < 1e-10);
    assert!((grid.dy - 0.1).abs() < 1e-10);
}

#[test]
fn test_grid_coordinates() {
    let grid = Grid2D::uniform(11, 11, 1.0, 1.0);
    assert!((grid.x_coord(0) - 0.0).abs() < 1e-10);
    assert!((grid.x_coord(10) - 1.0).abs() < 1e-10);
    assert!((grid.y_coord(5) - 0.5).abs() < 1e-10);
}

#[test]
fn test_grid_boundary_detection() {
    let grid = Grid2D::uniform(5, 5, 1.0, 1.0);
    assert!(grid.is_boundary(0, 2));
    assert!(grid.is_boundary(4, 2));
    assert!(grid.is_boundary(2, 0));
    assert!(grid.is_boundary(2, 4));
    assert!(grid.is_interior(2, 2));
    assert!(!grid.is_interior(0, 0));
}

#[test]
fn test_grid_non_square() {
    let grid = Grid2D::uniform(21, 11, 2.0, 1.0);
    assert_eq!(grid.nx, 21);
    assert_eq!(grid.ny, 11);
    assert!((grid.dx - 0.1).abs() < 1e-10);
    assert!((grid.dy - 0.1).abs() < 1e-10);
    assert!((grid.x_coord(20) - 2.0).abs() < 1e-10);
    assert!((grid.y_coord(10) - 1.0).abs() < 1e-10);
}

// ============================================================================
// BOUNDARY CONDITIONS TESTS
// ============================================================================

#[test]
fn test_lid_driven_cavity_bc() {
    use ndarray::Array2;

    let bc = BoundaryConditions::lid_driven_cavity(1.0);
    let mut u = Array2::<f64>::zeros((5, 5));
    let mut v = Array2::<f64>::zeros((5, 5));

    bc.apply_velocity_boundary(&mut u, &mut v);

    // Check top boundary has lid velocity
    for i in 0..5 {
        assert_eq!(u[[i, 4]], 1.0);
        assert_eq!(v[[i, 4]], 0.0);
    }

    // Check other boundaries are no-slip
    for i in 0..5 {
        assert_eq!(u[[i, 0]], 0.0);
        assert_eq!(v[[i, 0]], 0.0);
    }
}

#[test]
fn test_channel_flow_bc() {
    use ndarray::Array2;

    let bc = BoundaryConditions::channel_flow(2.0);
    let mut u = Array2::<f64>::zeros((5, 5));
    let mut v = Array2::<f64>::zeros((5, 5));

    bc.apply_velocity_boundary(&mut u, &mut v);

    // Interior inlet points should have inlet velocity
    for j in 1..4 {
        assert_eq!(u[[0, j]], 2.0);
        assert_eq!(v[[0, j]], 0.0);
    }
}

#[test]
fn test_no_slip_bc() {
    use ndarray::Array2;

    let bc = BoundaryConditions::all_no_slip();
    let mut u = Array2::<f64>::from_elem((5, 5), 1.0);
    let mut v = Array2::<f64>::from_elem((5, 5), 1.0);

    bc.apply_velocity_boundary(&mut u, &mut v);

    // All boundaries should be zero
    for i in 0..5 {
        assert_eq!(u[[i, 0]], 0.0); // bottom
        assert_eq!(v[[i, 0]], 0.0);
        assert_eq!(u[[i, 4]], 0.0); // top
        assert_eq!(v[[i, 4]], 0.0);
    }
    for j in 0..5 {
        assert_eq!(u[[0, j]], 0.0); // left
        assert_eq!(v[[0, j]], 0.0);
        assert_eq!(u[[4, j]], 0.0); // right
        assert_eq!(v[[4, j]], 0.0);
    }
}

// ============================================================================
// ANALYTICAL SOLUTIONS TESTS
// ============================================================================

#[test]
fn test_analytical_couette_flow() {
    use serde_json::json;

    let solver = AnalyticalSolutions::new();
    let params = json!({
        "plate_velocity": 1.0,
        "plate_separation": 1.0,
        "num_points": 11,
        "viscosity": 0.001
    });

    let result = solver.solve("couette_flow", params);
    assert!(result.is_ok(), "Should solve Couette flow: {:?}", result);

    let output = result.unwrap();
    assert_eq!(output["solution_type"], "couette_flow");
    assert!(output["numerical_data"]["y_coordinates"].is_array());
    assert!(output["numerical_data"]["u_velocity"].is_array());

    // Check that velocity at bottom wall is 0 and at top wall is plate_velocity
    let velocities = output["numerical_data"]["u_velocity"].as_array().unwrap();
    assert!((velocities[0].as_f64().unwrap() - 0.0).abs() < 1e-10);
    assert!((velocities[10].as_f64().unwrap() - 1.0).abs() < 1e-10);
}

#[test]
fn test_analytical_poiseuille_flow() {
    use serde_json::json;

    let solver = AnalyticalSolutions::new();
    let params = json!({
        "pressure_gradient": -1.0,
        "channel_height": 1.0,
        "viscosity": 0.1,
        "num_points": 11
    });

    let result = solver.solve("poiseuille_flow", params);
    assert!(result.is_ok(), "Should solve Poiseuille flow: {:?}", result);

    let output = result.unwrap();
    assert_eq!(output["solution_type"], "poiseuille_flow");

    // Check that velocity is maximum at center and zero at walls
    let velocities = output["numerical_data"]["u_velocity"].as_array().unwrap();
    assert!((velocities[0].as_f64().unwrap() - 0.0).abs() < 1e-10); // bottom wall
    assert!((velocities[10].as_f64().unwrap() - 0.0).abs() < 1e-10); // top wall

    // Maximum should be at center (index 5)
    let max_vel = velocities[5].as_f64().unwrap();
    assert!(max_vel > 0.0, "Center velocity should be positive");
}

#[test]
fn test_analytical_stagnation_point_flow() {
    use serde_json::json;

    let solver = AnalyticalSolutions::new();
    let params = json!({
        "strain_rate": 1.0,
        "domain_size": 2.0,
        "num_points": 5,
        "density": 1.0
    });

    let result = solver.solve("stagnation_point_flow", params);
    assert!(result.is_ok(), "Should solve stagnation point flow: {:?}", result);

    let output = result.unwrap();
    assert_eq!(output["solution_type"], "stagnation_point_flow");
    assert_eq!(output["flow_characteristics"]["flow_type"], "irrotational");
    assert_eq!(output["flow_characteristics"]["vorticity"], 0.0);
}

#[test]
fn test_analytical_potential_vortex() {
    use serde_json::json;

    let solver = AnalyticalSolutions::new();
    let params = json!({
        "circulation": 1.0,
        "core_radius": 0.1,
        "radius_max": 2.0,
        "num_points": 50
    });

    let result = solver.solve("potential_vortex", params);
    assert!(result.is_ok(), "Should solve potential vortex: {:?}", result);

    let output = result.unwrap();
    assert_eq!(output["solution_type"], "potential_vortex");
    assert_eq!(output["flow_characteristics"]["flow_type"], "irrotational");
}

#[test]
fn test_analytical_rankine_vortex() {
    use serde_json::json;

    let solver = AnalyticalSolutions::new();
    let params = json!({
        "circulation": 1.0,
        "core_radius": 0.5,
        "radius_max": 3.0,
        "num_points": 100
    });

    let result = solver.solve("rankine_vortex", params);
    assert!(result.is_ok(), "Should solve Rankine vortex: {:?}", result);

    let output = result.unwrap();
    assert_eq!(output["solution_type"], "rankine_vortex");
    assert!(output["analytical_solution"]["core_vorticity"].as_f64().is_some());
}

// ============================================================================
// FLOW ANALYZER TESTS
// ============================================================================

#[test]
fn test_flow_analyzer_velocity_magnitude() {
    use serde_json::json;

    let analyzer = FlowAnalyzer::new();

    // Create a simple 3x3 velocity field
    let velocity_field = json!({
        "u_velocity": [[1.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 1.0]],
        "v_velocity": [[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
    });

    let grid_spacing = json!({
        "dx": 0.1,
        "dy": 0.1
    });

    let result = analyzer.analyze(velocity_field, grid_spacing, vec!["velocity_magnitude".to_string()]);
    assert!(result.is_ok(), "Should calculate velocity magnitude: {:?}", result);

    let output = result.unwrap();
    assert!(output["analysis_results"]["velocity_magnitude"].is_array());
}

#[test]
fn test_flow_analyzer_vorticity() {
    use serde_json::json;

    let analyzer = FlowAnalyzer::new();

    // Create a 5x5 velocity field (minimum needed for vorticity)
    let velocity_field = json!({
        "u_velocity": [
            [0.0, 0.1, 0.2, 0.1, 0.0],
            [0.1, 0.2, 0.3, 0.2, 0.1],
            [0.2, 0.3, 0.4, 0.3, 0.2],
            [0.1, 0.2, 0.3, 0.2, 0.1],
            [0.0, 0.1, 0.2, 0.1, 0.0]
        ],
        "v_velocity": [
            [0.0, -0.1, -0.2, -0.1, 0.0],
            [-0.1, -0.2, -0.3, -0.2, -0.1],
            [-0.2, -0.3, -0.4, -0.3, -0.2],
            [-0.1, -0.2, -0.3, -0.2, -0.1],
            [0.0, -0.1, -0.2, -0.1, 0.0]
        ]
    });

    let grid_spacing = json!({
        "dx": 0.1,
        "dy": 0.1
    });

    let result = analyzer.analyze(velocity_field, grid_spacing, vec!["vorticity".to_string()]);
    assert!(result.is_ok(), "Should calculate vorticity: {:?}", result);
}

#[test]
fn test_flow_analyzer_divergence() {
    use serde_json::json;

    let analyzer = FlowAnalyzer::new();

    // Create incompressible flow field (divergence should be near zero)
    let velocity_field = json!({
        "u_velocity": [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0, 2.0, 2.0],
            [1.0, 1.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 0.0]
        ],
        "v_velocity": [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0]
        ]
    });

    let grid_spacing = json!({
        "dx": 0.1,
        "dy": 0.1
    });

    let result = analyzer.analyze(velocity_field, grid_spacing, vec!["divergence".to_string()]);
    assert!(result.is_ok(), "Should calculate divergence: {:?}", result);
}

#[test]
fn test_flow_analyzer_kinetic_energy() {
    use serde_json::json;

    let analyzer = FlowAnalyzer::new();

    // Create a uniform flow field
    let velocity_field = json!({
        "u_velocity": [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],
        "v_velocity": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    });

    let grid_spacing = json!({
        "dx": 0.1,
        "dy": 0.1
    });

    let result = analyzer.analyze(velocity_field, grid_spacing, vec!["kinetic_energy".to_string()]);
    assert!(result.is_ok(), "Should calculate kinetic energy: {:?}", result);

    let output = result.unwrap();
    let ke = &output["analysis_results"]["kinetic_energy"];
    assert!(ke["total"].as_f64().is_some());
    // For uniform flow with u=1, v=0: KE = 0.5 * (1^2 + 0^2) = 0.5 per point
    let total = ke["total"].as_f64().unwrap();
    assert!((total - 0.5).abs() < 1e-10);
}

// ============================================================================
// ERROR HANDLING TESTS
// ============================================================================

#[test]
fn test_qns_requires_range() {
    let dispatcher = create_default_dispatcher();

    let params = HashMap::new();
    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::QuantumNavierStokes1D,
        None, // No range
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "QNS 1D should require range");
}

#[test]
fn test_ns_3d_requires_range() {
    let dispatcher = create_default_dispatcher();

    let params = HashMap::new();
    let request = ToolRequest::Simulate(create_fluid_input(
        FluidSim::NavierStokes3D,
        None, // No range
        None,
        params,
    ));

    let result = dispatcher.dispatch(request);
    assert!(result.is_err(), "NS 3D should require range");
}
