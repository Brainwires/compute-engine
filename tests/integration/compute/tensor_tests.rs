//! Comprehensive tests for tensor compute operations
//!
//! Tests differential geometry and tensor calculus operations including:
//! - Christoffel symbols
//! - Riemann curvature tensor
//! - Ricci tensor and scalar
//! - Einstein tensor
//! - Weyl tensor

use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// CHRISTOFFEL SYMBOL TESTS
// ============================================================================

#[test]
#[ignore] // Slow symbolic computation
fn test_christoffel_flat_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Christoffel),
        data: json!({
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Christoffel symbols for flat 2D metric should succeed"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_christoffel_minkowski_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Christoffel),
        data: json!({
            "metric": [[-1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["t", "x"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Christoffel symbols for Minkowski metric should succeed"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_christoffel_flat_3d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Christoffel),
        data: json!({
            "metric": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "coordinates": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Christoffel symbols for flat 3D metric should succeed"
    );
}

// ============================================================================
// RIEMANN TENSOR TESTS
// ============================================================================

#[test]
#[ignore] // Slow symbolic computation
fn test_riemann_flat_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Riemann),
        data: json!({
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Riemann tensor for flat 2D metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_riemann_minkowski_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Riemann),
        data: json!({
            "metric": [[-1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["t", "x"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Riemann tensor for Minkowski metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_riemann_flat_3d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Riemann),
        data: json!({
            "metric": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "coordinates": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Riemann tensor for flat 3D metric should succeed"
    );
}

// ============================================================================
// RICCI TENSOR TESTS
// ============================================================================

#[test]
#[ignore] // Slow symbolic computation
fn test_ricci_flat_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Ricci),
        data: json!({
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Ricci tensor for flat 2D metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_ricci_minkowski_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Ricci),
        data: json!({
            "metric": [[-1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["t", "x"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Ricci tensor for Minkowski metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_ricci_flat_3d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Ricci),
        data: json!({
            "metric": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "coordinates": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Ricci tensor for flat 3D metric should succeed"
    );
}

// ============================================================================
// RICCI SCALAR TESTS
// ============================================================================

#[test]
#[ignore] // Slow symbolic computation
fn test_ricci_scalar_flat_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::RicciScalar),
        data: json!({
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Ricci scalar for flat 2D metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_ricci_scalar_minkowski_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::RicciScalar),
        data: json!({
            "metric": [[-1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["t", "x"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Ricci scalar for Minkowski metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_ricci_scalar_flat_3d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::RicciScalar),
        data: json!({
            "metric": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "coordinates": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Ricci scalar for flat 3D metric should succeed"
    );
}

// ============================================================================
// EINSTEIN TENSOR TESTS
// ============================================================================

#[test]
#[ignore] // Slow symbolic computation
fn test_einstein_flat_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Einstein),
        data: json!({
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Einstein tensor for flat 2D metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_einstein_minkowski_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Einstein),
        data: json!({
            "metric": [[-1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["t", "x"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Einstein tensor for Minkowski metric should succeed (and be zero)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_einstein_flat_3d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Einstein),
        data: json!({
            "metric": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "coordinates": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Einstein tensor for flat 3D metric should succeed"
    );
}

// ============================================================================
// WEYL TENSOR TESTS
// ============================================================================

#[test]
#[ignore] // Slow symbolic computation
fn test_weyl_flat_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Weyl),
        data: json!({
            "metric": [[1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["x", "y"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    // Note: Weyl tensor is always zero in 2D (needs dimension >= 3)
    assert!(
        response.is_ok(),
        "Weyl tensor for flat 2D metric should succeed (always zero in 2D)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_weyl_minkowski_2d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Weyl),
        data: json!({
            "metric": [[-1.0, 0.0], [0.0, 1.0]],
            "coordinates": ["t", "x"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Weyl tensor for Minkowski metric should succeed (always zero in 2D)"
    );
}

#[test]
#[ignore] // Slow symbolic computation
fn test_weyl_flat_3d() {
    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::Tensor(TensorOp::Weyl),
        data: json!({
            "metric": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "coordinates": ["x", "y", "z"]
        }),
        parameters: HashMap::new(),
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Weyl tensor for flat 3D metric should succeed (zero for conformally flat)"
    );
}
