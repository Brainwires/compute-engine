use super::super::{ComputationRequest, ComputationResponse};

pub fn process_calculus_request(request: &ComputationRequest) -> ComputationResponse {
    use crate::advanced_calculus::*;

    // Helper macro to convert result types to ComputationResponse
    macro_rules! handle_op {
        ($handler:expr) => {{
            let result = $handler;
            ComputationResponse {
                success: result.success,
                module: request.module.clone(),
                operation: request.operation.clone(),
                result: result.result,
                error: result.error,
            }
        }};
    }

    match request.operation.as_str() {
        "fractional_derivative" => handle_op!(handle_fractional_derivative(&request.parameters)),
        "fractional_integral" => handle_op!(handle_fractional_integral(&request.parameters)),
        "fractional_calculus" => handle_op!(handle_fractional_calculus(&request.parameters)),
        "riemann_zeta" => handle_op!(handle_riemann_zeta(&request.parameters)),
        "elliptic_integral" => handle_op!(handle_elliptic_integral(&request.parameters)),
        "hypergeometric" => handle_op!(handle_hypergeometric(&request.parameters)),
        "jacobi_theta" => handle_op!(handle_jacobi_theta(&request.parameters)),
        "bessel_function" => handle_op!(handle_bessel_function(&request.parameters)),
        "legendre_polynomial" => handle_op!(handle_legendre_polynomial(&request.parameters)),
        "special_functions" => handle_op!(handle_special_functions(&request.parameters)),
        "euler_lagrange" => handle_op!(handle_euler_lagrange(&request.parameters)),
        "variational_calculus" => handle_op!(handle_variational_calculus(&request.parameters)),
        "ito_integral" => handle_op!(handle_ito_integral(&request.parameters)),
        "stratonovich_integral" => handle_op!(handle_stratonovich_integral(&request.parameters)),
        "sde_solution" => handle_op!(handle_sde_solution(&request.parameters)),
        "stochastic_calculus" => handle_op!(handle_stochastic_calculus_ops(&request.parameters)),
        "symbolic_integral" => handle_op!(handle_symbolic_integral(&request.parameters)),
        "definite_integral" => handle_op!(handle_definite_integral(&request.parameters)),
        "improper_integral" => handle_op!(handle_improper_integral(&request.parameters)),
        _ => ComputationResponse {
            success: false,
            module: request.module.clone(),
            operation: request.operation.clone(),
            result: None,
            error: Some(format!("Unknown calculus operation: {}", request.operation)),
        },
    }
}
