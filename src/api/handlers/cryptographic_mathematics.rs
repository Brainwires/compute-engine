//! Cryptographic mathematics operations handler

use crate::api::types::{ComputationRequest, ComputationResponse};
use crate::specialized::cryptographic_mathematics::*;
use num_bigint::BigInt;
use serde_json::json;
use std::str::FromStr;

pub fn handle(request: &ComputationRequest) -> ComputationResponse {
    let result = match request.operation.as_str() {
        "rsa_keygen" => {
            let bits: u32 = match request.parameters.get("bit_length") {
                Some(v) => match v.as_u64() {
                    Some(b) => b as u32,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid bit_length parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing bit_length parameter".to_string(),
                    );
                }
            };
            let (n, e, d) = generate_rsa_keypair(bits);
            Ok(json!({"n": n.to_string(), "e": e.to_string(), "d": d.to_string()}))
        }
        "generate_prime" => {
            let bits: u32 = match request.parameters.get("bit_length") {
                Some(v) => match v.as_u64() {
                    Some(b) => b as u32,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid bit_length parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing bit_length parameter".to_string(),
                    );
                }
            };
            let prime = generate_prime(bits);
            Ok(json!({"prime": prime.to_string()}))
        }
        "modular_exponentiation" => {
            let base_str = match request.parameters.get("base") {
                Some(v) => match v.as_str() {
                    Some(s) => s,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid base parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing base parameter".to_string(),
                    );
                }
            };
            let exp_str = match request.parameters.get("exponent") {
                Some(v) => match v.as_str() {
                    Some(s) => s,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid exponent parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing exponent parameter".to_string(),
                    );
                }
            };
            let mod_str = match request.parameters.get("modulus") {
                Some(v) => match v.as_str() {
                    Some(s) => s,
                    None => {
                        return ComputationResponse::error(
                            request.module.clone(),
                            request.operation.clone(),
                            "Invalid modulus parameter".to_string(),
                        );
                    }
                },
                None => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        "Missing modulus parameter".to_string(),
                    );
                }
            };

            let base = match BigInt::from_str(base_str) {
                Ok(b) => b,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid base: {}", e),
                    );
                }
            };
            let exp = match BigInt::from_str(exp_str) {
                Ok(e) => e,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid exponent: {}", e),
                    );
                }
            };
            let modulus = match BigInt::from_str(mod_str) {
                Ok(m) => m,
                Err(e) => {
                    return ComputationResponse::error(
                        request.module.clone(),
                        request.operation.clone(),
                        format!("Invalid modulus: {}", e),
                    );
                }
            };

            let result_val = mod_exp(&base, &exp, &modulus);
            Ok(json!({"result": result_val.to_string()}))
        }
        _ => {
            return ComputationResponse::error(
                request.module.clone(),
                request.operation.clone(),
                format!("Unknown operation: {}", request.operation),
            );
        }
    };

    match result {
        Ok(result_value) => ComputationResponse::success(
            request.module.clone(),
            request.operation.clone(),
            result_value,
        ),
        Err(error_msg) => {
            ComputationResponse::error(request.module.clone(), request.operation.clone(), error_msg)
        }
    }
}

// TODO: Fix test compilation error
// #[cfg(test)]
// #[path = "../../../tests/unit/api/handlers/cryptographic_mathematics_handler_tests.rs"]
// mod tests;
