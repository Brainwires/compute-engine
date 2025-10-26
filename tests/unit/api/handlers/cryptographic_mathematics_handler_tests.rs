//! Unit tests for cryptographic_mathematics API handler

use super::handle; // handle is in the parent module (cryptographic_mathematics)
use crate::api::types::ComputationRequest;
use num_bigint::BigInt;
use serde_json::json;
use std::collections::HashMap;
use std::str::FromStr;

// Helper function to create a request
fn create_request(operation: &str, parameters: serde_json::Value) -> ComputationRequest {
    let params_map = if let serde_json::Value::Object(map) = parameters {
        map.into_iter().collect()
    } else {
        HashMap::new()
    };

    ComputationRequest {
        module: "cryptographic_mathematics".to_string(),
        operation: operation.to_string(),
        parameters: params_map,
    }
}

// ============================================================================
// RSA Keygen Tests
// ============================================================================

#[test]
fn test_rsa_keygen_success_64_bits() {
    let request = create_request("rsa_keygen", json!({"bit_length": 64}));
    let response = handle(&request);

    assert!(response.success, "RSA keygen should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    assert!(result.get("n").is_some(), "Should have modulus n");
    assert!(result.get("e").is_some(), "Should have public exponent e");
    assert!(result.get("d").is_some(), "Should have private exponent d");

    // Verify n, e, d are valid BigInt strings
    let n_str = result["n"].as_str().unwrap();
    let e_str = result["e"].as_str().unwrap();
    let d_str = result["d"].as_str().unwrap();

    assert!(BigInt::from_str(n_str).is_ok(), "n should be valid BigInt");
    assert!(BigInt::from_str(e_str).is_ok(), "e should be valid BigInt");
    assert!(BigInt::from_str(d_str).is_ok(), "d should be valid BigInt");
}

#[test]
fn test_rsa_keygen_success_128_bits() {
    let request = create_request("rsa_keygen", json!({"bit_length": 128}));
    let response = handle(&request);

    assert!(response.success, "RSA keygen should succeed");
    assert!(response.error.is_none(), "Should have no error");
}

#[test]
fn test_rsa_keygen_missing_bit_length() {
    let request = create_request("rsa_keygen", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail without bit_length");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Missing bit_length parameter",
        "Should report missing parameter"
    );
}

#[test]
fn test_rsa_keygen_invalid_bit_length_type() {
    let request = create_request("rsa_keygen", json!({"bit_length": "invalid"}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid type");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Invalid bit_length parameter",
        "Should report invalid parameter"
    );
}

#[test]
fn test_rsa_keygen_invalid_bit_length_float() {
    let request = create_request("rsa_keygen", json!({"bit_length": 64.5}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with float bit_length");
    assert!(response.error.is_some(), "Should have error message");
}

// ============================================================================
// Generate Prime Tests
// ============================================================================

#[test]
fn test_generate_prime_success_32_bits() {
    let request = create_request("generate_prime", json!({"bit_length": 32}));
    let response = handle(&request);

    assert!(response.success, "Generate prime should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    assert!(result.get("prime").is_some(), "Should have prime field");

    let prime_str = result["prime"].as_str().unwrap();
    let prime = BigInt::from_str(prime_str).unwrap();

    // Basic sanity check: prime should be positive
    assert!(prime > BigInt::from(0), "Prime should be positive");
}

#[test]
fn test_generate_prime_success_64_bits() {
    let request = create_request("generate_prime", json!({"bit_length": 64}));
    let response = handle(&request);

    assert!(response.success, "Generate prime should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    let prime_str = result["prime"].as_str().unwrap();
    assert!(BigInt::from_str(prime_str).is_ok(), "Prime should be valid BigInt");
}

#[test]
fn test_generate_prime_missing_bit_length() {
    let request = create_request("generate_prime", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail without bit_length");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Missing bit_length parameter",
        "Should report missing parameter"
    );
}

#[test]
fn test_generate_prime_invalid_bit_length_type() {
    let request = create_request("generate_prime", json!({"bit_length": null}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with null bit_length");
    assert!(response.error.is_some(), "Should have error message");
    assert_eq!(
        response.error.unwrap(),
        "Invalid bit_length parameter",
        "Should report invalid parameter"
    );
}

// ============================================================================
// Modular Exponentiation Tests
// ============================================================================

#[test]
fn test_modular_exponentiation_success_small_numbers() {
    // 2^3 mod 5 = 8 mod 5 = 3
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "exponent": "3",
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Modular exponentiation should succeed");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    assert_eq!(result["result"].as_str().unwrap(), "3", "2^3 mod 5 = 3");
}

#[test]
fn test_modular_exponentiation_success_large_numbers() {
    // 123^456 mod 789
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "123",
            "exponent": "456",
            "modulus": "789"
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Should handle large numbers");
    assert!(response.error.is_none(), "Should have no error");

    let result = response.result.unwrap();
    let result_str = result["result"].as_str().unwrap();
    assert!(BigInt::from_str(result_str).is_ok(), "Result should be valid BigInt");
}

#[test]
fn test_modular_exponentiation_zero_exponent() {
    // Any number to the power of 0 is 1
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "999",
            "exponent": "0",
            "modulus": "100"
        }),
    );
    let response = handle(&request);

    assert!(response.success, "Should succeed with zero exponent");
    assert_eq!(response.result.unwrap()["result"].as_str().unwrap(), "1");
}

#[test]
fn test_modular_exponentiation_missing_base() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "exponent": "3",
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail without base");
    assert_eq!(
        response.error.unwrap(),
        "Missing base parameter",
        "Should report missing base"
    );
}

#[test]
fn test_modular_exponentiation_missing_exponent() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail without exponent");
    assert_eq!(
        response.error.unwrap(),
        "Missing exponent parameter",
        "Should report missing exponent"
    );
}

#[test]
fn test_modular_exponentiation_missing_modulus() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "exponent": "3"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail without modulus");
    assert_eq!(
        response.error.unwrap(),
        "Missing modulus parameter",
        "Should report missing modulus"
    );
}

#[test]
fn test_modular_exponentiation_invalid_base_type() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": 123,  // Number instead of string
            "exponent": "3",
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with non-string base");
    assert_eq!(
        response.error.unwrap(),
        "Invalid base parameter",
        "Should report invalid base"
    );
}

#[test]
fn test_modular_exponentiation_invalid_exponent_type() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "exponent": true,  // Boolean instead of string
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with non-string exponent");
    assert_eq!(
        response.error.unwrap(),
        "Invalid exponent parameter",
        "Should report invalid exponent"
    );
}

#[test]
fn test_modular_exponentiation_invalid_modulus_type() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "exponent": "3",
            "modulus": []  // Array instead of string
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with non-string modulus");
    assert_eq!(
        response.error.unwrap(),
        "Invalid modulus parameter",
        "Should report invalid modulus"
    );
}

#[test]
fn test_modular_exponentiation_invalid_base_format() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "not_a_number",
            "exponent": "3",
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid base format");
    assert!(
        response.error.unwrap().starts_with("Invalid base:"),
        "Should report invalid base format"
    );
}

#[test]
fn test_modular_exponentiation_invalid_exponent_format() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "exponent": "3.14",
            "modulus": "5"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid exponent format");
    assert!(
        response.error.unwrap().starts_with("Invalid exponent:"),
        "Should report invalid exponent format"
    );
}

#[test]
fn test_modular_exponentiation_invalid_modulus_format() {
    let request = create_request(
        "modular_exponentiation",
        json!({
            "base": "2",
            "exponent": "3",
            "modulus": "abc"
        }),
    );
    let response = handle(&request);

    assert!(!response.success, "Should fail with invalid modulus format");
    assert!(
        response.error.unwrap().starts_with("Invalid modulus:"),
        "Should report invalid modulus format"
    );
}

// ============================================================================
// Unknown Operation Tests
// ============================================================================

#[test]
fn test_unknown_operation() {
    let request = create_request("unknown_operation", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with unknown operation");
    assert!(response.error.is_some(), "Should have error message");
    assert!(
        response
            .error
            .unwrap()
            .contains("Unknown operation: unknown_operation"),
        "Should report unknown operation"
    );
}

#[test]
fn test_empty_operation() {
    let request = create_request("", json!({}));
    let response = handle(&request);

    assert!(!response.success, "Should fail with empty operation");
    assert!(response.error.is_some(), "Should have error message");
}
