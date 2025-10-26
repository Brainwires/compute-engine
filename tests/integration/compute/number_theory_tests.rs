//! Comprehensive tests for number theory compute operations
//!
//! Tests number theory and cryptographic operations including:
//! - Prime generation and testing
//! - Modular arithmetic (exponentiation, inverse)
//! - GCD and LCM
//! - RSA cryptography
//! - Hashing functions

use computational_engine::engine::*;
use serde_json::json;
use std::collections::HashMap;

// ============================================================================
// PRIME GENERATION TESTS
// ============================================================================

#[test]
fn test_generate_prime_small() {
    let mut parameters = HashMap::new();
    parameters.insert("bit_length".to_string(), json!(8));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::GeneratePrime),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "8-bit prime generation should succeed");
}

#[test]
fn test_generate_prime_medium() {
    let mut parameters = HashMap::new();
    parameters.insert("bit_length".to_string(), json!(64));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::GeneratePrime),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "64-bit prime generation should succeed");
}

#[test]
fn test_generate_prime_large() {
    let mut parameters = HashMap::new();
    parameters.insert("bit_length".to_string(), json!(128));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::GeneratePrime),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "128-bit prime generation should succeed");
}

// ============================================================================
// MODULAR EXPONENTIATION TESTS
// ============================================================================

#[test]
fn test_modular_exponentiation_small() {
    let mut parameters = HashMap::new();
    parameters.insert("base".to_string(), json!("5"));
    parameters.insert("exponent".to_string(), json!("3"));
    parameters.insert("modulus".to_string(), json!("13"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ModExp),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Small modular exponentiation should succeed");
}

#[test]
fn test_modular_exponentiation_large() {
    let mut parameters = HashMap::new();
    parameters.insert("base".to_string(), json!("123456789"));
    parameters.insert("exponent".to_string(), json!("987654321"));
    parameters.insert("modulus".to_string(), json!("1000000007"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ModExp),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Large modular exponentiation should succeed");
}

// ============================================================================
// MODULAR INVERSE TESTS
// ============================================================================

#[test]
fn test_modular_inverse_exists() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("3"));
    parameters.insert("modulus".to_string(), json!("11"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ModInv),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Modular inverse should succeed when inverse exists"
    );
}

#[test]
fn test_modular_inverse_not_exists() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("6"));
    parameters.insert("modulus".to_string(), json!("9"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ModInv),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_err(),
        "Modular inverse should fail when inverse doesn't exist"
    );
}

// ============================================================================
// GCD TESTS
// ============================================================================

#[test]
fn test_gcd_coprime() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("15"));
    parameters.insert("b".to_string(), json!("28"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::GCD),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "GCD of coprime numbers should succeed");
}

#[test]
fn test_gcd_common_divisor() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("48"));
    parameters.insert("b".to_string(), json!("18"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::GCD),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "GCD with common divisor should succeed");
}

#[test]
fn test_gcd_large_numbers() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("123456789012345"));
    parameters.insert("b".to_string(), json!("987654321098765"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::GCD),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "GCD of large numbers should succeed");
}

// ============================================================================
// LCM TESTS
// ============================================================================

#[test]
fn test_lcm_small() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("4"));
    parameters.insert("b".to_string(), json!("6"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::LCM),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "LCM of small numbers should succeed");
}

#[test]
fn test_lcm_coprime() {
    let mut parameters = HashMap::new();
    parameters.insert("a".to_string(), json!("7"));
    parameters.insert("b".to_string(), json!("13"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::LCM),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "LCM of coprime numbers should succeed");
}

// ============================================================================
// EULER TOTIENT TESTS
// ============================================================================

#[test]
fn test_euler_totient_prime() {
    let mut parameters = HashMap::new();
    parameters.insert("n".to_string(), json!("17"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::EulerTotient),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Euler totient of prime should succeed");
}

#[test]
fn test_euler_totient_composite() {
    let mut parameters = HashMap::new();
    parameters.insert("n".to_string(), json!("36"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::EulerTotient),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Euler totient of composite should succeed");
}

// ============================================================================
// CARMICHAEL LAMBDA TESTS
// ============================================================================

#[test]
fn test_carmichael_lambda_prime() {
    let mut parameters = HashMap::new();
    parameters.insert("n".to_string(), json!("11"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::CarmichaelLambda),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "Carmichael lambda of prime should succeed");
}

#[test]
fn test_carmichael_lambda_composite() {
    let mut parameters = HashMap::new();
    parameters.insert("n".to_string(), json!("15"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::CarmichaelLambda),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Carmichael lambda of composite should succeed"
    );
}

// ============================================================================
// RSA KEYPAIR TESTS
// ============================================================================

#[test]
fn test_rsa_keypair_small() {
    let mut parameters = HashMap::new();
    parameters.insert("bit_length".to_string(), json!(512));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::RSAKeypair),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "512-bit RSA keypair generation should succeed");
}

#[test]
fn test_rsa_keypair_medium() {
    let mut parameters = HashMap::new();
    parameters.insert("bit_length".to_string(), json!(1024));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::RSAKeypair),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "1024-bit RSA keypair generation should succeed"
    );
}

// ============================================================================
// RSA ENCRYPT/DECRYPT TESTS
// ============================================================================

#[test]
fn test_rsa_encrypt() {
    let mut parameters = HashMap::new();
    parameters.insert("message".to_string(), json!("42"));
    parameters.insert("e".to_string(), json!("65537"));
    parameters.insert("n".to_string(), json!("3233")); // Example modulus

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::RSAEncrypt),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "RSA encryption should succeed");
}

#[test]
fn test_rsa_decrypt() {
    let mut parameters = HashMap::new();
    parameters.insert("ciphertext".to_string(), json!("2557")); // Example ciphertext
    parameters.insert("d".to_string(), json!("2753")); // Example private exponent
    parameters.insert("n".to_string(), json!("3233")); // Example modulus

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::RSADecrypt),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "RSA decryption should succeed");
}

// ============================================================================
// HASHING TESTS
// ============================================================================

#[test]
fn test_sha256_empty() {
    let mut parameters = HashMap::new();
    parameters.insert("input".to_string(), json!(""));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::SHA256),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "SHA-256 of empty string should succeed");
}

#[test]
fn test_sha256_simple() {
    let mut parameters = HashMap::new();
    parameters.insert("input".to_string(), json!("Hello, World!"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::SHA256),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "SHA-256 of simple message should succeed");
}

#[test]
fn test_sha256_long() {
    let mut parameters = HashMap::new();
    parameters.insert(
        "input".to_string(),
        json!("The quick brown fox jumps over the lazy dog"),
    );

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::SHA256),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "SHA-256 of longer message should succeed");
}

#[test]
fn test_sha3_256_empty() {
    let mut parameters = HashMap::new();
    parameters.insert("input".to_string(), json!(""));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::SHA3_256),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "SHA3-256 of empty string should succeed");
}

#[test]
fn test_sha3_256_simple() {
    let mut parameters = HashMap::new();
    parameters.insert("input".to_string(), json!("Hello, World!"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::SHA3_256),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "SHA3-256 of simple message should succeed");
}

// ============================================================================
// CHINESE REMAINDER THEOREM TESTS
// ============================================================================

#[test]
fn test_chinese_remainder_two_congruences() {
    let mut parameters = HashMap::new();
    parameters.insert("remainders".to_string(), json!(["2", "3"]));
    parameters.insert("moduli".to_string(), json!(["3", "5"]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ChineseRemainder),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Chinese Remainder Theorem with 2 congruences should succeed"
    );
}

#[test]
fn test_chinese_remainder_three_congruences() {
    let mut parameters = HashMap::new();
    parameters.insert("remainders".to_string(), json!(["2", "3", "2"]));
    parameters.insert("moduli".to_string(), json!(["3", "5", "7"]));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ChineseRemainder),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(
        response.is_ok(),
        "Chinese Remainder Theorem with 3 congruences should succeed"
    );
}

// ============================================================================
// ELLIPTIC CURVE TESTS
// ============================================================================

#[test]
fn test_ec_point_add() {
    let mut parameters = HashMap::new();
    parameters.insert("curve_a".to_string(), json!("0"));
    parameters.insert("curve_b".to_string(), json!("7"));
    parameters.insert("p".to_string(), json!("17"));
    parameters.insert("x1".to_string(), json!("5"));
    parameters.insert("y1".to_string(), json!("1"));
    parameters.insert("x2".to_string(), json!("6"));
    parameters.insert("y2".to_string(), json!("3"));

    let request = ToolRequest::Compute(ComputeInput {
        operation: ComputeOp::NumberTheory(NumberTheoryOp::ECPointAdd),
        data: json!({}),
        parameters,
    });

    let dispatcher = computational_engine::implementations::create_default_dispatcher();
    let response = dispatcher.dispatch(request);
    assert!(response.is_ok(), "EC point addition should succeed");
}
