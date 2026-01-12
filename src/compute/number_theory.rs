//! Number theory operations
//!
//! Provides prime generation, modular arithmetic, RSA cryptography, and related operations

use crate::engine::*;
use crate::compute::cryptographic_mathematics::*;
use num_bigint::BigInt;
use std::str::FromStr;

/// Compute number theory operations
pub fn compute_number_theory(
    op: &NumberTheoryOp,
    input: &ComputeInput,
) -> ToolResult<ComputeOutput> {
    let result_json = match op {
        NumberTheoryOp::GeneratePrime => {
            let bits: u32 = input
                .parameters
                .get("bit_length")
                .and_then(|v| v.as_u64())
                .map(|v| v as u32)
                .ok_or("Missing or invalid bit_length parameter")?;

            let prime = generate_prime(bits);
            serde_json::json!({
                "prime": prime.to_string(),
                "bit_length": bits
            })
        }
        NumberTheoryOp::ModExp => {
            let base = input
                .parameters
                .get("base")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid base parameter")?;
            let exp = input
                .parameters
                .get("exponent")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid exponent parameter")?;
            let modulus = input
                .parameters
                .get("modulus")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid modulus parameter")?;

            let result = mod_exp(&base, &exp, &modulus);
            serde_json::json!({
                "result": result.to_string(),
                "operation": "modular_exponentiation"
            })
        }
        NumberTheoryOp::ModInv => {
            let a = input
                .parameters
                .get("a")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid a parameter")?;
            let m = input
                .parameters
                .get("modulus")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid modulus parameter")?;

            let result = mod_inverse(&a, &m).ok_or("Modular inverse does not exist")?;
            serde_json::json!({
                "inverse": result.to_string(),
                "operation": "modular_inverse"
            })
        }
        NumberTheoryOp::GCD => {
            let a = input
                .parameters
                .get("a")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid a parameter")?;
            let b = input
                .parameters
                .get("b")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid b parameter")?;

            let (gcd, x, y) = extended_gcd(&a, &b);
            serde_json::json!({
                "gcd": gcd.to_string(),
                "bezout_x": x.to_string(),
                "bezout_y": y.to_string(),
                "operation": "extended_gcd"
            })
        }
        NumberTheoryOp::LCM => {
            let a = input
                .parameters
                .get("a")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid a parameter")?;
            let b = input
                .parameters
                .get("b")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid b parameter")?;

            let (gcd, _, _) = extended_gcd(&a, &b);
            let lcm = (&a * &b) / gcd;
            serde_json::json!({
                "lcm": lcm.to_string(),
                "operation": "lcm"
            })
        }
        NumberTheoryOp::RSAKeypair => {
            let bits: u32 = input
                .parameters
                .get("bit_length")
                .and_then(|v| v.as_u64())
                .map(|v| v as u32)
                .unwrap_or(2048);

            if bits < 512 || bits > 4096 {
                return Err("RSA key size must be between 512 and 4096 bits".to_string());
            }

            let (n, e, d) = generate_rsa_keypair(bits);
            serde_json::json!({
                "public_key": {
                    "n": n.to_string(),
                    "e": e.to_string()
                },
                "private_key": {
                    "n": n.to_string(),
                    "d": d.to_string()
                }
            })
        }
        NumberTheoryOp::RSAEncrypt => {
            let message = input
                .parameters
                .get("message")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid message parameter")?;
            let e = input
                .parameters
                .get("e")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid e parameter")?;
            let n = input
                .parameters
                .get("n")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid n parameter")?;

            let ciphertext = rsa_encrypt(&message, &e, &n);
            serde_json::json!({
                "ciphertext": ciphertext.to_string()
            })
        }
        NumberTheoryOp::RSADecrypt => {
            let ciphertext = input
                .parameters
                .get("ciphertext")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid ciphertext parameter")?;
            let d = input
                .parameters
                .get("d")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid d parameter")?;
            let n = input
                .parameters
                .get("n")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid n parameter")?;

            let plaintext = rsa_decrypt(&ciphertext, &d, &n);
            serde_json::json!({
                "plaintext": plaintext.to_string()
            })
        }
        NumberTheoryOp::SHA256 => {
            let input_str = input
                .parameters
                .get("input")
                .and_then(|v| v.as_str())
                .ok_or("Missing input parameter")?;

            let hash = sha256(input_str);
            serde_json::json!({
                "hash": hash
            })
        }
        NumberTheoryOp::SHA3_256 => {
            let input_str = input
                .parameters
                .get("input")
                .and_then(|v| v.as_str())
                .ok_or("Missing input parameter")?;

            let hash = sha3_256(input_str);
            serde_json::json!({
                "hash": hash
            })
        }
        NumberTheoryOp::ChineseRemainder => {
            let remainders_arr = input
                .parameters
                .get("remainders")
                .and_then(|v| v.as_array())
                .ok_or("Missing or invalid remainders parameter")?;
            let moduli_arr = input
                .parameters
                .get("moduli")
                .and_then(|v| v.as_array())
                .ok_or("Missing or invalid moduli parameter")?;

            let remainders: Result<Vec<BigInt>, _> = remainders_arr
                .iter()
                .map(|v| {
                    v.as_str()
                        .and_then(|s| BigInt::from_str(s).ok())
                        .ok_or("Invalid remainder")
                })
                .collect();
            let moduli: Result<Vec<BigInt>, _> = moduli_arr
                .iter()
                .map(|v| {
                    v.as_str()
                        .and_then(|s| BigInt::from_str(s).ok())
                        .ok_or("Invalid modulus")
                })
                .collect();

            let remainders = remainders?;
            let moduli = moduli?;

            let result = chinese_remainder_theorem(&remainders, &moduli)
                .ok_or("No solution exists for given system")?;

            serde_json::json!({
                "result": result.to_string()
            })
        }
        NumberTheoryOp::DiscreteLog => {
            let base = input
                .parameters
                .get("base")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid base parameter")?;
            let target = input
                .parameters
                .get("target")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid target parameter")?;
            let modulus = input
                .parameters
                .get("modulus")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid modulus parameter")?;
            let max_exp: u32 = input
                .parameters
                .get("max_exponent")
                .and_then(|v| v.as_u64())
                .map(|v| v as u32)
                .unwrap_or(1000000);

            let result = discrete_log_bsgs(&base, &target, &modulus, max_exp)
                .ok_or("Discrete logarithm not found within max_exponent bound")?;

            serde_json::json!({
                "exponent": result.to_string()
            })
        }
        NumberTheoryOp::PrimalityTest => {
            let n = input
                .parameters
                .get("n")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid n parameter")?;
            let rounds: u32 = input
                .parameters
                .get("rounds")
                .and_then(|v| v.as_u64())
                .map(|v| v as u32)
                .unwrap_or(10);

            let is_prime = miller_rabin_test(&n, rounds);
            let confidence = 1.0 - (0.25_f64).powi(rounds as i32);

            serde_json::json!({
                "is_prime": is_prime,
                "confidence": format!("{:.6}", confidence)
            })
        }

        NumberTheoryOp::EulerTotient => {
            let n = input
                .parameters
                .get("n")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid n parameter")?;

            let result = euler_totient(&n);
            serde_json::json!({
                "totient": result.to_string()
            })
        }
        NumberTheoryOp::CarmichaelLambda => {
            let n = input
                .parameters
                .get("n")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid n parameter")?;

            let result = carmichael_lambda(&n);
            serde_json::json!({
                "lambda": result.to_string()
            })
        }
        NumberTheoryOp::ECPointAdd => {
            let x1 = input
                .parameters
                .get("x1")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid x1 parameter")?;
            let y1 = input
                .parameters
                .get("y1")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid y1 parameter")?;
            let x2 = input
                .parameters
                .get("x2")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid x2 parameter")?;
            let y2 = input
                .parameters
                .get("y2")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid y2 parameter")?;
            let curve_a = input
                .parameters
                .get("curve_a")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid curve_a parameter")?;
            let p = input
                .parameters
                .get("p")
                .and_then(|v| v.as_str())
                .and_then(|s| BigInt::from_str(s).ok())
                .ok_or("Missing or invalid p parameter")?;

            let (x3_opt, y3_opt) = elliptic_curve_point_add(&x1, &y1, &x2, &y2, &curve_a, &p);

            match (x3_opt, y3_opt) {
                (Some(x3), Some(y3)) => serde_json::json!({
                    "x3": x3.to_string(),
                    "y3": y3.to_string()
                }),
                _ => serde_json::json!({
                    "point": "Infinity"
                }),
            }
        }
    };

    Ok(ComputeOutput {
        result: result_json,
        additional: None,
        metadata: None,
    })
}
