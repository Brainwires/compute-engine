use num_bigint::{BigInt, BigUint, RandBigInt};
use num_traits::{Num, One, Signed, Zero};
use rand::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest as Sha2Digest, Sha256};
use sha3::Sha3_256;
use std::collections::HashMap;

// Legacy request/response structs for backwards compatibility
#[derive(Debug, Serialize, Deserialize)]
#[allow(dead_code)]
struct CryptographicRequest {
    operation: String,
    parameters: HashMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
#[allow(dead_code)]
struct CryptographicResult {
    success: bool,
    result: Option<serde_json::Value>,
    error: Option<String>,
}

// Modular exponentiation: base^exp mod modulus
pub fn mod_exp(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    if modulus == &BigInt::one() {
        return BigInt::zero();
    }

    let mut result = BigInt::one();
    let mut base = base % modulus;
    let mut exp = exp.clone();

    while exp > BigInt::zero() {
        if &exp % 2 == BigInt::one() {
            result = (&result * &base) % modulus;
        }
        exp >>= 1;
        base = (&base * &base) % modulus;
    }

    result
}

// Extended Euclidean Algorithm
pub fn extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if b == &BigInt::zero() {
        return (a.clone(), BigInt::one(), BigInt::zero());
    }

    let (gcd, x1, y1) = extended_gcd(b, &(a % b));
    let x = y1.clone();
    let y = x1 - (a / b) * &y1;

    (gcd, x, y)
}

// Modular multiplicative inverse
pub fn mod_inverse(a: &BigInt, m: &BigInt) -> Option<BigInt> {
    let (gcd, x, _) = extended_gcd(a, m);
    if gcd != BigInt::one() {
        None
    } else {
        Some(((x % m) + m) % m)
    }
}

// Chinese Remainder Theorem
pub fn chinese_remainder_theorem(remainders: &[BigInt], moduli: &[BigInt]) -> Option<BigInt> {
    if remainders.len() != moduli.len() || remainders.is_empty() {
        return None;
    }

    let product = moduli.iter().fold(BigInt::one(), |acc, m| acc * m);
    let mut result = BigInt::zero();

    for (i, remainder) in remainders.iter().enumerate() {
        let mi = &product / &moduli[i];
        let yi = mod_inverse(&mi, &moduli[i])?;
        result += remainder * &mi * &yi;
    }

    Some(result % &product)
}

// Miller-Rabin primality test
pub fn miller_rabin_test(n: &BigInt, k: u32) -> bool {
    if n <= &BigInt::one() {
        return false;
    }
    if n <= &BigInt::from(3) {
        return n == &BigInt::from(2) || n == &BigInt::from(3);
    }
    if n % 2 == BigInt::zero() {
        return false;
    }

    // Write n-1 as d * 2^r
    let n_minus_1: BigInt = n - 1;
    let mut r = 0u32;
    let mut d = n_minus_1.clone();

    while &d % 2 == BigInt::zero() {
        d >>= 1;
        r += 1;
    }

    let mut rng = thread_rng();

    for _ in 0..k {
        let a = rng.gen_bigint_range(&BigInt::from(2), &(&n_minus_1 + 1));
        let mut x = mod_exp(&a, &d, n);

        if x == BigInt::one() || x == n_minus_1 {
            continue;
        }

        let mut composite = true;
        for _ in 0..r - 1 {
            x = mod_exp(&x, &BigInt::from(2), n);
            if x == n_minus_1 {
                composite = false;
                break;
            }
        }

        if composite {
            return false;
        }
    }

    true
}

// Generate large prime number
pub fn generate_prime(bits: u32) -> BigInt {
    let mut rng = thread_rng();

    loop {
        let mut candidate = rng.gen_biguint(bits as u64);

        // Ensure odd number
        if &candidate % 2u32 == BigUint::zero() {
            candidate += 1u32;
        }

        let candidate_bigint = BigInt::from(candidate);
        if miller_rabin_test(&candidate_bigint, 10) {
            return candidate_bigint;
        }
    }
}

// RSA key generation
pub fn generate_rsa_keypair(bits: u32) -> (BigInt, BigInt, BigInt) {
    let half_bits = bits / 2;
    let p = generate_prime(half_bits);
    let q = generate_prime(half_bits);

    let n = &p * &q;
    let phi_n = (&p - 1) * (&q - 1);

    let e = BigInt::from(65537); // Common choice for public exponent
    let d = mod_inverse(&e, &phi_n).expect("Failed to compute private exponent");

    (n, e, d)
}

// RSA encryption
pub fn rsa_encrypt(message: &BigInt, e: &BigInt, n: &BigInt) -> BigInt {
    mod_exp(message, e, n)
}

// RSA decryption
pub fn rsa_decrypt(ciphertext: &BigInt, d: &BigInt, n: &BigInt) -> BigInt {
    mod_exp(ciphertext, d, n)
}

// Discrete logarithm (baby-step giant-step for small groups)
pub fn discrete_log_bsgs(
    base: &BigInt,
    target: &BigInt,
    modulus: &BigInt,
    max_exp: u32,
) -> Option<BigInt> {
    let m = ((max_exp as f64).sqrt().ceil() as u32) + 1;
    let mut baby_steps = HashMap::new();

    // Baby steps
    let mut gamma = BigInt::one();
    for j in 0..m {
        if gamma == *target {
            return Some(BigInt::from(j));
        }
        baby_steps.insert(gamma.clone(), j);
        gamma = (&gamma * base) % modulus;
    }

    // Giant steps
    let base_inv_m = mod_inverse(&mod_exp(base, &BigInt::from(m), modulus), modulus)?;
    let mut y = target.clone();

    for i in 0..m {
        if let Some(&j) = baby_steps.get(&y) {
            let result = BigInt::from(i) * BigInt::from(m) + BigInt::from(j);
            if result < BigInt::from(max_exp) {
                return Some(result);
            }
        }
        y = (&y * &base_inv_m) % modulus;
    }

    None
}

// Elliptic curve point addition (simplified for y^2 = x^3 + ax + b over Fp)
#[derive(Debug, Clone, PartialEq)]
enum ECPoint {
    Infinity,
    Point(BigInt, BigInt),
}

fn ec_add(p1: &ECPoint, p2: &ECPoint, a: &BigInt, modulus: &BigInt) -> ECPoint {
    match (p1, p2) {
        (ECPoint::Infinity, p) | (p, ECPoint::Infinity) => p.clone(),
        (ECPoint::Point(x1, y1), ECPoint::Point(x2, y2)) => {
            if x1 == x2 {
                if y1 == y2 {
                    // Point doubling
                    let numerator = (3 * x1 * x1 + a) % modulus;
                    let denominator = (2 * y1) % modulus;
                    let slope =
                        (&numerator * &mod_inverse(&denominator, modulus).unwrap()) % modulus;

                    let x3 = (&slope * &slope - 2 * x1) % modulus;
                    let x3 = ((x3 % modulus) + modulus) % modulus;

                    let y3 = (&slope * (x1 - &x3) - y1) % modulus;
                    let y3 = ((y3 % modulus) + modulus) % modulus;

                    ECPoint::Point(x3, y3)
                } else {
                    ECPoint::Infinity
                }
            } else {
                // Point addition
                let numerator = (y2 - y1) % modulus;
                let numerator = ((numerator % modulus) + modulus) % modulus;

                let denominator = (x2 - x1) % modulus;
                let denominator = ((denominator % modulus) + modulus) % modulus;

                let slope = (&numerator * &mod_inverse(&denominator, modulus).unwrap()) % modulus;

                let x3 = (&slope * &slope - x1 - x2) % modulus;
                let x3 = ((x3 % modulus) + modulus) % modulus;

                let y3 = (&slope * (x1 - &x3) - y1) % modulus;
                let y3 = ((y3 % modulus) + modulus) % modulus;

                ECPoint::Point(x3, y3)
            }
        }
    }
}

// Euler's totient function φ(n)
// Counts integers from 1 to n that are coprime with n
pub fn euler_totient(n: &BigInt) -> BigInt {
    if n <= &BigInt::one() {
        return BigInt::one();
    }

    let mut result = n.clone();
    let mut n_copy = n.clone();
    let mut p = BigInt::from(2);

    // Factor n and apply formula: φ(n) = n * ∏(1 - 1/p) for all prime factors p
    while &p * &p <= n_copy {
        if &n_copy % &p == BigInt::zero() {
            // Remove all factors of p
            while &n_copy % &p == BigInt::zero() {
                n_copy /= &p;
            }
            // Apply φ(p^k) = p^(k-1) * (p-1)
            result -= &result / &p;
        }
        p += 1;
    }

    // If n_copy > 1, then it's a prime factor
    if n_copy > BigInt::one() {
        result -= &result / &n_copy;
    }

    result
}

// Carmichael's lambda function λ(n)
// Returns the smallest positive integer m such that a^m ≡ 1 (mod n) for all a coprime to n
pub fn carmichael_lambda(n: &BigInt) -> BigInt {
    if n <= &BigInt::one() {
        return BigInt::one();
    }

    let mut result = BigInt::one();
    let mut n_copy = n.clone();
    let mut p = BigInt::from(2);

    // Factor n and compute λ for each prime power
    while &p * &p <= n_copy {
        if &n_copy % &p == BigInt::zero() {
            let mut pk = p.clone();
            n_copy /= &p;

            // Count the exponent k
            while &n_copy % &p == BigInt::zero() {
                pk *= &p;
                n_copy /= &p;
            }

            // λ(p^k) = φ(p^k) for odd primes, λ(2^k) = 2^(k-2) for k >= 3
            let lambda_pk = if &p == &BigInt::from(2) && &pk >= &BigInt::from(8) {
                &pk / 4
            } else {
                euler_totient(&pk)
            };

            // λ(n) = lcm of all λ(p^k)
            result = lcm(&result, &lambda_pk);
        }
        p += 1;
    }

    // If n_copy > 1, then it's a prime
    if n_copy > BigInt::one() {
        let lambda_p = &n_copy - 1;
        result = lcm(&result, &lambda_p);
    }

    result
}

// Least common multiple
fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    if a == &BigInt::zero() || b == &BigInt::zero() {
        return BigInt::zero();
    }
    let (gcd, _, _) = extended_gcd(a, b);
    (a * b).abs() / gcd
}

// Public wrapper for elliptic curve point addition
pub fn elliptic_curve_point_add(
    x1: &BigInt,
    y1: &BigInt,
    x2: &BigInt,
    y2: &BigInt,
    a: &BigInt,
    modulus: &BigInt,
) -> (Option<BigInt>, Option<BigInt>) {
    let p1 = ECPoint::Point(x1.clone(), y1.clone());
    let p2 = ECPoint::Point(x2.clone(), y2.clone());

    match ec_add(&p1, &p2, a, modulus) {
        ECPoint::Infinity => (None, None),
        ECPoint::Point(x, y) => (Some(x), Some(y)),
    }
}

// SHA256 hashing
pub fn sha256(input: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(input.as_bytes());
    let result = hasher.finalize();
    hex::encode(result)
}

// SHA3-256 hashing
pub fn sha3_256(input: &str) -> String {
    let mut hasher = Sha3_256::new();
    hasher.update(input.as_bytes());
    let result = hasher.finalize();
    hex::encode(result)
}

// Legacy request processor for backwards compatibility
#[allow(dead_code)]
fn process_request(input: &str) -> String {
    let request: CryptographicRequest = match serde_json::from_str(input) {
        Ok(req) => req,
        Err(e) => {
            return serde_json::to_string(&CryptographicResult {
                success: false,
                result: None,
                error: Some(format!("Invalid JSON: {}", e)),
            })
            .unwrap();
        }
    };

    let result = match request.operation.as_str() {
        "modular_exponentiation" => {
            let base_str = request
                .parameters
                .get("base")
                .and_then(|v| v.as_str())
                .unwrap_or("0");
            let exp_str = request
                .parameters
                .get("exponent")
                .and_then(|v| v.as_str())
                .unwrap_or("0");
            let mod_str = request
                .parameters
                .get("modulus")
                .and_then(|v| v.as_str())
                .unwrap_or("1");

            match (
                BigInt::from_str_radix(base_str, 10),
                BigInt::from_str_radix(exp_str, 10),
                BigInt::from_str_radix(mod_str, 10),
            ) {
                (Ok(base), Ok(exp), Ok(modulus)) => {
                    let result = mod_exp(&base, &exp, &modulus);
                    CryptographicResult {
                        success: true,
                        result: Some(serde_json::json!({ "result": result.to_string() })),
                        error: None,
                    }
                }
                _ => CryptographicResult {
                    success: false,
                    result: None,
                    error: Some("Invalid number format".to_string()),
                },
            }
        }

        "extended_gcd" => {
            let a_str = request
                .parameters
                .get("a")
                .and_then(|v| v.as_str())
                .unwrap_or("0");
            let b_str = request
                .parameters
                .get("b")
                .and_then(|v| v.as_str())
                .unwrap_or("0");

            match (
                BigInt::from_str_radix(a_str, 10),
                BigInt::from_str_radix(b_str, 10),
            ) {
                (Ok(a), Ok(b)) => {
                    let (gcd, x, y) = extended_gcd(&a, &b);
                    CryptographicResult {
                        success: true,
                        result: Some(serde_json::json!({
                            "gcd": gcd.to_string(),
                            "x": x.to_string(),
                            "y": y.to_string()
                        })),
                        error: None,
                    }
                }
                _ => CryptographicResult {
                    success: false,
                    result: None,
                    error: Some("Invalid number format".to_string()),
                },
            }
        }

        "chinese_remainder_theorem" => {
            let remainders_json = request
                .parameters
                .get("remainders")
                .cloned()
                .unwrap_or_default();
            let moduli_json = request
                .parameters
                .get("moduli")
                .cloned()
                .unwrap_or_default();

            let remainders: Result<Vec<BigInt>, _> =
                serde_json::from_value::<Vec<String>>(remainders_json)
                    .unwrap_or_default()
                    .iter()
                    .map(|s| BigInt::from_str_radix(s, 10))
                    .collect();

            let moduli: Result<Vec<BigInt>, _> = serde_json::from_value::<Vec<String>>(moduli_json)
                .unwrap_or_default()
                .iter()
                .map(|s| BigInt::from_str_radix(s, 10))
                .collect();

            match (remainders, moduli) {
                (Ok(remainders), Ok(moduli)) => {
                    match chinese_remainder_theorem(&remainders, &moduli) {
                        Some(result) => CryptographicResult {
                            success: true,
                            result: Some(serde_json::json!({ "result": result.to_string() })),
                            error: None,
                        },
                        None => CryptographicResult {
                            success: false,
                            result: None,
                            error: Some("No solution exists".to_string()),
                        },
                    }
                }
                _ => CryptographicResult {
                    success: false,
                    result: None,
                    error: Some("Invalid number format".to_string()),
                },
            }
        }

        "generate_prime" => {
            let bits = request
                .parameters
                .get("bits")
                .and_then(|v| v.as_u64())
                .unwrap_or(1024) as u32;
            if bits < 8 || bits > 4096 {
                CryptographicResult {
                    success: false,
                    result: None,
                    error: Some("Bits must be between 8 and 4096".to_string()),
                }
            } else {
                let prime = generate_prime(bits);
                CryptographicResult {
                    success: true,
                    result: Some(serde_json::json!({ "prime": prime.to_string() })),
                    error: None,
                }
            }
        }

        "generate_rsa_keypair" => {
            let bits = request
                .parameters
                .get("bits")
                .and_then(|v| v.as_u64())
                .unwrap_or(2048) as u32;
            if bits < 512 || bits > 4096 {
                CryptographicResult {
                    success: false,
                    result: None,
                    error: Some("RSA key size must be between 512 and 4096 bits".to_string()),
                }
            } else {
                let (n, e, d) = generate_rsa_keypair(bits);
                CryptographicResult {
                    success: true,
                    result: Some(serde_json::json!({
                        "public_key": {
                            "n": n.to_string(),
                            "e": e.to_string()
                        },
                        "private_key": {
                            "n": n.to_string(),
                            "d": d.to_string()
                        }
                    })),
                    error: None,
                }
            }
        }

        "primality_test" => {
            let n_str = request
                .parameters
                .get("n")
                .and_then(|v| v.as_str())
                .unwrap_or("0");
            let rounds = request
                .parameters
                .get("rounds")
                .and_then(|v| v.as_u64())
                .unwrap_or(10) as u32;

            match BigInt::from_str_radix(n_str, 10) {
                Ok(n) => {
                    let is_prime = miller_rabin_test(&n, rounds);
                    CryptographicResult {
                        success: true,
                        result: Some(serde_json::json!({
                            "is_prime": is_prime,
                            "confidence": format!("{:.6}", 1.0 - (0.25_f64).powi(rounds as i32))
                        })),
                        error: None,
                    }
                }
                _ => CryptographicResult {
                    success: false,
                    result: None,
                    error: Some("Invalid number format".to_string()),
                },
            }
        }

        "sha256" => {
            let input = request
                .parameters
                .get("input")
                .and_then(|v| v.as_str())
                .unwrap_or("");
            let hash = sha256(input);
            CryptographicResult {
                success: true,
                result: Some(serde_json::json!({ "hash": hash })),
                error: None,
            }
        }

        "sha3_256" => {
            let input = request
                .parameters
                .get("input")
                .and_then(|v| v.as_str())
                .unwrap_or("");
            let hash = sha3_256(input);
            CryptographicResult {
                success: true,
                result: Some(serde_json::json!({ "hash": hash })),
                error: None,
            }
        }

        _ => CryptographicResult {
            success: false,
            result: None,
            error: Some(format!("Unknown operation: {}", request.operation)),
        },
    };

    serde_json::to_string(&result).unwrap()
}

// Test module
#[cfg(test)]
#[path = "../../../tests/unit/specialized/cryptographic_mathematics_tests.rs"]
mod tests;

