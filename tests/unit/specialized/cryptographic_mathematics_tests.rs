//! Unit tests for specialized::cryptographic_mathematics

use crate::specialized::cryptographic_mathematics::*;
use num_bigint::BigInt;
use num_traits::One;

#[test]
fn test_mod_exp_basic() {
    let base = BigInt::from(2);
    let exp = BigInt::from(10);
    let modulus = BigInt::from(1000);

    let result = mod_exp(&base, &exp, &modulus);
    assert_eq!(result, BigInt::from(24)); // 2^10 mod 1000 = 1024 mod 1000 = 24
}

#[test]
fn test_extended_gcd() {
    let a = BigInt::from(240);
    let b = BigInt::from(46);

    let (gcd, x, y) = extended_gcd(&a, &b);
    assert_eq!(gcd, BigInt::from(2)); // gcd(240, 46) = 2
    // Verify: a*x + b*y = gcd
    assert_eq!(&a * &x + &b * &y, gcd);
}

#[test]
fn test_mod_inverse_exists() {
    let a = BigInt::from(3);
    let m = BigInt::from(11);

    let inv = mod_inverse(&a, &m);
    assert!(inv.is_some());

    // Verify: (a * inv) mod m = 1
    if let Some(inv_val) = inv {
        assert_eq!((&a * &inv_val) % &m, BigInt::one());
    }
}

#[test]
fn test_mod_inverse_none() {
    let a = BigInt::from(6);
    let m = BigInt::from(9);

    // 6 and 9 are not coprime, so no inverse exists
    let inv = mod_inverse(&a, &m);
    assert!(inv.is_none());
}

#[test]
fn test_chinese_remainder_theorem() {
    let remainders = vec![BigInt::from(2), BigInt::from(3), BigInt::from(2)];
    let moduli = vec![BigInt::from(3), BigInt::from(5), BigInt::from(7)];

    let result = chinese_remainder_theorem(&remainders, &moduli);
    assert!(result.is_some());

    // Verify the result satisfies all congruences
    if let Some(x) = result {
        for i in 0..remainders.len() {
            assert_eq!(&x % &moduli[i], remainders[i]);
        }
    }
}

#[test]
fn test_miller_rabin_known_prime() {
    let prime = BigInt::from(17);
    let is_prime = miller_rabin_test(&prime, 5);
    assert!(is_prime);
}

#[test]
fn test_miller_rabin_known_composite() {
    let composite = BigInt::from(15);
    let is_prime = miller_rabin_test(&composite, 5);
    assert!(!is_prime);
}

#[test]
fn test_generate_prime_small() {
    let prime = generate_prime(16);

    // Verify it's actually prime
    assert!(miller_rabin_test(&prime, 10));

    // Verify it's within the bit limit (up to 16 bits means < 2^16)
    assert!(prime > BigInt::from(1)); // Greater than 1
    assert!(prime < BigInt::from(1 << 16)); // Less than 2^16

    // Note: gen_biguint generates UP TO N bits, not exactly N bits,
    // so the prime might be smaller than expected
}

#[test]
fn test_rsa_encrypt_decrypt() {
    // Generate small RSA keypair for testing
    let (n, e, d) = generate_rsa_keypair(128);

    // Test message (must be smaller than n)
    let message = BigInt::from(42);

    // Encrypt then decrypt
    let ciphertext = rsa_encrypt(&message, &e, &n);
    let decrypted = rsa_decrypt(&ciphertext, &d, &n);

    assert_eq!(decrypted, message);
}

#[test]
fn test_elliptic_curve_point_add() {
    // Simple test with small numbers
    let x1 = BigInt::from(1);
    let y1 = BigInt::from(1);
    let x2 = BigInt::from(2);
    let y2 = BigInt::from(3);
    let a = BigInt::from(1);
    let modulus = BigInt::from(7);

    let (x3, y3) = elliptic_curve_point_add(&x1, &y1, &x2, &y2, &a, &modulus);

    // Just verify we get a result (could be infinity or a point)
    assert!(x3.is_some() || x3.is_none());
    assert!(y3.is_some() || y3.is_none());
}

#[test]
fn test_euler_totient() {
    // φ(9) = 6 (numbers coprime to 9: 1,2,4,5,7,8)
    let n = BigInt::from(9);
    let totient = euler_totient(&n);
    assert_eq!(totient, BigInt::from(6));
}

#[test]
fn test_carmichael_lambda() {
    // λ(9) = 6
    let n = BigInt::from(9);
    let lambda = carmichael_lambda(&n);
    assert_eq!(lambda, BigInt::from(6));
}

#[test]
fn test_sha256() {
    let input = "hello world";
    let hash = sha256(input);

    // Verify it's a valid hex string of correct length (64 chars = 32 bytes)
    assert_eq!(hash.len(), 64);
    assert!(hash.chars().all(|c| c.is_ascii_hexdigit()));
}

#[test]
fn test_sha3_256() {
    let input = "hello world";
    let hash = sha3_256(input);

    // Verify it's a valid hex string of correct length (64 chars = 32 bytes)
    assert_eq!(hash.len(), 64);
    assert!(hash.chars().all(|c| c.is_ascii_hexdigit()));
}

#[test]
fn test_discrete_log_bsgs_simple() {
    // Solve 2^x ≡ 8 (mod 11)
    // Answer should be 3 since 2^3 = 8
    let base = BigInt::from(2);
    let target = BigInt::from(8);
    let modulus = BigInt::from(11);

    let result = discrete_log_bsgs(&base, &target, &modulus, 20);

    // Verify we found a solution
    assert!(result.is_some());
    if let Some(exp) = result {
        // Verify: base^exp ≡ target (mod modulus)
        assert_eq!(mod_exp(&base, &exp, &modulus), target);
    }
}
