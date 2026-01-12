// Unit tests for specialized::cryptographic_mathematics::lib
use computational_engine::compute::cryptographic_mathematics::lib::*;

use super::*;
    use num_bigint::BigInt;
    use num_traits::{One, Zero};

    #[test]
    fn test_mod_exp_basic() {
        let base = BigInt::from(2);
        let exp = BigInt::from(10);
        let modulus = BigInt::from(1000);
        let result = mod_exp(&base, &exp, &modulus);
        assert_eq!(result, BigInt::from(24));
    }

    #[test]
    fn test_extended_gcd() {
        let a = BigInt::from(240);
        let b = BigInt::from(46);
        let (gcd, x, y) = extended_gcd(&a, &b);
        assert_eq!(gcd, BigInt::from(2));
        assert_eq!(&a * &x + &b * &y, gcd);
    }

    #[test]
    fn test_mod_inverse() {
        let a = BigInt::from(3);
        let m = BigInt::from(11);
        let inv = mod_inverse(&a, &m).unwrap();
        assert_eq!(inv, BigInt::from(4));
        assert_eq!((&a * &inv) % &m, BigInt::one());
    }

    #[test]
    fn test_chinese_remainder_theorem() {
        let remainders = vec![BigInt::from(2), BigInt::from(3), BigInt::from(2)];
        let moduli = vec![BigInt::from(3), BigInt::from(5), BigInt::from(7)];
        let result = chinese_remainder_theorem(&remainders, &moduli).unwrap();
        assert_eq!(result, BigInt::from(23));
    }

    #[test]
    fn test_miller_rabin_primes() {
        assert!(miller_rabin_test(&BigInt::from(2), 10));
        assert!(miller_rabin_test(&BigInt::from(3), 10));
        assert!(miller_rabin_test(&BigInt::from(7), 10));
        assert!(miller_rabin_test(&BigInt::from(97), 10));
    }

    #[test]
    fn test_miller_rabin_composites() {
        assert!(!miller_rabin_test(&BigInt::from(4), 10));
        assert!(!miller_rabin_test(&BigInt::from(9), 10));
        assert!(!miller_rabin_test(&BigInt::from(100), 10));
    }

    #[test]
    fn test_rsa_encrypt_decrypt() {
        let (n, e, d) = generate_rsa_keypair(512);
        let message = BigInt::from(42);
        let ciphertext = rsa_encrypt(&message, &e, &n);
        let decrypted = rsa_decrypt(&ciphertext, &d, &n);
        assert_eq!(decrypted, message);
    }

    #[test]
    fn test_sha256_known_hash() {
        let hash = sha256("hello world");
        assert_eq!(
            hash,
            "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
        );
    }

    #[test]
    fn test_sha3_256_known_hash() {
        let hash = sha3_256("hello world");
        assert_eq!(
            hash,
            "644bcc7e564373040999aac89e7622f3ca71fba1d972fd94a31c3bfbf24e3938"
        );
    }

    #[test]
    fn test_discrete_log_small() {
        let base = BigInt::from(2);
        let target = BigInt::from(8);
        let modulus = BigInt::from(11);
        let result = discrete_log_bsgs(&base, &target, &modulus, 10).unwrap();
        assert_eq!(result, BigInt::from(3));
    }
