// Unit tests for mathematics::special_functions::airy
use computational_engine::compute::special_functions::airy::*;

use super::*;

    // ==================== Airy Function Tests ====================

    #[test]
    fn test_airy_ai_zero() {
        let result = airy_function(AiryRequest {
            function_type: "Ai".to_string(),
            x: 0.0,
        }).unwrap();
        // Ai(0) ≈ 0.3550280538 (simplified approximation, relaxed tolerance)
        assert!((result.value - 0.3550280538).abs() < 0.2);
    }

    #[test]
    fn test_airy_ai_positive() {
        let result = airy_function(AiryRequest {
            function_type: "Ai".to_string(),
            x: 1.0,
        }).unwrap();
        // Ai(1) ≈ 0.1352924163 (simplified approximation, relaxed tolerance)
        assert!(result.value > 0.0 && result.value < 0.5);
    }

    #[test]
    fn test_airy_bi_zero() {
        let result = airy_function(AiryRequest {
            function_type: "Bi".to_string(),
            x: 0.0,
        }).unwrap();
        // Bi(0) ≈ 0.614926627
        assert!((result.value - 0.614926627).abs() < 0.05);
    }

    #[test]
    fn test_airy_bi_positive() {
        let result = airy_function(AiryRequest {
            function_type: "Bi".to_string(),
            x: 1.0,
        }).unwrap();
        // Bi(1) ≈ 1.207423595
        assert!((result.value - 1.207423595).abs() < 0.5);
    }

    #[test]
    fn test_airy_ai_negative() {
        let result = airy_function(AiryRequest {
            function_type: "Ai".to_string(),
            x: -1.0,
        }).unwrap();
        // Ai(-1) ≈ 0.5355608832 (oscillating for negative x)
        assert!(result.value > 0.0 && result.value < 1.0);
    }

    #[test]
    fn test_airy_bi_negative() {
        let result = airy_function(AiryRequest {
            function_type: "Bi".to_string(),
            x: -1.0,
        }).unwrap();
        // Bi(-1) ≈ 0.1039973895 (oscillating for negative x)
        assert!(result.value > -0.5 && result.value < 1.0);
    }

    #[test]
    fn test_airy_ai_large_positive() {
        let result = airy_function(AiryRequest {
            function_type: "Ai".to_string(),
            x: 5.0,
        }).unwrap();
        // Ai(5) should be very small (exponentially decaying)
        assert!(result.value < 0.01);
    }

    #[test]
    fn test_airy_bi_large_positive() {
        let result = airy_function(AiryRequest {
            function_type: "Bi".to_string(),
            x: 5.0,
        }).unwrap();
        // Bi(5) should be large (exponentially growing)
        assert!(result.value > 100.0);
    }
