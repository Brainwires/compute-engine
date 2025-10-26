// Unit tests for api_legacy
use computational_engine::api_legacy::*;

use super::*;

    #[test]
    fn test_list_operations() {
        let ops = list_all_operations();
        assert!(ops.contains_key("tensor_calculus"));
        assert!(ops.contains_key("advanced_calculus"));

        let total: usize = ops.values().map(|v| v.len()).sum();
        eprintln!("Total operations listed: {}", total);
        assert!(total >= 110, "Should have at least 110 operations");
    }

    #[test]
    fn test_json_request() {
        let request_json = r#"{
            "module": "tensor_calculus",
            "operation": "solve_vacuum_einstein",
            "parameters": {
                "system": "spherical"
            }
        }"#;

        let response = process_json_request(request_json);
        assert!(!response.is_empty());
    }
