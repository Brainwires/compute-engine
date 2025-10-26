// Unit tests for api::mod
use computational_engine::api::mod::*;

use super::*;
    use serde_json::json;

    #[test]
    fn test_process_json_request() {
        let json_req = json!({
            "module": "linear_algebra",
            "operation": "compute_qr",
            "parameters": {
                "matrix": [[1.0, 2.0], [3.0, 4.0]]
            }
        });

        let response_str = process_json_request(&json_req.to_string());
        assert!(!response_str.is_empty());
    }
