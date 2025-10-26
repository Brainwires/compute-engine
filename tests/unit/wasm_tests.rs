// Unit tests for wasm
use computational_engine::wasm::*;

use super::*;

    #[test]
    fn test_engine_creation() {
        let engine = ComputationalEngine::new();
        assert_eq!(engine.version(), crate::VERSION);
    }

    #[test]
    fn test_json_processing() {
        let engine = ComputationalEngine::new();
        let request = r#"{"tool":"solve","input":{"equations":["x^2 - 4 = 0"]}}"#;
        let response = engine.process_json(request);
        assert!(!response.is_empty());
    }
