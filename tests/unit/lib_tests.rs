// Unit tests for lib
use computational_engine::lib::*;

use super::*;

    #[test]
    fn test_version() {
        assert!(!VERSION.is_empty());
        assert_eq!(NAME, "brainwires-compute-engine");
    }
