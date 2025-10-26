// Unit tests for specialized::game_theory::extensive_form
use computational_engine::specialized::game_theory::*;

use super::*;

    #[test]
    fn test_extensive_game_creation() {
        let game = ExtensiveGame::new();
        assert_eq!(game.nodes.len(), 0);
    }
