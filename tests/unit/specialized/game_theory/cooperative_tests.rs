// Unit tests for specialized::game_theory::cooperative
use computational_engine::solve::specialized::game_theory::*;

use super::*;

    #[test]
    fn test_cooperative_game() {
        let mut game = CooperativeGame::new(3);
        game.set_value(&[0], 0.0);
        game.set_value(&[1], 0.0);
        game.set_value(&[2], 0.0);
        game.set_value(&[0, 1, 2], 12.0);

        let shapley = game.shapley_value();
        assert_eq!(shapley.len(), 3);
    }
