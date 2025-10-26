// Unit tests for specialized::game_theory::evolutionary
use computational_engine::specialized::game_theory::evolutionary::*;

use super::*;

    #[test]
    fn test_evolutionary_game() {
        // Hawk-Dove game
        let payoffs = vec![
            vec![0.0, 3.0],  // Hawk vs Hawk, Hawk vs Dove
            vec![1.0, 2.0],  // Dove vs Hawk, Dove vs Dove
        ];

        let mut game = EvolutionaryGame::new(payoffs, vec![0.5, 0.5]);

        // Simulate dynamics
        for _ in 0..10 {
            game.step(0.1);
        }

        // Population should sum to 1
        let total: f64 = game.population.iter().sum();
        assert!((total - 1.0).abs() < 1e-6);
    }
