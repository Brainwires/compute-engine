// Unit tests for specialized::game_theory::normal_form
use computational_engine::solve::specialized::game_theory::*;

use super::*;

    #[test]
    fn test_prisoners_dilemma_nash() {
        let game = prisoners_dilemma();
        let equilibria = game.find_nash_equilibria();

        // (Defect, Defect) should be the unique Nash equilibrium
        assert_eq!(equilibria.len(), 1);
        assert_eq!(equilibria[0], vec![1, 1]);
    }

    #[test]
    fn test_prisoners_dilemma_dominant() {
        let game = prisoners_dilemma();
        let dominant = game.find_dominant_strategies();

        // Both players have Defect as dominant strategy
        assert_eq!(dominant[0], Some(1));
        assert_eq!(dominant[1], Some(1));
    }

    #[test]
    fn test_battle_of_sexes_nash() {
        let game = battle_of_sexes();
        let equilibria = game.find_nash_equilibria();

        // Two pure Nash equilibria: (Opera, Opera) and (Football, Football)
        assert_eq!(equilibria.len(), 2);
    }

    #[test]
    fn test_pareto_optimal() {
        let game = prisoners_dilemma();

        // (Cooperate, Cooperate) is Pareto optimal
        assert!(game.is_pareto_optimal(&vec![0, 0]));

        // (Defect, Defect) is not Pareto optimal
        assert!(!game.is_pareto_optimal(&vec![1, 1]));
    }

    #[test]
    fn test_game_creation() {
        let game = NormalFormGame::new(2, vec![2, 3]);
        assert_eq!(game.num_players, 2);
        assert_eq!(game.payoffs.len(), 6); // 2 * 3 profiles
    }
