//! Unit tests for specialized::game_theory

use crate::solve::specialized::game_theory::normal_form::*;

#[test]
fn test_prisoners_dilemma() {
    let game = prisoners_dilemma();

    assert_eq!(game.num_players, 2);
}

#[test]
fn test_battle_of_sexes() {
    let game = battle_of_sexes();

    assert_eq!(game.num_players, 2);
}
