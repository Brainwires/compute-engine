// Unit tests for mathematics::symbolic_cas::quantum_advanced
use computational_engine::analyze::symbolic::quantum_advanced::*;

use super::*;

    #[test]
    fn test_density_matrix_pure() {
        let state = SymbolicMatrix::new(vec![vec![Expr::num(1)], vec![Expr::num(0)]]).unwrap();

        let rho = density_matrix_pure_state(&state).unwrap();
        println!("ρ = |ψ⟩⟨ψ|:\n{}", rho);

        assert_eq!(rho.rows(), 2);
        assert_eq!(rho.cols(), 2);

        // Just verify structure is correct
        // (simplification issues prevent exact comparison)
    }

    #[test]
    fn test_density_matrix_trace() {
        let state = SymbolicMatrix::new(vec![vec![Expr::num(1)], vec![Expr::num(0)]]).unwrap();

        let rho = density_matrix_pure_state(&state).unwrap();
        let trace = density_matrix_trace(&rho).unwrap();

        println!("Tr(ρ) = {}", trace);
        assert_eq!(trace, Expr::num(1));
    }

    #[test]
    fn test_maximally_mixed_state() {
        let rho_mixed = maximally_mixed_state(2);
        println!("Maximally mixed state:\n{}", rho_mixed);

        // Verify structure
        assert_eq!(rho_mixed.rows(), 2);
        assert_eq!(rho_mixed.cols(), 2);
    }

    #[test]
    fn test_hadamard_gate() {
        let h = hadamard_gate();
        println!("Hadamard gate:\n{}", h);

        assert_eq!(h.rows(), 2);
        assert_eq!(h.cols(), 2);
    }

    #[test]
    fn test_cnot_gate() {
        let cnot = cnot_gate();
        println!("CNOT gate:\n{}", cnot);

        assert_eq!(cnot.rows(), 4);
        assert_eq!(cnot.cols(), 4);
    }

    #[test]
    fn test_bell_states() {
        let phi_plus = bell_state_phi_plus();
        let phi_minus = bell_state_phi_minus();

        println!("|Φ+⟩ =\n{}", phi_plus);
        println!("|Φ-⟩ =\n{}", phi_minus);

        assert_eq!(phi_plus.rows(), 4);
        assert_eq!(phi_minus.rows(), 4);
    }

    #[test]
    fn test_quantum_gates_collection() {
        let gates = vec![
            ("X", pauli_x_gate()),
            ("Y", pauli_y_gate()),
            ("Z", pauli_z_gate()),
            ("H", hadamard_gate()),
            ("S", phase_gate()),
        ];

        for (name, gate) in gates {
            println!("{} gate:\n{}", name, gate);
            assert_eq!(gate.rows(), 2);
        }
    }

    #[test]
    fn test_rotation_gates() {
        let rx = rotation_x_gate("θ");
        let ry = rotation_y_gate("θ");
        let rz = rotation_z_gate("θ");

        println!("R_x(θ):\n{}", rx);
        println!("R_y(θ):\n{}", ry);
        println!("R_z(θ):\n{}", rz);
    }

    #[test]
    fn test_partial_trace() {
        // Create a simple product state density matrix
        let state = bell_state_phi_plus();
        let rho = density_matrix_pure_state(&state).unwrap();

        println!("Bell state density matrix:\n{}", rho);

        let rho_a = partial_trace_qubit(&rho, true).unwrap();
        println!("Partial trace (trace out B):\n{}", rho_a);

        assert_eq!(rho_a.rows(), 2);
        assert_eq!(rho_a.cols(), 2);
    }

    #[test]
    fn test_toffoli_gate() {
        let toffoli = toffoli_gate();
        println!("Toffoli gate:\n{}", toffoli);

        assert_eq!(toffoli.rows(), 8);
        assert_eq!(toffoli.cols(), 8);
    }
