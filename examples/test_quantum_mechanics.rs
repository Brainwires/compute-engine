use computational_engine::analyze::symbolic::{
    Expr, SymbolicMatrix, angular_momentum_x, angular_momentum_y, angular_momentum_z,
    annihilation_operator_symbolic, anticommutator, commutator, creation_operator_symbolic,
    dirac_gamma_0, dirac_gamma_1, dirac_gamma_2, dirac_gamma_3, expectation_value, pauli_x,
    pauli_y, pauli_z, time_evolution_operator,
};

fn main() {
    println!("=== Quantum Mechanics Operations Demo ===\n");

    // Test 1: Pauli Matrices
    println!("Test 1: Pauli Spin Matrices");
    let sigma_x = pauli_x();
    let sigma_y = pauli_y();
    let sigma_z = pauli_z();

    println!("σ_x (Pauli X):");
    println!("{}\n", sigma_x);

    println!("σ_y (Pauli Y):");
    println!("{}\n", sigma_y);

    println!("σ_z (Pauli Z):");
    println!("{}\n", sigma_z);

    println!("---\n");

    // Test 2: Pauli Matrix Properties
    println!("Test 2: Pauli Matrix Properties");

    // σ_x² = I
    let sigma_x_squared = sigma_x.mul(&sigma_x).unwrap();
    println!("σ_x² (should be identity):");
    println!("{}\n", sigma_x_squared);

    // [σ_x, σ_y] = 2i σ_z
    let comm_xy = commutator(&sigma_x, &sigma_y).unwrap();
    println!("[σ_x, σ_y] = 2iσ_z:");
    println!("{}\n", comm_xy);

    // {σ_x, σ_y} = 0
    let anticomm_xy = anticommutator(&sigma_x, &sigma_y).unwrap();
    println!("{{σ_x, σ_y}} = 0 (anticommutator):");
    println!("{}\n", anticomm_xy);

    println!("---\n");

    // Test 3: Dirac Gamma Matrices
    println!("Test 3: Dirac Gamma Matrices (for Dirac equation)");
    let gamma_0 = dirac_gamma_0();
    let gamma_1 = dirac_gamma_1();

    println!("γ^0:");
    println!("{}\n", gamma_0);

    println!("γ^1:");
    println!("{}\n", gamma_1);

    // {γ^μ, γ^ν} = 2g^μν I
    let anticomm_gamma = anticommutator(&gamma_0, &gamma_0).unwrap();
    println!("{{γ^0, γ^0}} (should be 2I):");
    println!("{}\n", anticomm_gamma);

    println!("---\n");

    // Test 4: Angular Momentum Operators
    println!("Test 4: Angular Momentum Operators");
    let l_x = angular_momentum_x(None);
    let l_y = angular_momentum_y(None);
    let l_z = angular_momentum_z(None);

    println!("L_x = (ℏ/2)σ_x:");
    println!("{}\n", l_x);

    println!("L_y = (ℏ/2)σ_y:");
    println!("{}\n", l_y);

    println!("L_z = (ℏ/2)σ_z:");
    println!("{}\n", l_z);

    // [L_x, L_y] = iℏ L_z
    let comm_l = commutator(&l_x, &l_y).unwrap();
    println!("[L_x, L_y] = iℏL_z:");
    println!("{}\n", comm_l);

    println!("---\n");

    // Test 5: Ladder Operators
    println!("Test 5: Harmonic Oscillator Ladder Operators");
    let a_dagger = creation_operator_symbolic();
    let a = annihilation_operator_symbolic();

    println!("Creation operator a†:");
    println!("{}\n", a_dagger);

    println!("Annihilation operator a:");
    println!("{}\n", a);

    println!("---\n");

    // Test 6: Time Evolution
    println!("Test 6: Time Evolution Operator");
    let u_t = time_evolution_operator("H", "t");
    println!("U(t) = exp(-iHt/ℏ):");
    println!("{}\n", u_t);

    println!("---\n");

    // Test 7: Expectation Values
    println!("Test 7: Expectation Values");

    // Spin-up state |↑> = [1, 0]
    let spin_up = SymbolicMatrix::new(vec![vec![Expr::num(1)], vec![Expr::num(0)]]).unwrap();

    // Spin-down state |↓> = [0, 1]
    let spin_down = SymbolicMatrix::new(vec![vec![Expr::num(0)], vec![Expr::num(1)]]).unwrap();

    println!("Spin-up state |↑> = [1, 0]ᵀ");
    println!("Spin-down state |↓> = [0, 1]ᵀ\n");

    // <↑|σ_z|↑>
    let exp_val_up_z = expectation_value(&spin_up, &sigma_z).unwrap();
    println!("<↑|σ_z|↑> = {} (should be 1)", exp_val_up_z);

    // <↓|σ_z|↓>
    let exp_val_down_z = expectation_value(&spin_down, &sigma_z).unwrap();
    println!("<↓|σ_z|↓> = {} (should be -1)", exp_val_down_z);

    // <↑|σ_x|↑>
    let exp_val_up_x = expectation_value(&spin_up, &sigma_x).unwrap();
    println!("<↑|σ_x|↑> = {} (should be 0)", exp_val_up_x);

    println!("\n---\n");

    // Test 8: Commutator Relations
    println!("Test 8: Fundamental Commutation Relations");

    // Create position and momentum operators symbolically
    println!("Position operator: x");
    println!("Momentum operator: p");
    println!("Canonical commutation: [x, p] = iℏ");
    println!("(Symbolic representation - actual computation requires operator algebra)\n");

    println!("---\n");

    println!("✅ Quantum mechanics operations complete!\n");
    println!("Key capabilities demonstrated:");
    println!("• Pauli matrices (σ_x, σ_y, σ_z)");
    println!("• Dirac gamma matrices (γ^μ)");
    println!("• Commutator and anticommutator algebra");
    println!("• Angular momentum operators (L_x, L_y, L_z)");
    println!("• Ladder operators (a†, a)");
    println!("• Time evolution operator U(t)");
    println!("• Expectation value calculations");
    println!("• Spin state operations");
}
