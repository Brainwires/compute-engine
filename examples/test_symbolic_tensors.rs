use computational_engine::mathematics::symbolic_cas::{
    Expr, IndexType, SymbolicMatrix, SymbolicTensor, euclidean_metric, minkowski_metric,
    schwarzschild_metric,
};

fn main() {
    println!("=== Symbolic Tensor Operations Demo ===\n");

    // Test 1: Scalar tensor
    println!("Test 1: Scalar (rank-0 tensor)");
    let scalar = SymbolicTensor::scalar(Expr::sym("φ"));
    println!("{}\n", scalar);

    // Test 2: Vector tensor
    println!("Test 2: Contravariant vector (rank-1 tensor)");
    let vector = SymbolicTensor::vector(
        vec![
            Expr::sym("v^0"),
            Expr::sym("v^1"),
            Expr::sym("v^2"),
            Expr::sym("v^3"),
        ],
        IndexType::Contravariant,
    )
    .unwrap();
    println!("{}\n", vector);

    // Test 3: Covariant vector
    println!("Test 3: Covariant vector");
    let covector = SymbolicTensor::vector(
        vec![Expr::sym("p_0"), Expr::sym("p_1"), Expr::sym("p_2")],
        IndexType::Covariant,
    )
    .unwrap();
    println!("{}\n", covector);

    // Test 4: Matrix as rank-2 tensor
    println!("Test 4: Rank-2 tensor (electromagnetic field tensor)");
    let f_matrix = SymbolicMatrix::new(vec![
        vec![
            Expr::num(0),
            Expr::sym("E_x"),
            Expr::sym("E_y"),
            Expr::sym("E_z"),
        ],
        vec![
            Expr::mul(Expr::num(-1), Expr::sym("E_x")),
            Expr::num(0),
            Expr::sym("B_z"),
            Expr::mul(Expr::num(-1), Expr::sym("B_y")),
        ],
        vec![
            Expr::mul(Expr::num(-1), Expr::sym("E_y")),
            Expr::mul(Expr::num(-1), Expr::sym("B_z")),
            Expr::num(0),
            Expr::sym("B_x"),
        ],
        vec![
            Expr::mul(Expr::num(-1), Expr::sym("E_z")),
            Expr::sym("B_y"),
            Expr::mul(Expr::num(-1), Expr::sym("B_x")),
            Expr::num(0),
        ],
    ])
    .unwrap();

    let f_tensor =
        SymbolicTensor::from_matrix(&f_matrix, [IndexType::Contravariant, IndexType::Covariant]);
    println!("Electromagnetic field tensor F^μ_ν:");
    println!("{}\n", f_tensor);

    // Test 5: Tensor contraction (trace)
    println!("Test 5: Tensor contraction (trace of 2×2 matrix)");
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::sym("T^0_0"), Expr::sym("T^0_1")],
        vec![Expr::sym("T^1_0"), Expr::sym("T^1_1")],
    ])
    .unwrap();

    let stress_tensor =
        SymbolicTensor::from_matrix(&mat, [IndexType::Contravariant, IndexType::Covariant]);
    let trace = stress_tensor.contract(0, 1).unwrap();
    println!("Trace of T^μ_ν:");
    println!("{}\n", trace);

    // Test 6: Outer product (tensor product)
    println!("Test 6: Outer product of two vectors");
    let u = SymbolicTensor::vector(
        vec![Expr::sym("u^0"), Expr::sym("u^1")],
        IndexType::Contravariant,
    )
    .unwrap();

    let v = SymbolicTensor::vector(
        vec![Expr::sym("v_0"), Expr::sym("v_1")],
        IndexType::Covariant,
    )
    .unwrap();

    let tensor_product = u.outer_product(&v).unwrap();
    println!("u^i ⊗ v_j:");
    println!("{}\n", tensor_product);

    // Test 7: Minkowski metric
    println!("Test 7: Minkowski metric tensor (special relativity)");
    let minkowski = minkowski_metric(true).unwrap();
    println!("η^μν (signature +,-,-,-):");
    println!("{}\n", minkowski);

    // Test 8: Raising and lowering indices
    println!("Test 8: Raising/lowering indices with metric");
    let minkowski_2d = SymbolicMatrix::new(vec![
        vec![Expr::num(1), Expr::num(0)],
        vec![Expr::num(0), Expr::num(-1)],
    ])
    .unwrap();

    let covariant_momentum =
        SymbolicTensor::vector(vec![Expr::sym("E"), Expr::sym("p_x")], IndexType::Covariant)
            .unwrap();

    let contravariant_momentum = covariant_momentum.raise_index(&minkowski_2d, 0).unwrap();
    println!("Original (covariant): p_μ");
    println!("{}", covariant_momentum);
    println!("Raised (contravariant): p^μ = η^μν p_ν");
    println!("{}\n", contravariant_momentum);

    // Test 9: Schwarzschild metric
    println!("Test 9: Schwarzschild metric (black hole spacetime)");
    let schwarzschild = schwarzschild_metric().unwrap();
    println!("Schwarzschild metric g_μν:");
    println!("{}\n", schwarzschild);

    // Test 10: Determinant of metric
    println!("Test 10: Determinant of metric tensor");
    let simple_metric = SymbolicMatrix::new(vec![
        vec![Expr::sym("g_00"), Expr::num(0)],
        vec![Expr::num(0), Expr::sym("g_11")],
    ])
    .unwrap();

    let det = simple_metric.determinant().unwrap();
    println!("det(g) for diagonal 2×2 metric:");
    println!("det(g) = {}\n", det);

    // Test 11: Euclidean metric
    println!("Test 11: Euclidean metric (flat space)");
    let euclidean = euclidean_metric(3);
    println!("δ_ij (3D Euclidean metric):");
    println!("{}\n", euclidean);

    println!("✅ Symbolic tensor operations complete!");
    println!("\nKey capabilities demonstrated:");
    println!("• Arbitrary rank tensors (scalars, vectors, matrices, higher)");
    println!("• Covariant/contravariant indices");
    println!("• Tensor contraction (trace)");
    println!("• Outer products (tensor products)");
    println!("• Raising/lowering indices with metric tensors");
    println!("• Physical metric tensors (Minkowski, Schwarzschild, etc.)");
}
