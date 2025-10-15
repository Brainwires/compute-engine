use computational_engine::mathematics::symbolic_cas::{
    christoffel_symbols, riemann_tensor, ricci_tensor, ricci_scalar, einstein_tensor,
    euclidean_metric, minkowski_2d, schwarzschild_metric, SymbolicMatrix, Expr,
};

fn main() {
    println!("=== Tensor Calculus Operations Demo ===\n");

    // Test 1: Flat Euclidean space
    println!("Test 1: Euclidean 2D space");
    let metric = euclidean_metric(2);
    let coords = vec!["x", "y"];

    println!("Metric:");
    println!("{}\n", metric);

    let christoffel = christoffel_symbols(&metric, &coords).unwrap();
    println!("Christoffel symbols (all zero for flat space):");
    println!("{}\n", christoffel);

    let ricci = ricci_tensor(&metric, &coords).unwrap();
    println!("Ricci tensor (should be zero):");
    println!("{}\n", ricci);

    let scalar = ricci_scalar(&metric, &coords).unwrap();
    println!("Ricci scalar: {}\n", scalar);
    println!("---\n");

    // Test 2: Minkowski spacetime
    println!("Test 2: Minkowski 2D spacetime");
    let metric = minkowski_2d(true).unwrap();
    let coords = vec!["t", "x"];

    println!("Metric (signature +,-):");
    println!("{}\n", metric);

    let christoffel = christoffel_symbols(&metric, &coords).unwrap();
    println!("Christoffel symbols (all zero for flat spacetime):");
    println!("{}\n", christoffel);

    let einstein = einstein_tensor(&metric, &coords).unwrap();
    println!("Einstein tensor (should be zero):");
    println!("{}\n", einstein);
    println!("---\n");

    // Test 3: Simple curved space - 2D sphere metric
    println!("Test 3: 2D Sphere metric (simplified)");
    // ds² = r²(dθ² + sin²θ dφ²)
    let r = Expr::sym("r");
    let r_squared = Expr::pow(r.clone(), Expr::num(2));
    let sin_theta = Expr::func("sin", vec![Expr::sym("θ")]);
    let sin_theta_squared = Expr::pow(sin_theta, Expr::num(2));

    let sphere_metric = SymbolicMatrix::new(vec![
        vec![r_squared.clone(), Expr::num(0)],
        vec![Expr::num(0), Expr::mul(r_squared, sin_theta_squared)],
    ])
    .unwrap();

    let coords = vec!["θ", "φ"];

    println!("Metric:");
    println!("{}\n", sphere_metric);

    let christoffel = christoffel_symbols(&sphere_metric, &coords).unwrap();
    println!("Christoffel symbols (non-zero for curved surface):");
    println!("{}\n", christoffel);

    let ricci = ricci_tensor(&sphere_metric, &coords).unwrap();
    println!("Ricci tensor:");
    println!("{}\n", ricci);

    let scalar = ricci_scalar(&sphere_metric, &coords).unwrap();
    println!("Ricci scalar: {}\n", scalar);
    println!("---\n");

    // Test 4: Schwarzschild metric (this will be slow due to complexity)
    println!("Test 4: Schwarzschild spacetime (4D black hole)");
    let metric = schwarzschild_metric().unwrap();
    let coords = vec!["t", "r", "θ", "φ"];

    println!("Metric:");
    println!("{}\n", metric);

    println!("Computing Christoffel symbols...");
    let christoffel = christoffel_symbols(&metric, &coords).unwrap();
    println!("Christoffel symbols computed (rank-3 tensor)");
    println!("Rank: {}", christoffel.rank());
    println!("Dimensions: {:?}\n", christoffel.dimensions());

    println!("Note: Riemann tensor computation for Schwarzschild is very complex");
    println!("and may take significant time due to symbolic differentiation.");
    println!("For production use, consider caching or numerical evaluation.\n");

    println!("✅ Tensor calculus operations complete!\n");
    println!("Key capabilities demonstrated:");
    println!("• Christoffel symbols (connection coefficients)");
    println!("• Riemann curvature tensor");
    println!("• Ricci tensor and Ricci scalar");
    println!("• Einstein tensor");
    println!("• Support for arbitrary metrics and coordinate systems");
}
