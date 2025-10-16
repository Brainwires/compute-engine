use computational_engine::mathematics::symbolic_cas::{
    Expr, SymbolicMatrix, characteristic_polynomial, eigenvalues_2x2, matrix_inverse,
};

fn main() {
    println!("=== Symbolic Matrix Operations Demo ===\n");

    // Test 1: Create a symbolic matrix
    println!("Test 1: Create a symbolic matrix");
    let mat = SymbolicMatrix::new(vec![
        vec![Expr::sym("a"), Expr::sym("b")],
        vec![Expr::sym("c"), Expr::sym("d")],
    ])
    .unwrap();
    println!("Matrix A:");
    println!("{}\n", mat);

    // Test 2: Matrix transpose
    println!("Test 2: Matrix transpose");
    let transposed = mat.transpose();
    println!("A^T:");
    println!("{}\n", transposed);

    // Test 3: Determinant
    println!("Test 3: Determinant");
    let det = mat.determinant().unwrap();
    println!("det(A) = {}\n", det);

    // Test 4: Trace
    println!("Test 4: Trace");
    let trace = mat.trace().unwrap();
    println!("tr(A) = {}\n", trace);

    // Test 5: Matrix multiplication
    println!("Test 5: Matrix multiplication (A × A)");
    let mat_squared = mat.mul(&mat).unwrap();
    println!("A²:");
    println!("{}\n", mat_squared);

    // Test 6: Characteristic polynomial
    println!("Test 6: Characteristic polynomial");
    let char_poly = characteristic_polynomial(&mat).unwrap();
    println!("det(A - λI) = {}\n", char_poly);

    // Test 7: Eigenvalues for 2x2 matrix
    println!("Test 7: Eigenvalues (symbolic)");
    let eigenvalues = eigenvalues_2x2(&mat).unwrap();
    println!("λ₁ = {}", eigenvalues[0]);
    println!("λ₂ = {}\n", eigenvalues[1]);

    // Test 8: Matrix inverse
    println!("Test 8: Matrix inverse");
    let inv = matrix_inverse(&mat).unwrap();
    println!("A⁻¹:");
    println!("{}\n", inv);

    // Test 9: Numeric example
    println!("Test 9: Numeric 2x2 matrix");
    let numeric_mat = SymbolicMatrix::new(vec![
        vec![Expr::num(3), Expr::num(1)],
        vec![Expr::num(1), Expr::num(3)],
    ])
    .unwrap();
    println!("Matrix B:");
    println!("{}", numeric_mat);

    let numeric_eigenvalues = eigenvalues_2x2(&numeric_mat).unwrap();
    println!("Eigenvalues:");
    println!("λ₁ = {}", numeric_eigenvalues[0]);
    println!("λ₂ = {}\n", numeric_eigenvalues[1]);

    // Test 10: Mixed symbolic and numeric
    println!("Test 10: Mixed symbolic and numeric matrix");
    let mixed_mat = SymbolicMatrix::new(vec![
        vec![Expr::num(2), Expr::sym("x")],
        vec![Expr::sym("x"), Expr::num(2)],
    ])
    .unwrap();
    println!("Matrix C:");
    println!("{}", mixed_mat);

    let mixed_det = mixed_mat.determinant().unwrap();
    println!("det(C) = {}\n", mixed_det);

    let mixed_eigenvalues = eigenvalues_2x2(&mixed_mat).unwrap();
    println!("Eigenvalues:");
    println!("λ₁ = {}", mixed_eigenvalues[0]);
    println!("λ₂ = {}\n", mixed_eigenvalues[1]);

    println!("✅ Symbolic matrix operations are working!");
}
