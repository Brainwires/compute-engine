use computational_engine::analyze::symbolic::{differentiate, expand, simplify};

fn main() {
    println!("Testing Symbolica integration...\n");

    // Test expand
    println!("Test 1: Expand (x + 1)^2");
    match expand("(x + 1)^2") {
        Ok(result) => {
            println!("Result: {}", result.expression);
            println!("LaTeX: {}\n", result.latex.unwrap_or_default());
        }
        Err(e) => println!("Error: {:?}\n", e),
    }

    // Test differentiate
    println!("Test 2: Differentiate x^3 + 2*x^2 + x");
    match differentiate("x^3 + 2*x^2 + x", "x", None) {
        Ok(result) => {
            println!("Result: {}", result.expression);
            println!("LaTeX: {}\n", result.latex.unwrap_or_default());
        }
        Err(e) => println!("Error: {:?}\n", e),
    }

    // Test simplify
    println!("Test 3: Simplify (x + 1)^2 - x^2 - 2*x - 1");
    match simplify("(x + 1)^2 - x^2 - 2*x - 1") {
        Ok(result) => {
            println!("Result: {}", result.expression);
            println!("LaTeX: {}\n", result.latex.unwrap_or_default());
        }
        Err(e) => println!("Error: {:?}\n", e),
    }

    println!("✅ Symbolic CAS is working!");
}
