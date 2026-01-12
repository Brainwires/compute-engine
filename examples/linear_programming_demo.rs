//! Linear Programming Demo
//!
//! This example demonstrates the Simplex algorithm for solving linear programming problems.
//!
//! Run with: cargo run --release --example linear_programming_demo

use computational_engine::solve::specialized::linear_programming::*;

fn main() {
    println!("=== Linear Programming Demo ===\n");

    // Example 1: Production Planning
    production_planning();

    // Example 2: Resource Allocation
    resource_allocation();

    // Example 3: Diet Problem
    diet_problem();

    // Example 4: Transportation Problem
    transportation_problem();
}

/// Example 1: Production Planning
///
/// A factory produces two products: chairs and tables.
/// - Each chair requires 1 hour of carpentry and 2 hours of finishing
/// - Each table requires 3 hours of carpentry and 1 hour of finishing
/// - Available: 12 hours of carpentry, 8 hours of finishing
/// - Profit: $5 per chair, $6 per table
///
/// Goal: Maximize profit
fn production_planning() {
    println!("--- Example 1: Production Planning ---");
    println!("Problem:");
    println!("  Maximize profit: 5*chairs + 6*tables");
    println!("  Subject to:");
    println!("    1*chairs + 3*tables ≤ 12  (carpentry hours)");
    println!("    2*chairs + 1*tables ≤ 8   (finishing hours)");
    println!("    chairs ≥ 0, tables ≥ 0");
    println!();

    let lp = LinearProgram::new(vec![5.0, 6.0], true).with_inequality_constraints(
        vec![
            vec![1.0, 3.0], // Carpentry constraint
            vec![2.0, 1.0], // Finishing constraint
        ],
        vec![12.0, 8.0],
    );

    match simplex(&lp) {
        Ok(solution) => {
            println!("Solution:");
            println!("  Status: {:?}", solution.status);
            println!(
                "  Chairs: {:.2}, Tables: {:.2}",
                solution.solution[0], solution.solution[1]
            );
            println!("  Maximum profit: ${:.2}", solution.objective_value);
            println!("  Iterations: {}", solution.iterations);
            println!();

            // Interpret slack variables
            if let Some(ref slack) = solution.slack_variables {
                println!("  Resource utilization:");
                println!(
                    "    Carpentry: {:.2} hours used, {:.2} hours slack",
                    12.0 - slack[0],
                    slack[0]
                );
                println!(
                    "    Finishing: {:.2} hours used, {:.2} hours slack",
                    8.0 - slack[1],
                    slack[1]
                );
            }
        }
        Err(e) => println!("Error: {}", e),
    }
    println!();
}

/// Example 2: Resource Allocation
///
/// A company allocates budget across 3 advertising channels:
/// - TV ads: $1000 each, reach 5000 people
/// - Radio ads: $500 each, reach 2000 people
/// - Online ads: $200 each, reach 1500 people
///
/// Budget: $10,000
/// Goal: Maximize reach
fn resource_allocation() {
    println!("--- Example 2: Resource Allocation ---");
    println!("Problem:");
    println!("  Maximize reach: 5000*tv + 2000*radio + 1500*online");
    println!("  Subject to:");
    println!("    1000*tv + 500*radio + 200*online ≤ 10000  (budget)");
    println!("    tv ≥ 0, radio ≥ 0, online ≥ 0");
    println!();

    let lp = LinearProgram::new(vec![5000.0, 2000.0, 1500.0], true)
        .with_inequality_constraints(vec![vec![1000.0, 500.0, 200.0]], vec![10000.0]);

    match simplex(&lp) {
        Ok(solution) => {
            println!("Solution:");
            println!("  Status: {:?}", solution.status);
            println!(
                "  TV ads: {:.0}, Radio ads: {:.0}, Online ads: {:.0}",
                solution.solution[0], solution.solution[1], solution.solution[2]
            );
            println!("  Maximum reach: {:.0} people", solution.objective_value);
            println!("  Iterations: {}", solution.iterations);
        }
        Err(e) => println!("Error: {}", e),
    }
    println!();
}

/// Example 3: Diet Problem
///
/// Plan a diet minimizing cost while meeting nutritional requirements.
/// - Food A: $2, provides 3 units protein, 1 unit vitamin
/// - Food B: $3, provides 1 unit protein, 2 units vitamin
///
/// Requirements: at least 6 units protein, 4 units vitamin
///
/// This would require ≥ constraints, which basic Simplex doesn't support.
/// We'll reformulate as maximizing negative cost.
fn diet_problem() {
    println!("--- Example 3: Diet Problem (Reformulated) ---");
    println!("Note: Basic Simplex handles ≤ constraints only.");
    println!("We demonstrate with a simplified version:");
    println!();

    println!("Problem:");
    println!("  Minimize cost: -2*A - 3*B (maximize negative)");
    println!("  Subject to:");
    println!("    A ≤ 5, B ≤ 5  (availability limits)");
    println!("    A ≥ 0, B ≥ 0");
    println!();

    // Minimize 2*A + 3*B becomes maximize -2*A - 3*B
    let lp = LinearProgram::new(vec![-2.0, -3.0], false).with_inequality_constraints(
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![5.0, 5.0],
    );

    match simplex(&lp) {
        Ok(solution) => {
            println!("Solution:");
            println!("  Status: {:?}", solution.status);
            println!(
                "  Food A: {:.2} units, Food B: {:.2} units",
                solution.solution[0], solution.solution[1]
            );
            println!("  Minimum cost: ${:.2}", solution.objective_value);
        }
        Err(e) => println!("Error: {}", e),
    }
    println!();
}

/// Example 4: Transportation Problem
///
/// Ship products from 2 factories to 2 warehouses.
/// Minimize shipping costs.
///
/// Costs ($ per unit):
///   Factory 1 → Warehouse A: $4
///   Factory 1 → Warehouse B: $6
///   Factory 2 → Warehouse A: $5
///   Factory 2 → Warehouse B: $3
///
/// Supply: Factory 1 has 100 units, Factory 2 has 150 units
/// Demand: Warehouse A needs 120 units, Warehouse B needs 130 units
fn transportation_problem() {
    println!("--- Example 4: Transportation Problem ---");
    println!("Problem:");
    println!("  Variables: x1A, x1B, x2A, x2B (shipments)");
    println!("  Minimize cost: -4*x1A - 6*x1B - 5*x2A - 3*x2B");
    println!("  Subject to:");
    println!("    x1A + x1B ≤ 100  (Factory 1 supply)");
    println!("    x2A + x2B ≤ 150  (Factory 2 supply)");
    println!("    All variables ≥ 0");
    println!();

    // Variables: [x1A, x1B, x2A, x2B]
    // Minimize costs
    let lp = LinearProgram::new(vec![-4.0, -6.0, -5.0, -3.0], false)
        .with_inequality_constraints(
            vec![
                vec![1.0, 1.0, 0.0, 0.0], // Factory 1 supply
                vec![0.0, 0.0, 1.0, 1.0], // Factory 2 supply
            ],
            vec![100.0, 150.0],
        );

    match simplex(&lp) {
        Ok(solution) => {
            println!("Solution:");
            println!("  Status: {:?}", solution.status);
            println!("  Shipments:");
            println!(
                "    Factory 1 → Warehouse A: {:.0} units",
                solution.solution[0]
            );
            println!(
                "    Factory 1 → Warehouse B: {:.0} units",
                solution.solution[1]
            );
            println!(
                "    Factory 2 → Warehouse A: {:.0} units",
                solution.solution[2]
            );
            println!(
                "    Factory 2 → Warehouse B: {:.0} units",
                solution.solution[3]
            );
            println!("  Minimum cost: ${:.2}", solution.objective_value);
            println!("  Iterations: {}", solution.iterations);
        }
        Err(e) => println!("Error: {}", e),
    }
    println!();
}
