//! Expression tree representation for mathematical functions

use rand::Rng;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Operator {
    Add,
    Subtract,
    Multiply,
    Divide,
    Power,
    Sin,
    Cos,
    Exp,
    Log,
    Abs,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExpressionNode {
    Constant(f64),
    Variable(String), // e.g., "x"
    Binary {
        op: Operator,
        left: Box<ExpressionNode>,
        right: Box<ExpressionNode>,
    },
    Unary {
        op: Operator,
        child: Box<ExpressionNode>,
    },
}

impl ExpressionNode {
    /// Evaluate the expression tree given variable values
    pub fn evaluate(&self, variables: &std::collections::HashMap<String, f64>) -> Result<f64, String> {
        match self {
            ExpressionNode::Constant(val) => Ok(*val),
            ExpressionNode::Variable(name) => {
                variables.get(name)
                    .copied()
                    .ok_or_else(|| format!("Variable '{}' not found", name))
            }
            ExpressionNode::Binary { op, left, right } => {
                let left_val = left.evaluate(variables)?;
                let right_val = right.evaluate(variables)?;

                let result = match op {
                    Operator::Add => left_val + right_val,
                    Operator::Subtract => left_val - right_val,
                    Operator::Multiply => left_val * right_val,
                    Operator::Divide => {
                        if right_val.abs() < 1e-10 {
                            return Err("Division by zero".to_string());
                        }
                        left_val / right_val
                    }
                    Operator::Power => {
                        if left_val < 0.0 && right_val.fract() != 0.0 {
                            return Err("Negative number to fractional power".to_string());
                        }
                        left_val.powf(right_val)
                    }
                    _ => return Err(format!("Invalid binary operator: {:?}", op)),
                };

                if result.is_nan() || result.is_infinite() {
                    return Err("Result is NaN or infinite".to_string());
                }

                Ok(result)
            }
            ExpressionNode::Unary { op, child } => {
                let child_val = child.evaluate(variables)?;

                let result = match op {
                    Operator::Sin => child_val.sin(),
                    Operator::Cos => child_val.cos(),
                    Operator::Exp => {
                        if child_val > 100.0 {
                            return Err("Exp overflow".to_string());
                        }
                        child_val.exp()
                    }
                    Operator::Log => {
                        if child_val <= 0.0 {
                            return Err("Log of non-positive number".to_string());
                        }
                        child_val.ln()
                    }
                    Operator::Abs => child_val.abs(),
                    _ => return Err(format!("Invalid unary operator: {:?}", op)),
                };

                if result.is_nan() || result.is_infinite() {
                    return Err("Result is NaN or infinite".to_string());
                }

                Ok(result)
            }
        }
    }

    /// Count the number of nodes in the tree (complexity measure)
    pub fn complexity(&self) -> usize {
        match self {
            ExpressionNode::Constant(_) | ExpressionNode::Variable(_) => 1,
            ExpressionNode::Binary { left, right, .. } => {
                1 + left.complexity() + right.complexity()
            }
            ExpressionNode::Unary { child, .. } => 1 + child.complexity(),
        }
    }

    /// Generate a random expression tree
    pub fn random<R: Rng>(
        rng: &mut R,
        variables: &[String],
        max_depth: usize,
        current_depth: usize,
    ) -> Self {
        if current_depth >= max_depth || (current_depth > 0 && rng.gen_bool(0.3)) {
            // Leaf node
            if rng.gen_bool(0.5) && !variables.is_empty() {
                ExpressionNode::Variable(variables[rng.gen_range(0..variables.len())].clone())
            } else {
                ExpressionNode::Constant(rng.gen_range(-10.0..10.0))
            }
        } else {
            // Internal node
            if rng.gen_bool(0.7) {
                // Binary operator
                let op = match rng.gen_range(0..5) {
                    0 => Operator::Add,
                    1 => Operator::Subtract,
                    2 => Operator::Multiply,
                    3 => Operator::Divide,
                    _ => Operator::Power,
                };
                ExpressionNode::Binary {
                    op,
                    left: Box::new(Self::random(rng, variables, max_depth, current_depth + 1)),
                    right: Box::new(Self::random(rng, variables, max_depth, current_depth + 1)),
                }
            } else {
                // Unary operator
                let op = match rng.gen_range(0..5) {
                    0 => Operator::Sin,
                    1 => Operator::Cos,
                    2 => Operator::Exp,
                    3 => Operator::Log,
                    _ => Operator::Abs,
                };
                ExpressionNode::Unary {
                    op,
                    child: Box::new(Self::random(rng, variables, max_depth, current_depth + 1)),
                }
            }
        }
    }

    /// Convert to human-readable string
    pub fn to_string(&self) -> String {
        match self {
            ExpressionNode::Constant(val) => {
                if val.fract() == 0.0 && val.abs() < 1e6 {
                    format!("{}", *val as i64)
                } else {
                    format!("{:.4}", val)
                }
            }
            ExpressionNode::Variable(name) => name.clone(),
            ExpressionNode::Binary { op, left, right } => {
                let op_str = match op {
                    Operator::Add => "+",
                    Operator::Subtract => "-",
                    Operator::Multiply => "*",
                    Operator::Divide => "/",
                    Operator::Power => "^",
                    _ => "?",
                };
                format!("({} {} {})", left.to_string(), op_str, right.to_string())
            }
            ExpressionNode::Unary { op, child } => {
                let op_str = match op {
                    Operator::Sin => "sin",
                    Operator::Cos => "cos",
                    Operator::Exp => "exp",
                    Operator::Log => "log",
                    Operator::Abs => "abs",
                    _ => "?",
                };
                format!("{}({})", op_str, child.to_string())
            }
        }
    }
}
