//! Evolution and mutation rules for expression trees

use super::expression::*;
use rand::Rng;

impl ExpressionNode {
    /// Mutate the expression tree according to the four rules
    pub fn mutate<R: Rng>(&self, rng: &mut R, variables: &[String]) -> Self {
        match rng.gen_range(0..4) {
            0 => self.mutate_constant(rng),
            1 => self.mutate_replace_subtree(rng, variables),
            2 => self.mutate_delete_node(rng),
            _ => self.mutate_simplify_constants(),
        }
    }

    /// Rule 1: Increase or decrease a random constant
    fn mutate_constant<R: Rng>(&self, rng: &mut R) -> Self {
        match self {
            ExpressionNode::Constant(val) => {
                let delta = rng.gen_range(-1.0..1.0);
                ExpressionNode::Constant(val + delta)
            }
            ExpressionNode::Variable(name) => ExpressionNode::Variable(name.clone()),
            ExpressionNode::Binary { op, left, right } => {
                if rng.gen_bool(0.5) {
                    ExpressionNode::Binary {
                        op: op.clone(),
                        left: Box::new(left.mutate_constant(rng)),
                        right: right.clone(),
                    }
                } else {
                    ExpressionNode::Binary {
                        op: op.clone(),
                        left: left.clone(),
                        right: Box::new(right.mutate_constant(rng)),
                    }
                }
            }
            ExpressionNode::Unary { op, child } => ExpressionNode::Unary {
                op: op.clone(),
                child: Box::new(child.mutate_constant(rng)),
            },
        }
    }

    /// Rule 2: Replace a random node with a new random subtree
    fn mutate_replace_subtree<R: Rng>(&self, rng: &mut R, variables: &[String]) -> Self {
        if rng.gen_bool(0.2) {
            // Replace this node
            return ExpressionNode::random(rng, variables, 3, 0);
        }

        match self {
            ExpressionNode::Constant(_) | ExpressionNode::Variable(_) => self.clone(),
            ExpressionNode::Binary { op, left, right } => {
                if rng.gen_bool(0.5) {
                    ExpressionNode::Binary {
                        op: op.clone(),
                        left: Box::new(left.mutate_replace_subtree(rng, variables)),
                        right: right.clone(),
                    }
                } else {
                    ExpressionNode::Binary {
                        op: op.clone(),
                        left: left.clone(),
                        right: Box::new(right.mutate_replace_subtree(rng, variables)),
                    }
                }
            }
            ExpressionNode::Unary { op, child } => ExpressionNode::Unary {
                op: op.clone(),
                child: Box::new(child.mutate_replace_subtree(rng, variables)),
            },
        }
    }

    /// Rule 3: Delete a node by replacing it with one of its children
    fn mutate_delete_node<R: Rng>(&self, rng: &mut R) -> Self {
        if rng.gen_bool(0.2) {
            match self {
                ExpressionNode::Binary { left, right, .. } => {
                    return if rng.gen_bool(0.5) {
                        (**left).clone()
                    } else {
                        (**right).clone()
                    };
                }
                ExpressionNode::Unary { child, .. } => {
                    return (**child).clone();
                }
                _ => {}
            }
        }

        match self {
            ExpressionNode::Constant(_) | ExpressionNode::Variable(_) => self.clone(),
            ExpressionNode::Binary { op, left, right } => {
                if rng.gen_bool(0.5) {
                    ExpressionNode::Binary {
                        op: op.clone(),
                        left: Box::new(left.mutate_delete_node(rng)),
                        right: right.clone(),
                    }
                } else {
                    ExpressionNode::Binary {
                        op: op.clone(),
                        left: left.clone(),
                        right: Box::new(right.mutate_delete_node(rng)),
                    }
                }
            }
            ExpressionNode::Unary { op, child } => ExpressionNode::Unary {
                op: op.clone(),
                child: Box::new(child.mutate_delete_node(rng)),
            },
        }
    }

    /// Rule 4: Simplify nodes with only constants
    fn mutate_simplify_constants(&self) -> Self {
        match self {
            ExpressionNode::Constant(_) | ExpressionNode::Variable(_) => self.clone(),
            ExpressionNode::Binary { op, left, right } => {
                let left_simplified = left.mutate_simplify_constants();
                let right_simplified = right.mutate_simplify_constants();

                // Check if both children are constants
                if let (ExpressionNode::Constant(l), ExpressionNode::Constant(r)) =
                    (&left_simplified, &right_simplified)
                {
                    let vars = std::collections::HashMap::new();
                    let temp_expr = ExpressionNode::Binary {
                        op: op.clone(),
                        left: Box::new(ExpressionNode::Constant(*l)),
                        right: Box::new(ExpressionNode::Constant(*r)),
                    };
                    if let Ok(result) = temp_expr.evaluate(&vars) {
                        return ExpressionNode::Constant(result);
                    }
                }

                ExpressionNode::Binary {
                    op: op.clone(),
                    left: Box::new(left_simplified),
                    right: Box::new(right_simplified),
                }
            }
            ExpressionNode::Unary { op, child } => {
                let child_simplified = child.mutate_simplify_constants();

                // Check if child is a constant
                if let ExpressionNode::Constant(c) = child_simplified {
                    let vars = std::collections::HashMap::new();
                    let temp_expr = ExpressionNode::Unary {
                        op: op.clone(),
                        child: Box::new(ExpressionNode::Constant(c)),
                    };
                    if let Ok(result) = temp_expr.evaluate(&vars) {
                        return ExpressionNode::Constant(result);
                    }
                }

                ExpressionNode::Unary {
                    op: op.clone(),
                    child: Box::new(child_simplified),
                }
            }
        }
    }
}
