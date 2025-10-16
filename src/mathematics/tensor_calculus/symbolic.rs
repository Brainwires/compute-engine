use crate::TensorError;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum SymbolicExpr {
    Variable(String),
    Constant(f64),
    Add(Box<SymbolicExpr>, Box<SymbolicExpr>),
    Subtract(Box<SymbolicExpr>, Box<SymbolicExpr>),
    Multiply(Box<SymbolicExpr>, Box<SymbolicExpr>),
    Divide(Box<SymbolicExpr>, Box<SymbolicExpr>),
    Power(Box<SymbolicExpr>, Box<SymbolicExpr>),
    Function(String, Vec<SymbolicExpr>),
    Zero,
    One,
}

impl fmt::Display for SymbolicExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SymbolicExpr::Variable(var) => write!(f, "{}", var),
            SymbolicExpr::Constant(val) => {
                if val.fract() == 0.0 {
                    write!(f, "{}", *val as i64)
                } else {
                    write!(f, "{}", val)
                }
            }
            SymbolicExpr::Add(left, right) => write!(f, "({} + {})", left, right),
            SymbolicExpr::Subtract(left, right) => write!(f, "({} - {})", left, right),
            SymbolicExpr::Multiply(left, right) => write!(f, "({} * {})", left, right),
            SymbolicExpr::Divide(left, right) => write!(f, "({} / {})", left, right),
            SymbolicExpr::Power(base, exp) => write!(f, "{}^{}", base, exp),
            SymbolicExpr::Function(name, args) => {
                if args.is_empty() {
                    write!(f, "{}()", name)
                } else {
                    write!(
                        f,
                        "{}({})",
                        name,
                        args.iter()
                            .map(|arg| arg.to_string())
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                }
            }
            SymbolicExpr::Zero => write!(f, "0"),
            SymbolicExpr::One => write!(f, "1"),
        }
    }
}

impl SymbolicExpr {
    pub fn parse(input: &str) -> Result<Self, TensorError> {
        let input = input.trim();

        if input.is_empty() {
            return Ok(SymbolicExpr::Zero);
        }

        // Handle special cases first
        match input {
            "0" => return Ok(SymbolicExpr::Zero),
            "1" => return Ok(SymbolicExpr::One),
            _ => {}
        }

        // Handle constants
        if let Ok(val) = input.parse::<f64>() {
            return Ok(SymbolicExpr::Constant(val));
        }

        // Handle simple variables
        if input.chars().all(|c| c.is_alphanumeric() || c == '_') {
            return Ok(SymbolicExpr::Variable(input.to_string()));
        }

        // Handle powers (r^2, etc.)
        if let Some(captures) = Regex::new(r"^([a-zA-Z_][a-zA-Z0-9_]*)\^(\d+)$")
            .unwrap()
            .captures(input)
        {
            let base = SymbolicExpr::Variable(captures.get(1).unwrap().as_str().to_string());
            let exp =
                SymbolicExpr::Constant(captures.get(2).unwrap().as_str().parse::<f64>().unwrap());
            return Ok(SymbolicExpr::Power(Box::new(base), Box::new(exp)));
        }

        // Handle functions like sin(x), cos(theta), etc.
        if let Some(captures) = Regex::new(r"^([a-zA-Z_][a-zA-Z0-9_]*)\(([^)]*)\)$")
            .unwrap()
            .captures(input)
        {
            let func_name = captures.get(1).unwrap().as_str().to_string();
            let args_str = captures.get(2).unwrap().as_str();

            let mut args = Vec::new();
            if !args_str.is_empty() {
                for arg in args_str.split(',') {
                    args.push(Self::parse(arg.trim())?);
                }
            }

            return Ok(SymbolicExpr::Function(func_name, args));
        }

        // For now, treat complex expressions as variables
        // In a more complete implementation, we'd parse arithmetic expressions
        Ok(SymbolicExpr::Variable(input.to_string()))
    }

    pub fn simplify(&self) -> Self {
        match self {
            SymbolicExpr::Add(left, right) => {
                let left = left.simplify();
                let right = right.simplify();

                match (&left, &right) {
                    (SymbolicExpr::Zero, expr) | (expr, SymbolicExpr::Zero) => expr.clone(),
                    (SymbolicExpr::Constant(a), SymbolicExpr::Constant(b)) => {
                        SymbolicExpr::Constant(a + b)
                    }
                    _ => SymbolicExpr::Add(Box::new(left), Box::new(right)),
                }
            }
            SymbolicExpr::Subtract(left, right) => {
                let left = left.simplify();
                let right = right.simplify();

                match (&left, &right) {
                    (expr, SymbolicExpr::Zero) => expr.clone(),
                    (SymbolicExpr::Zero, expr) => {
                        SymbolicExpr::Subtract(Box::new(SymbolicExpr::Zero), Box::new(expr.clone()))
                    }
                    (SymbolicExpr::Constant(a), SymbolicExpr::Constant(b)) => {
                        SymbolicExpr::Constant(a - b)
                    }
                    _ => SymbolicExpr::Subtract(Box::new(left), Box::new(right)),
                }
            }
            SymbolicExpr::Multiply(left, right) => {
                let left = left.simplify();
                let right = right.simplify();

                match (&left, &right) {
                    (SymbolicExpr::Zero, _) | (_, SymbolicExpr::Zero) => SymbolicExpr::Zero,
                    (SymbolicExpr::One, expr) | (expr, SymbolicExpr::One) => expr.clone(),
                    (SymbolicExpr::Constant(a), SymbolicExpr::Constant(b)) => {
                        SymbolicExpr::Constant(a * b)
                    }
                    _ => SymbolicExpr::Multiply(Box::new(left), Box::new(right)),
                }
            }
            SymbolicExpr::Divide(left, right) => {
                let left = left.simplify();
                let right = right.simplify();

                match (&left, &right) {
                    (SymbolicExpr::Zero, _) => SymbolicExpr::Zero,
                    (expr, SymbolicExpr::One) => expr.clone(),
                    (SymbolicExpr::Constant(a), SymbolicExpr::Constant(b)) if *b != 0.0 => {
                        SymbolicExpr::Constant(a / b)
                    }
                    _ => SymbolicExpr::Divide(Box::new(left), Box::new(right)),
                }
            }
            SymbolicExpr::Power(base, exp) => {
                let base = base.simplify();
                let exp = exp.simplify();

                match (&base, &exp) {
                    (_, SymbolicExpr::Zero) => SymbolicExpr::One,
                    (expr, SymbolicExpr::One) => expr.clone(),
                    (SymbolicExpr::Zero, _) => SymbolicExpr::Zero,
                    (SymbolicExpr::One, _) => SymbolicExpr::One,
                    (SymbolicExpr::Constant(a), SymbolicExpr::Constant(b)) => {
                        SymbolicExpr::Constant(a.powf(*b))
                    }
                    _ => SymbolicExpr::Power(Box::new(base), Box::new(exp)),
                }
            }
            _ => self.clone(),
        }
    }

    pub fn derivative(&self, var: &str) -> Self {
        match self {
            SymbolicExpr::Variable(v) => {
                if v == var {
                    SymbolicExpr::One
                } else {
                    SymbolicExpr::Zero
                }
            }
            SymbolicExpr::Constant(_) => SymbolicExpr::Zero,
            SymbolicExpr::Add(left, right) => SymbolicExpr::Add(
                Box::new(left.derivative(var)),
                Box::new(right.derivative(var)),
            ),
            SymbolicExpr::Subtract(left, right) => SymbolicExpr::Subtract(
                Box::new(left.derivative(var)),
                Box::new(right.derivative(var)),
            ),
            SymbolicExpr::Multiply(left, right) => {
                // Product rule: (fg)' = f'g + fg'
                SymbolicExpr::Add(
                    Box::new(SymbolicExpr::Multiply(
                        Box::new(left.derivative(var)),
                        Box::new((**right).clone()),
                    )),
                    Box::new(SymbolicExpr::Multiply(
                        Box::new((**left).clone()),
                        Box::new(right.derivative(var)),
                    )),
                )
            }
            SymbolicExpr::Divide(left, right) => {
                // Quotient rule: (f/g)' = (f'g - fg')/g^2
                SymbolicExpr::Divide(
                    Box::new(SymbolicExpr::Subtract(
                        Box::new(SymbolicExpr::Multiply(
                            Box::new(left.derivative(var)),
                            Box::new((**right).clone()),
                        )),
                        Box::new(SymbolicExpr::Multiply(
                            Box::new((**left).clone()),
                            Box::new(right.derivative(var)),
                        )),
                    )),
                    Box::new(SymbolicExpr::Power(
                        Box::new((**right).clone()),
                        Box::new(SymbolicExpr::Constant(2.0)),
                    )),
                )
            }
            SymbolicExpr::Power(base, exp) => {
                match (&**exp, &**base) {
                    // Constant exponent: (f^n)' = n*f^(n-1)*f'
                    (SymbolicExpr::Constant(n), _) => SymbolicExpr::Multiply(
                        Box::new(SymbolicExpr::Multiply(
                            Box::new(SymbolicExpr::Constant(*n)),
                            Box::new(SymbolicExpr::Power(
                                Box::new((**base).clone()),
                                Box::new(SymbolicExpr::Constant(n - 1.0)),
                            )),
                        )),
                        Box::new(base.derivative(var)),
                    ),
                    // For more complex cases, we'd need logarithmic differentiation
                    _ => SymbolicExpr::Zero, // Simplified for now
                }
            }
            SymbolicExpr::Function(name, args) => {
                // Basic derivatives for common functions
                match name.as_str() {
                    "sin" if args.len() == 1 => SymbolicExpr::Multiply(
                        Box::new(SymbolicExpr::Function("cos".to_string(), args.clone())),
                        Box::new(args[0].derivative(var)),
                    ),
                    "cos" if args.len() == 1 => SymbolicExpr::Multiply(
                        Box::new(SymbolicExpr::Subtract(
                            Box::new(SymbolicExpr::Zero),
                            Box::new(SymbolicExpr::Function("sin".to_string(), args.clone())),
                        )),
                        Box::new(args[0].derivative(var)),
                    ),
                    _ => SymbolicExpr::Zero, // Unknown function, assume constant for now
                }
            }
            SymbolicExpr::Zero => SymbolicExpr::Zero,
            SymbolicExpr::One => SymbolicExpr::Zero,
        }
    }

    pub fn is_zero(&self) -> bool {
        matches!(self, SymbolicExpr::Zero)
            || matches!(self, SymbolicExpr::Constant(val) if *val == 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple() {
        assert_eq!(
            SymbolicExpr::parse("x").unwrap(),
            SymbolicExpr::Variable("x".to_string())
        );
        assert_eq!(
            SymbolicExpr::parse("42").unwrap(),
            SymbolicExpr::Constant(42.0)
        );
        assert_eq!(SymbolicExpr::parse("0").unwrap(), SymbolicExpr::Zero);
        assert_eq!(SymbolicExpr::parse("1").unwrap(), SymbolicExpr::One);
    }

    #[test]
    fn test_parse_power() {
        match SymbolicExpr::parse("r^2").unwrap() {
            SymbolicExpr::Power(base, exp) => {
                assert_eq!(*base, SymbolicExpr::Variable("r".to_string()));
                assert_eq!(*exp, SymbolicExpr::Constant(2.0));
            }
            _ => panic!("Expected power expression"),
        }
    }

    #[test]
    fn test_simplify() {
        let expr = SymbolicExpr::Add(
            Box::new(SymbolicExpr::Constant(2.0)),
            Box::new(SymbolicExpr::Constant(3.0)),
        );
        assert_eq!(expr.simplify(), SymbolicExpr::Constant(5.0));
    }

    #[test]
    fn test_derivative() {
        let expr = SymbolicExpr::Variable("x".to_string());
        assert_eq!(expr.derivative("x"), SymbolicExpr::One);
        assert_eq!(expr.derivative("y"), SymbolicExpr::Zero);
    }
}
