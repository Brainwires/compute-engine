//! Core symbolic expression representation
//!
//! Provides a custom expression tree for symbolic mathematics

use super::{SymbolicError, SymbolicResult};
use std::collections::HashMap;
use std::fmt;

/// Symbolic expression tree
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// Constant number (integer or rational)
    Number(Rational),

    /// Variable (x, y, z, etc.)
    Symbol(String),

    /// Addition: a + b
    Add(Box<Expr>, Box<Expr>),

    /// Multiplication: a * b
    Mul(Box<Expr>, Box<Expr>),

    /// Power: a^b
    Pow(Box<Expr>, Box<Expr>),

    /// Function call: f(args)
    Function(String, Vec<Expr>),
}

/// Rational number representation
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rational {
    pub numerator: i64,
    pub denominator: i64,
}

impl Rational {
    pub fn new(n: i64, d: i64) -> Result<Self, String> {
        if d == 0 {
            return Err("Denominator cannot be zero".to_string());
        }
        let gcd = gcd(n.abs(), d.abs());
        let sign = if (n < 0) ^ (d < 0) { -1 } else { 1 };
        Ok(Rational {
            numerator: sign * n.abs() / gcd,
            denominator: d.abs() / gcd,
        })
    }

    pub fn from_int(n: i64) -> Self {
        Rational { numerator: n, denominator: 1 }
    }

    pub fn to_f64(&self) -> f64 {
        self.numerator as f64 / self.denominator as f64
    }

    pub fn is_zero(&self) -> bool {
        self.numerator == 0
    }

    pub fn is_one(&self) -> bool {
        self.numerator == self.denominator
    }

    pub fn is_integer(&self) -> bool {
        self.denominator == 1
    }
}

impl std::ops::Add for Rational {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        // Denominators are always non-zero for valid Rationals
        Rational::new(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator
        ).expect("BUG: denominator became zero in Add")
    }
}

impl std::ops::Mul for Rational {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        // Denominators are always non-zero for valid Rationals
        Rational::new(
            self.numerator * other.numerator,
            self.denominator * other.denominator
        ).expect("BUG: denominator became zero in Mul")
    }
}

impl std::ops::Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self {
        // Denominator is always non-zero for valid Rationals
        Rational::new(-self.numerator, self.denominator)
            .expect("BUG: denominator became zero in Neg")
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

impl Expr {
    /// Create a number from an integer
    pub fn num(n: i64) -> Self {
        Expr::Number(Rational::from_int(n))
    }

    /// Create a rational number
    pub fn rational(n: i64, d: i64) -> SymbolicResult<Self> {
        match Rational::new(n, d) {
            Ok(r) => Ok(Expr::Number(r)),
            Err(e) => Err(SymbolicError::ParseError(e)),
        }
    }

    /// Create a rational number (panics if denominator is zero)
    ///
    /// This is a convenience method for internal use where denominators are known to be non-zero.
    /// For user input or external data, use `rational()` which returns Result.
    pub fn rational_unchecked(n: i64, d: i64) -> Self {
        Expr::Number(Rational::new(n, d).expect("BUG: rational_unchecked called with zero denominator"))
    }

    /// Create a symbol/variable
    pub fn sym(s: impl Into<String>) -> Self {
        Expr::Symbol(s.into())
    }

    /// Create an addition
    pub fn add(a: Expr, b: Expr) -> Self {
        Expr::Add(Box::new(a), Box::new(b))
    }

    /// Create a multiplication
    pub fn mul(a: Expr, b: Expr) -> Self {
        Expr::Mul(Box::new(a), Box::new(b))
    }

    /// Create a power
    pub fn pow(base: Expr, exp: Expr) -> Self {
        Expr::Pow(Box::new(base), Box::new(exp))
    }

    /// Create a function call
    pub fn func(name: impl Into<String>, args: Vec<Expr>) -> Self {
        Expr::Function(name.into(), args)
    }

    /// Evaluate the expression with given variable values
    pub fn evaluate(&self, vars: &HashMap<String, f64>) -> SymbolicResult<f64> {
        match self {
            Expr::Number(r) => Ok(r.to_f64()),
            Expr::Symbol(s) => vars.get(s)
                .copied()
                .ok_or_else(|| SymbolicError::EvaluationError(format!("Unknown variable: {}", s))),
            Expr::Add(a, b) => Ok(a.evaluate(vars)? + b.evaluate(vars)?),
            Expr::Mul(a, b) => Ok(a.evaluate(vars)? * b.evaluate(vars)?),
            Expr::Pow(base, exp) => Ok(base.evaluate(vars)?.powf(exp.evaluate(vars)?)),
            Expr::Function(name, args) => {
                match name.as_str() {
                    "sin" if args.len() == 1 => Ok(args[0].evaluate(vars)?.sin()),
                    "cos" if args.len() == 1 => Ok(args[0].evaluate(vars)?.cos()),
                    "exp" if args.len() == 1 => Ok(args[0].evaluate(vars)?.exp()),
                    "ln" if args.len() == 1 => Ok(args[0].evaluate(vars)?.ln()),
                    "sqrt" if args.len() == 1 => Ok(args[0].evaluate(vars)?.sqrt()),
                    _ => Err(SymbolicError::EvaluationError(format!("Unknown function: {}", name)))
                }
            }
        }
    }

    /// Check if expression is zero
    pub fn is_zero(&self) -> bool {
        matches!(self, Expr::Number(r) if r.is_zero())
    }

    /// Check if expression is one
    pub fn is_one(&self) -> bool {
        matches!(self, Expr::Number(r) if r.is_one())
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Expr::Number(r) => write!(f, "{}", r),
            Expr::Symbol(s) => write!(f, "{}", s),
            Expr::Add(a, b) => write!(f, "({} + {})", a, b),
            Expr::Mul(a, b) => write!(f, "({} * {})", a, b),
            Expr::Pow(base, exp) => write!(f, "{}^{}", base, exp),
            Expr::Function(name, args) => {
                write!(f, "{}(", name)?;
                for (i, arg) in args.iter().enumerate() {
                    if i > 0 { write!(f, ", ")?; }
                    write!(f, "{}", arg)?;
                }
                write!(f, ")")
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rational() {
        let r = Rational::new(6, 9).unwrap();
        assert_eq!(r.numerator, 2);
        assert_eq!(r.denominator, 3);
    }

    #[test]
    fn test_rational_zero_denominator() {
        let result = Rational::new(5, 0);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Denominator cannot be zero");
    }

    #[test]
    fn test_expr_display() {
        let expr = Expr::add(Expr::sym("x"), Expr::num(1));
        assert_eq!(format!("{}", expr), "(x + 1)");
    }

    #[test]
    fn test_evaluate() {
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), 2.0);

        let expr = Expr::add(Expr::sym("x"), Expr::num(3));
        assert_eq!(expr.evaluate(&vars).unwrap(), 5.0);
    }
}
