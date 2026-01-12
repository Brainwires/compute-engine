//! Symbolic differentiation
//!
//! Computes derivatives of expressions

use super::expr::Expr;
use super::simplify::simplify;

/// Compute the derivative of an expression with respect to a variable
pub fn differentiate(expr: &Expr, var: &str) -> Expr {
    let derivative = diff_internal(expr, var);
    simplify(&derivative)
}

fn diff_internal(expr: &Expr, var: &str) -> Expr {
    match expr {
        // d/dx(c) = 0
        Expr::Number(_) => Expr::num(0),

        // d/dx(x) = 1, d/dx(y) = 0
        Expr::Symbol(s) => {
            if s == var {
                Expr::num(1)
            } else {
                Expr::num(0)
            }
        }

        // d/dx(a + b) = da/dx + db/dx
        Expr::Add(a, b) => Expr::add(diff_internal(a, var), diff_internal(b, var)),

        // Product rule: d/dx(a * b) = a * db/dx + da/dx * b
        Expr::Mul(a, b) => {
            let da = diff_internal(a, var);
            let db = diff_internal(b, var);
            Expr::add(Expr::mul((**a).clone(), db), Expr::mul(da, (**b).clone()))
        }

        // Power rule: d/dx(x^n) = n * x^(n-1)
        // Chain rule: d/dx(f^n) = n * f^(n-1) * df/dx
        Expr::Pow(base, exp) => {
            let df = diff_internal(base, var);

            // If exponent is constant: n * base^(n-1) * dbase/dx
            Expr::mul(
                (**exp).clone(),
                Expr::mul(
                    Expr::pow((**base).clone(), Expr::add((**exp).clone(), Expr::num(-1))),
                    df,
                ),
            )
        }

        // Function derivatives
        Expr::Function(name, args) => {
            if args.len() != 1 {
                // Multi-arg functions not yet supported
                return Expr::func(format!("d/d{}({})", var, name), args.clone());
            }

            let arg = &args[0];
            let darg = diff_internal(arg, var);

            match name.as_str() {
                // d/dx(sin(f)) = cos(f) * df/dx
                "sin" => Expr::mul(Expr::func("cos", vec![arg.clone()]), darg),

                // d/dx(cos(f)) = -sin(f) * df/dx
                "cos" => Expr::mul(
                    Expr::num(-1),
                    Expr::mul(Expr::func("sin", vec![arg.clone()]), darg),
                ),

                // d/dx(exp(f)) = exp(f) * df/dx
                "exp" => Expr::mul(Expr::func("exp", vec![arg.clone()]), darg),

                // d/dx(ln(f)) = df/dx / f
                "ln" => Expr::mul(darg, Expr::pow(arg.clone(), Expr::num(-1))),

                // d/dx(sqrt(f)) = df/dx / (2*sqrt(f))
                "sqrt" => Expr::mul(
                    darg,
                    Expr::pow(
                        Expr::mul(Expr::num(2), Expr::func("sqrt", vec![arg.clone()])),
                        Expr::num(-1),
                    ),
                ),

                _ => Expr::func(format!("d/d{}({})", var, name), vec![arg.clone()]),
            }
        }
    }
}

