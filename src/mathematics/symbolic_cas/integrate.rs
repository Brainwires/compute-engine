//! Symbolic integration
//!
//! Computes indefinite integrals of expressions

use super::expr::Expr;
use super::simplify::simplify;

/// Compute the indefinite integral of an expression with respect to a variable
pub fn integrate(expr: &Expr, var: &str) -> Expr {
    let integral = integrate_internal(expr, var);
    simplify(&integral)
}

fn integrate_internal(expr: &Expr, var: &str) -> Expr {
    match expr {
        // ∫ c dx = c*x
        Expr::Number(n) => Expr::mul(Expr::Number(n.clone()), Expr::sym(var)),

        // ∫ x dx = x^2/2, ∫ y dx = y*x
        Expr::Symbol(s) => {
            if s == var {
                // ∫ x dx = x^2/2
                Expr::mul(
                    Expr::rational_unchecked(1, 2),
                    Expr::pow(Expr::sym(var), Expr::num(2))
                )
            } else {
                // ∫ y dx = y*x (constant with respect to x)
                Expr::mul(Expr::sym(s), Expr::sym(var))
            }
        }

        // ∫ (a + b) dx = ∫ a dx + ∫ b dx
        Expr::Add(a, b) => {
            Expr::add(integrate_internal(a, var), integrate_internal(b, var))
        }

        // ∫ a * b dx - try to factor out constants
        Expr::Mul(a, b) => {
            // Check if 'a' is constant with respect to var
            if is_constant(a, var) {
                // ∫ c * f(x) dx = c * ∫ f(x) dx
                Expr::mul((**a).clone(), integrate_internal(b, var))
            } else if is_constant(b, var) {
                // ∫ f(x) * c dx = c * ∫ f(x) dx
                Expr::mul((**b).clone(), integrate_internal(a, var))
            } else {
                // Can't integrate product in general - would need integration by parts
                Expr::func(format!("∫({} * {}) d{}", a, b, var), vec![])
            }
        }

        // Power rule: ∫ x^n dx = x^(n+1)/(n+1) for n ≠ -1
        Expr::Pow(base, exp) => {
            // Check if base is the variable and exp is constant
            if let Expr::Symbol(s) = &**base {
                if s == var && is_constant(exp, var) {
                    // ∫ x^n dx = x^(n+1)/(n+1)
                    let new_exp = Expr::add((**exp).clone(), Expr::num(1));
                    return Expr::mul(
                        Expr::pow(new_exp.clone(), Expr::num(-1)),  // 1/(n+1)
                        Expr::pow(Expr::sym(var), new_exp)  // x^(n+1)
                    );
                }
            }

            // Can't integrate in general
            Expr::func(format!("∫({}^{}) d{}", base, exp, var), vec![])
        }

        // Function integrals
        Expr::Function(name, args) => {
            if args.len() != 1 {
                return Expr::func(format!("∫({}) d{}", name, var), args.clone());
            }

            let arg = &args[0];

            // Only handle simple cases where arg is exactly the variable
            if let Expr::Symbol(s) = arg {
                if s == var {
                    match name.as_str() {
                        // ∫ sin(x) dx = -cos(x)
                        "sin" => return Expr::mul(
                            Expr::num(-1),
                            Expr::func("cos", vec![Expr::sym(var)])
                        ),

                        // ∫ cos(x) dx = sin(x)
                        "cos" => return Expr::func("sin", vec![Expr::sym(var)]),

                        // ∫ exp(x) dx = exp(x)
                        "exp" => return Expr::func("exp", vec![Expr::sym(var)]),

                        _ => {}
                    }
                }
            }

            // Can't integrate general functions
            Expr::func(format!("∫{} d{}", name, var), args.clone())
        }

        _ => Expr::func(format!("∫({}) d{}", expr, var), vec![])
    }
}

/// Check if an expression is constant with respect to a variable
fn is_constant(expr: &Expr, var: &str) -> bool {
    match expr {
        Expr::Number(_) => true,
        Expr::Symbol(s) => s != var,
        Expr::Add(a, b) => is_constant(a, var) && is_constant(b, var),
        Expr::Mul(a, b) => is_constant(a, var) && is_constant(b, var),
        Expr::Pow(base, exp) => is_constant(base, var) && is_constant(exp, var),
        Expr::Function(_, args) => args.iter().all(|arg| is_constant(arg, var)),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mathematics::symbolic_cas::parser::parse;

    #[test]
    fn test_integrate_constant() {
        let expr = Expr::num(5);
        let integral = integrate(&expr, "x");
        println!("∫ 5 dx = {}", integral);
    }

    #[test]
    fn test_integrate_variable() {
        let expr = Expr::sym("x");
        let integral = integrate(&expr, "x");
        println!("∫ x dx = {}", integral);
    }

    #[test]
    fn test_integrate_power() {
        let expr = parse("x^2").unwrap();
        let integral = integrate(&expr, "x");
        println!("∫ x^2 dx = {}", integral);
    }

    #[test]
    fn test_integrate_polynomial() {
        let expr = parse("2*x + 3").unwrap();
        let integral = integrate(&expr, "x");
        println!("∫ (2*x + 3) dx = {}", integral);
    }
}
