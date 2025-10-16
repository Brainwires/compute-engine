//! Expression simplification
//!
//! Applies algebraic simplification rules

use super::expr::Expr;

/// Simplify an expression
pub fn simplify(expr: &Expr) -> Expr {
    let mut current = expr.clone();
    let mut changed = true;

    // Apply simplification rules repeatedly until no changes
    for _ in 0..10 {
        if !changed {
            break;
        }
        let new_expr = simplify_once(&current);
        changed = new_expr != current;
        current = new_expr;
    }

    current
}

fn simplify_once(expr: &Expr) -> Expr {
    match expr {
        // Simplify nested expressions first
        Expr::Add(a, b) => simplify_add(&simplify_once(a), &simplify_once(b)),
        Expr::Mul(a, b) => simplify_mul(&simplify_once(a), &simplify_once(b)),
        Expr::Pow(base, exp) => simplify_pow(&simplify_once(base), &simplify_once(exp)),
        Expr::Function(name, args) => {
            let simplified_args: Vec<Expr> = args.iter().map(simplify_once).collect();
            Expr::func(name.clone(), simplified_args)
        }
        _ => expr.clone(),
    }
}

fn simplify_add(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        // 0 + x = x
        (Expr::Number(r), x) | (x, Expr::Number(r)) if r.is_zero() => x.clone(),

        // n + m = (n+m)
        (Expr::Number(r1), Expr::Number(r2)) => Expr::Number(r1.clone() + r2.clone()),

        // x + x = 2*x
        (x, y) if x == y => Expr::mul(Expr::num(2), x.clone()),

        // Collect like terms: a*x + b*x = (a+b)*x
        (Expr::Mul(a1, x1), Expr::Mul(a2, x2)) if x1 == x2 => {
            let coef = simplify_add(a1, a2);
            simplify_mul(&coef, x1)
        }

        _ => Expr::add(a.clone(), b.clone()),
    }
}

fn simplify_mul(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        // 0 * x = 0
        (Expr::Number(r), _) | (_, Expr::Number(r)) if r.is_zero() => Expr::num(0),

        // 1 * x = x
        (Expr::Number(r), x) | (x, Expr::Number(r)) if r.is_one() => x.clone(),

        // n * m = (n*m)
        (Expr::Number(r1), Expr::Number(r2)) => Expr::Number(r1.clone() * r2.clone()),

        // x * x = x^2
        (x, y) if x == y => Expr::pow(x.clone(), Expr::num(2)),

        // x^a * x^b = x^(a+b)
        (Expr::Pow(base1, exp1), Expr::Pow(base2, exp2)) if base1 == base2 => {
            let new_exp = simplify_add(exp1, exp2);
            Expr::pow((**base1).clone(), new_exp)
        }

        // x * x^a = x^(a+1)
        (x, Expr::Pow(base, exp)) | (Expr::Pow(base, exp), x) if x == &**base => {
            let new_exp = simplify_add(exp, &Expr::num(1));
            Expr::pow((**base).clone(), new_exp)
        }

        _ => Expr::mul(a.clone(), b.clone()),
    }
}

fn simplify_pow(base: &Expr, exp: &Expr) -> Expr {
    match (base, exp) {
        // x^0 = 1
        (_, Expr::Number(r)) if r.is_zero() => Expr::num(1),

        // x^1 = x
        (x, Expr::Number(r)) if r.is_one() => x.clone(),

        // 0^x = 0 (for positive x)
        (Expr::Number(r), _) if r.is_zero() => Expr::num(0),

        // 1^x = 1
        (Expr::Number(r), _) if r.is_one() => Expr::num(1),

        // (x^a)^b = x^(a*b)
        (Expr::Pow(base2, exp2), exp1) => {
            let new_exp = simplify_mul(exp2, exp1);
            Expr::pow((**base2).clone(), new_exp)
        }

        _ => Expr::pow(base.clone(), exp.clone()),
    }
}

/// Expand an expression (distribute multiplication over addition)
pub fn expand(expr: &Expr) -> Expr {
    match expr {
        Expr::Mul(a, b) => expand_mul(&expand(a), &expand(b)),
        Expr::Pow(base, exp) => {
            if let Expr::Number(r) = &**exp {
                if r.is_integer() && r.numerator > 0 && r.numerator < 10 {
                    return expand_pow(&expand(base), r.numerator as u32);
                }
            }
            Expr::pow(expand(base), expand(exp))
        }
        Expr::Add(a, b) => Expr::add(expand(a), expand(b)),
        Expr::Function(name, args) => {
            let expanded_args: Vec<Expr> = args.iter().map(expand).collect();
            Expr::func(name.clone(), expanded_args)
        }
        _ => expr.clone(),
    }
}

fn expand_mul(a: &Expr, b: &Expr) -> Expr {
    match (a, b) {
        // (a + b) * c = a*c + b*c
        (Expr::Add(a1, a2), c) | (c, Expr::Add(a1, a2)) => {
            let left = expand_mul(a1, c);
            let right = expand_mul(a2, c);
            Expr::add(left, right)
        }
        _ => Expr::mul(a.clone(), b.clone()),
    }
}

fn expand_pow(base: &Expr, exp: u32) -> Expr {
    if exp == 0 {
        return Expr::num(1);
    }
    if exp == 1 {
        return base.clone();
    }

    let mut result = base.clone();
    for _ in 1..exp {
        result = expand_mul(&result, base);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mathematics::symbolic_cas::parser::parse;

    #[test]
    fn test_simplify_add_zero() {
        let expr = parse("x + 0").unwrap();
        let simplified = simplify(&expr);
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    #[test]
    fn test_simplify_mul_one() {
        let expr = parse("x * 1").unwrap();
        let simplified = simplify(&expr);
        assert!(matches!(simplified, Expr::Symbol(_)));
    }

    #[test]
    fn test_expand_square() {
        let expr = parse("(x + 1)^2").unwrap();
        let expanded = expand(&expr);
        println!("Expanded: {}", expanded);
    }
}
