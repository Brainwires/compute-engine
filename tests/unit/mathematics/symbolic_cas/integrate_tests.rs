// Unit tests for mathematics::symbolic_cas::integrate
use computational_engine::mathematics::symbolic_cas::integrate::*;

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
