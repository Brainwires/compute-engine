//! Common metric tensors for physics and general relativity

use super::{Expr, SymbolicMatrix, SymbolicResult};

/// Create Minkowski metric tensor for special relativity
/// η = diag(1, -1, -1, -1) in signature (+,-,-,-)
/// or η = diag(-1, 1, 1, 1) in signature (-,+,+,+)
pub fn minkowski_metric(signature_positive_time: bool) -> SymbolicResult<SymbolicMatrix> {
    let sign = if signature_positive_time { 1 } else { -1 };

    SymbolicMatrix::new(vec![
        vec![Expr::num(sign), Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(-sign), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(-sign), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), Expr::num(-sign)],
    ])
}

/// Create 2D Minkowski metric
pub fn minkowski_2d(signature_positive_time: bool) -> SymbolicResult<SymbolicMatrix> {
    let sign = if signature_positive_time { 1 } else { -1 };

    SymbolicMatrix::new(vec![
        vec![Expr::num(sign), Expr::num(0)],
        vec![Expr::num(0), Expr::num(-sign)],
    ])
}

/// Create Euclidean metric (identity matrix)
pub fn euclidean_metric(n: usize) -> SymbolicMatrix {
    SymbolicMatrix::identity(n)
}

/// Create Schwarzschild metric in Schwarzschild coordinates
/// ds² = -(1-2M/r)dt² + (1-2M/r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
pub fn schwarzschild_metric() -> SymbolicResult<SymbolicMatrix> {
    // Symbolic metric coefficients
    let r = Expr::sym("r");
    let m = Expr::sym("M");

    // 1 - 2M/r
    let f = Expr::add(
        Expr::num(1),
        Expr::mul(
            Expr::num(-2),
            Expr::mul(m.clone(), Expr::pow(r.clone(), Expr::num(-1)))
        )
    );

    // g_tt = -(1 - 2M/r)
    let g_tt = Expr::mul(Expr::num(-1), f.clone());

    // g_rr = (1 - 2M/r)^(-1)
    let g_rr = Expr::pow(f, Expr::num(-1));

    // g_θθ = r²
    let g_theta_theta = Expr::pow(r.clone(), Expr::num(2));

    // g_φφ = r²sin²θ
    let g_phi_phi = Expr::mul(
        Expr::pow(r, Expr::num(2)),
        Expr::pow(Expr::func("sin", vec![Expr::sym("θ")]), Expr::num(2))
    );

    SymbolicMatrix::new(vec![
        vec![g_tt, Expr::num(0), Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), g_rr, Expr::num(0), Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), g_theta_theta, Expr::num(0)],
        vec![Expr::num(0), Expr::num(0), Expr::num(0), g_phi_phi],
    ])
}

/// Create Kerr metric (rotating black hole) - simplified 2D version
pub fn kerr_metric_2d() -> SymbolicResult<SymbolicMatrix> {
    let r = Expr::sym("r");
    let m = Expr::sym("M");
    let a = Expr::sym("a"); // angular momentum parameter

    // Simplified: g_tt = -(1 - 2Mr/(r² + a²))
    let r_squared = Expr::pow(r.clone(), Expr::num(2));
    let a_squared = Expr::pow(a, Expr::num(2));

    let denom = Expr::add(r_squared.clone(), a_squared);
    let factor = Expr::mul(
        Expr::mul(Expr::num(2), m),
        Expr::mul(r.clone(), Expr::pow(denom.clone(), Expr::num(-1)))
    );

    let g_tt = Expr::mul(Expr::num(-1), Expr::add(Expr::num(1), Expr::mul(Expr::num(-1), factor)));

    // g_rr = simplified
    let g_rr = Expr::pow(denom.clone(), Expr::num(-1));

    SymbolicMatrix::new(vec![
        vec![g_tt, Expr::num(0)],
        vec![Expr::num(0), g_rr],
    ])
}

/// Create Friedmann-Lemaître-Robertson-Walker (FLRW) metric for cosmology
/// ds² = -dt² + a²(t)[dr²/(1-kr²) + r²dΩ²]
/// Simplified 2D version
pub fn flrw_metric_2d() -> SymbolicResult<SymbolicMatrix> {
    let a = Expr::sym("a"); // scale factor
    let k = Expr::sym("k"); // curvature parameter
    let r = Expr::sym("r");

    // g_tt = -1
    let g_tt = Expr::num(-1);

    // g_rr = a²/(1-kr²)
    let a_squared = Expr::pow(a, Expr::num(2));
    let one_minus_kr2 = Expr::add(
        Expr::num(1),
        Expr::mul(
            Expr::num(-1),
            Expr::mul(k, Expr::pow(r, Expr::num(2)))
        )
    );
    let g_rr = Expr::mul(a_squared, Expr::pow(one_minus_kr2, Expr::num(-1)));

    SymbolicMatrix::new(vec![
        vec![g_tt, Expr::num(0)],
        vec![Expr::num(0), g_rr],
    ])
}

/// Create a general diagonal metric from symbolic expressions
pub fn diagonal_metric(diagonal_elements: Vec<Expr>) -> SymbolicResult<SymbolicMatrix> {
    let n = diagonal_elements.len();
    let mut data = vec![vec![Expr::num(0); n]; n];

    for (i, elem) in diagonal_elements.into_iter().enumerate() {
        data[i][i] = elem;
    }

    SymbolicMatrix::new(data)
}

/// Create an anti-de Sitter (AdS) metric in Poincaré coordinates
/// ds² = (L²/z²)(-dt² + dx² + dy² + dz²)
/// Simplified 2D version
pub fn ads_metric_2d() -> SymbolicResult<SymbolicMatrix> {
    let l = Expr::sym("L"); // AdS radius
    let z = Expr::sym("z");

    // Conformal factor: L²/z²
    let conformal_factor = Expr::mul(
        Expr::pow(l, Expr::num(2)),
        Expr::pow(z, Expr::num(-2))
    );

    // g_tt = -L²/z²
    let g_tt = Expr::mul(Expr::num(-1), conformal_factor.clone());

    // g_zz = L²/z²
    let g_zz = conformal_factor;

    SymbolicMatrix::new(vec![
        vec![g_tt, Expr::num(0)],
        vec![Expr::num(0), g_zz],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minkowski_metric() {
        let metric = minkowski_metric(true).unwrap();
        assert_eq!(metric.rows(), 4);
        assert_eq!(metric.cols(), 4);

        // Check signature
        assert_eq!(metric.get(0, 0), Some(&Expr::num(1)));
        assert_eq!(metric.get(1, 1), Some(&Expr::num(-1)));
    }

    #[test]
    fn test_euclidean_metric() {
        let metric = euclidean_metric(3);
        assert_eq!(metric.rows(), 3);

        // Should be identity
        assert_eq!(metric.get(0, 0), Some(&Expr::num(1)));
        assert_eq!(metric.get(0, 1), Some(&Expr::num(0)));
    }

    #[test]
    fn test_schwarzschild_metric() {
        let metric = schwarzschild_metric().unwrap();
        assert_eq!(metric.rows(), 4);
        println!("Schwarzschild metric:");
        println!("{}", metric);
    }

    #[test]
    fn test_kerr_metric() {
        let metric = kerr_metric_2d().unwrap();
        println!("Kerr metric (2D):");
        println!("{}", metric);
    }

    #[test]
    fn test_flrw_metric() {
        let metric = flrw_metric_2d().unwrap();
        println!("FLRW metric (2D):");
        println!("{}", metric);
    }

    #[test]
    fn test_ads_metric() {
        let metric = ads_metric_2d().unwrap();
        println!("AdS metric (2D):");
        println!("{}", metric);
    }
}
