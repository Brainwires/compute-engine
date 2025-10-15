use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

#[derive(Debug, Serialize, Deserialize)]
pub struct SpecialFunctionsResult {
    pub success: bool,
    pub operation: String,
    pub result: Option<Value>,
    pub error: Option<String>,
}

// Riemann zeta function using Euler-Maclaurin formula
pub fn riemann_zeta(s: f64) -> f64 {
    if s == 1.0 {
        return f64::INFINITY;
    }
    
    let n_terms = 1000;
    let mut sum = 0.0;
    
    // Direct summation for the first part
    for k in 1..=n_terms {
        sum += 1.0 / (k as f64).powf(s);
    }
    
    // Euler-Maclaurin correction terms for better convergence
    let n = n_terms as f64;
    let correction = 1.0 / (2.0 * n.powf(s)) + 
                     s / (12.0 * n.powf(s + 1.0)) -
                     s * (s + 1.0) / (120.0 * n.powf(s + 3.0));
    
    sum + correction
}

// Elliptic integrals of first and second kind
pub fn elliptic_integral_first_kind(k: f64) -> f64 {
    if k.abs() >= 1.0 {
        return f64::NAN;
    }
    
    let n_steps = 10000;
    let dt = std::f64::consts::PI / (2.0 * n_steps as f64);
    let mut sum = 0.0;
    
    for i in 0..n_steps {
        let t = (i as f64 + 0.5) * dt;
        let sin_t = t.sin();
        let integrand = 1.0 / (1.0 - k * k * sin_t * sin_t).sqrt();
        sum += integrand * dt;
    }
    
    sum
}

pub fn elliptic_integral_second_kind(k: f64) -> f64 {
    if k.abs() >= 1.0 {
        return f64::NAN;
    }
    
    let n_steps = 10000;
    let dt = std::f64::consts::PI / (2.0 * n_steps as f64);
    let mut sum = 0.0;
    
    for i in 0..n_steps {
        let t = (i as f64 + 0.5) * dt;
        let sin_t = t.sin();
        let integrand = (1.0 - k * k * sin_t * sin_t).sqrt();
        sum += integrand * dt;
    }
    
    sum
}

// Hypergeometric function 2F1(a, b; c; z)
pub fn hypergeometric_2f1(a: f64, b: f64, c: f64, z: f64) -> f64 {
    if z.abs() >= 1.0 {
        return f64::NAN; // Series doesn't converge for |z| >= 1
    }
    
    let max_terms = 1000;
    let tolerance = 1e-12;
    let mut sum = 1.0;
    let mut term = 1.0;
    
    for n in 1..=max_terms {
        let n_f64 = n as f64;
        term *= (a + n_f64 - 1.0) * (b + n_f64 - 1.0) * z / ((c + n_f64 - 1.0) * n_f64);
        sum += term;
        
        if term.abs() < tolerance {
            break;
        }
    }
    
    sum
}

// Jacobi theta functions
pub fn jacobi_theta_3(q: f64, z: f64) -> f64 {
    let max_terms = 200;
    let mut sum = 1.0;
    
    for n in 1..=max_terms {
        let n_f64 = n as f64;
        let q_power = q.powf(n_f64 * n_f64);
        sum += 2.0 * q_power * (2.0 * n_f64 * z).cos();
        
        if q_power < 1e-15 {
            break;
        }
    }
    
    sum
}

pub fn jacobi_theta_4(q: f64, z: f64) -> f64 {
    let max_terms = 200;
    let mut sum = 1.0;
    
    for n in 1..=max_terms {
        let n_f64 = n as f64;
        let q_power = q.powf(n_f64 * n_f64);
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        sum += 2.0 * sign * q_power * (2.0 * n_f64 * z).cos();
        
        if q_power < 1e-15 {
            break;
        }
    }
    
    sum
}

// Bessel functions of the first kind (using series expansion)
pub fn bessel_j0(x: f64) -> f64 {
    if x.abs() < 1e-8 {
        return 1.0;
    }
    
    let max_terms = 100;
    let tolerance = 1e-12;
    let mut sum = 1.0;
    let mut term = 1.0;
    let x_half_squared = (x / 2.0).powi(2);
    
    for k in 1..=max_terms {
        let k_f64 = k as f64;
        term *= -x_half_squared / (k_f64 * k_f64);
        sum += term;
        
        if term.abs() < tolerance {
            break;
        }
    }
    
    sum
}

pub fn bessel_j1(x: f64) -> f64 {
    if x.abs() < 1e-8 {
        return 0.0;
    }
    
    let max_terms = 100;
    let tolerance = 1e-12;
    let x_half = x / 2.0;
    let mut sum = x_half;
    let mut term = x_half;
    let x_half_squared = x_half * x_half;
    
    for k in 1..=max_terms {
        let k_f64 = k as f64;
        term *= -x_half_squared / (k_f64 * (k_f64 + 1.0));
        sum += term;
        
        if term.abs() < tolerance {
            break;
        }
    }
    
    sum
}

// Legendre polynomials
pub fn legendre_polynomial(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    } else if n == 1 {
        return x;
    }
    
    // Use recurrence relation: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
    let mut p_prev = 1.0;  // P_0(x)
    let mut p_curr = x;    // P_1(x)
    
    for k in 1..n {
        let k_f64 = k as f64;
        let p_next = ((2.0 * k_f64 + 1.0) * x * p_curr - k_f64 * p_prev) / (k_f64 + 1.0);
        p_prev = p_curr;
        p_curr = p_next;
    }
    
    p_curr
}

// Handler functions
pub fn handle_riemann_zeta(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let s = params.get("s")
        .and_then(|v| v.as_f64())
        .unwrap_or(2.0);
    
    let result = riemann_zeta(s);
    
    SpecialFunctionsResult {
        success: true,
        operation: "riemann_zeta".to_string(),
        result: Some(serde_json::json!({
            "s": s,
            "zeta_value": result
        })),
        error: None,
    }
}

pub fn handle_elliptic_integral(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let k = params.get("k")
        .and_then(|v| v.as_f64())
        .unwrap_or(0.5);
    let kind = params.get("type")
        .and_then(|v| v.as_str())
        .unwrap_or("first");
    
    let result = match kind {
        "first" => elliptic_integral_first_kind(k),
        "second" => elliptic_integral_second_kind(k),
        _ => elliptic_integral_first_kind(k),
    };
    
    SpecialFunctionsResult {
        success: true,
        operation: "elliptic_integral".to_string(),
        result: Some(serde_json::json!({
            "k": k,
            "type": kind,
            "value": result
        })),
        error: None,
    }
}

pub fn handle_hypergeometric(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let a = params.get("a").and_then(|v| v.as_f64()).unwrap_or(1.0);
    let b = params.get("b").and_then(|v| v.as_f64()).unwrap_or(1.0);
    let c = params.get("c").and_then(|v| v.as_f64()).unwrap_or(2.0);
    let z = params.get("z").and_then(|v| v.as_f64()).unwrap_or(0.5);
    
    let result = hypergeometric_2f1(a, b, c, z);
    
    SpecialFunctionsResult {
        success: true,
        operation: "hypergeometric".to_string(),
        result: Some(serde_json::json!({
            "a": a, "b": b, "c": c, "z": z,
            "value": result
        })),
        error: None,
    }
}

pub fn handle_jacobi_theta(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let q = params.get("q").and_then(|v| v.as_f64()).unwrap_or(0.1);
    let z = params.get("z").and_then(|v| v.as_f64()).unwrap_or(0.0);
    let variant = params.get("variant").and_then(|v| v.as_u64()).unwrap_or(3);
    
    let result = match variant {
        3 => jacobi_theta_3(q, z),
        4 => jacobi_theta_4(q, z),
        _ => jacobi_theta_3(q, z),
    };
    
    SpecialFunctionsResult {
        success: true,
        operation: "jacobi_theta".to_string(),
        result: Some(serde_json::json!({
            "q": q, "z": z, "variant": variant,
            "value": result
        })),
        error: None,
    }
}

pub fn handle_bessel_function(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let x = params.get("x").and_then(|v| v.as_f64()).unwrap_or(1.0);
    let order = params.get("order").and_then(|v| v.as_u64()).unwrap_or(0);
    
    let result = match order {
        0 => bessel_j0(x),
        1 => bessel_j1(x),
        _ => bessel_j0(x), // Default to J_0
    };
    
    SpecialFunctionsResult {
        success: true,
        operation: "bessel_function".to_string(),
        result: Some(serde_json::json!({
            "x": x, "order": order,
            "value": result
        })),
        error: None,
    }
}

pub fn handle_legendre_polynomial(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let n = params.get("n").and_then(|v| v.as_u64()).unwrap_or(2) as usize;
    let x = params.get("x").and_then(|v| v.as_f64()).unwrap_or(0.5);
    
    let result = legendre_polynomial(n, x);
    
    SpecialFunctionsResult {
        success: true,
        operation: "legendre_polynomial".to_string(),
        result: Some(serde_json::json!({
            "n": n, "x": x,
            "value": result
        })),
        error: None,
    }
}

pub fn handle_special_functions(params: &HashMap<String, Value>) -> SpecialFunctionsResult {
    let function_type = params.get("function_type")
        .and_then(|v| v.as_str())
        .unwrap_or("riemann_zeta");
    
    match function_type {
        "riemann_zeta" => handle_riemann_zeta(params),
        "elliptic_integral" => handle_elliptic_integral(params),
        "hypergeometric" => handle_hypergeometric(params),
        "jacobi_theta" => handle_jacobi_theta(params),
        "bessel" => handle_bessel_function(params),
        "legendre" => handle_legendre_polynomial(params),
        _ => SpecialFunctionsResult {
            success: false,
            operation: "special_functions".to_string(),
            result: None,
            error: Some(format!("Unknown special function type: {}", function_type)),
        }
    }
}