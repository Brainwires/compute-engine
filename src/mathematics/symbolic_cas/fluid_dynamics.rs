//! Fluid dynamics operations
//!
//! Navier-Stokes equations, continuity equation, vorticity, and fluid flow analysis

use super::{Expr, SymbolicResult, SymbolicError};
use super::differentiate::differentiate as diff;
use super::simplify::simplify;

/// Continuity equation for incompressible flow: ∇·v = 0
///
/// ∂u/∂x + ∂v/∂y + ∂w/∂z = 0
///
/// Returns the divergence expression
pub fn continuity_equation_incompressible(
    velocity_x: &str,
    velocity_y: &str,
    velocity_z: &str,
    coords: &[&str],
) -> SymbolicResult<Expr> {
    if coords.len() != 3 {
        return Err(SymbolicError::InvalidOperation(
            "Need 3 coordinates (x, y, z)".to_string()
        ));
    }

    let u = Expr::sym(velocity_x);
    let v = Expr::sym(velocity_y);
    let w = Expr::sym(velocity_z);

    // ∂u/∂x
    let du_dx = diff(&u, coords[0]);

    // ∂v/∂y
    let dv_dy = diff(&v, coords[1]);

    // ∂w/∂z
    let dw_dz = diff(&w, coords[2]);

    // Sum
    let divergence = Expr::add(
        Expr::add(du_dx, dv_dy),
        dw_dz
    );

    Ok(simplify(&divergence))
}

/// Compute vorticity: ω = ∇ × v
///
/// For 2D flow (x-y plane):
/// ω_z = ∂v/∂x - ∂u/∂y
pub fn vorticity_2d(
    velocity_x: &str,
    velocity_y: &str,
) -> Expr {
    let u = Expr::sym(velocity_x);
    let v = Expr::sym(velocity_y);

    // ∂v/∂x
    let dv_dx = diff(&v, "x");

    // ∂u/∂y
    let du_dy = diff(&u, "y");

    // ω_z = ∂v/∂x - ∂u/∂y
    simplify(&Expr::add(dv_dx, Expr::mul(Expr::num(-1), du_dy)))
}

/// Compute stream function ψ for 2D incompressible flow
///
/// u = ∂ψ/∂y, v = -∂ψ/∂x
///
/// Returns symbolic expression for stream function
pub fn stream_function_velocity(psi: &str) -> (Expr, Expr) {
    let psi_expr = Expr::sym(psi);

    // u = ∂ψ/∂y
    let u = diff(&psi_expr, "y");

    // v = -∂ψ/∂x
    let v = Expr::mul(Expr::num(-1), diff(&psi_expr, "x"));

    (u, v)
}

/// Bernoulli equation for steady, incompressible, inviscid flow
///
/// p + (1/2)ρv² + ρgh = constant
///
/// Returns the Bernoulli constant expression
pub fn bernoulli_equation(
    pressure: &str,
    velocity: &str,
    height: &str,
    density: &str,
    gravity: &str,
) -> Expr {
    let p = Expr::sym(pressure);
    let v = Expr::sym(velocity);
    let h = Expr::sym(height);
    let rho = Expr::sym(density);
    let g = Expr::sym(gravity);

    // p + (1/2)ρv² + ρgh
    let kinetic_term = Expr::mul(
        Expr::mul(Expr::rational_unchecked(1, 2), rho.clone()),
        Expr::pow(v, Expr::num(2))
    );

    let potential_term = Expr::mul(Expr::mul(rho, g), h);

    simplify(&Expr::add(Expr::add(p, kinetic_term), potential_term))
}

/// Reynolds number: Re = ρvL/μ = vL/ν
///
/// Characterizes flow regime (laminar vs turbulent)
pub fn reynolds_number(
    velocity: &str,
    length_scale: &str,
    kinematic_viscosity: &str,
) -> Expr {
    let v = Expr::sym(velocity);
    let l = Expr::sym(length_scale);
    let nu = Expr::sym(kinematic_viscosity);

    // Re = vL/ν
    simplify(&Expr::mul(
        Expr::mul(v, l),
        Expr::pow(nu, Expr::num(-1))
    ))
}

/// Navier-Stokes momentum equation (symbolic form)
///
/// ∂v/∂t + (v·∇)v = -∇p/ρ + ν∇²v + f
///
/// Returns symbolic representation
pub fn navier_stokes_momentum_symbolic() -> Expr {
    // Symbolic representation: ∂v/∂t + (v·∇)v = -∇p/ρ + ν∇²v + f
    let lhs = Expr::add(
        Expr::func("∂/∂t", vec![Expr::sym("v")]),
        Expr::func("(v·∇)", vec![Expr::sym("v")])
    );

    let pressure_term = Expr::mul(
        Expr::num(-1),
        Expr::mul(
            Expr::func("∇", vec![Expr::sym("p")]),
            Expr::pow(Expr::sym("ρ"), Expr::num(-1))
        )
    );

    let viscosity_term = Expr::mul(
        Expr::sym("ν"),
        Expr::func("∇²", vec![Expr::sym("v")])
    );

    let rhs = Expr::add(
        Expr::add(pressure_term, viscosity_term),
        Expr::sym("f")
    );

    Expr::add(lhs, Expr::mul(Expr::num(-1), rhs))
}

/// Stokes flow (low Reynolds number): ∇p = μ∇²v
///
/// Symbolic representation for creeping flow
pub fn stokes_flow_equation() -> Expr {
    // ∇p = μ∇²v
    Expr::add(
        Expr::func("∇", vec![Expr::sym("p")]),
        Expr::mul(
            Expr::num(-1),
            Expr::mul(
                Expr::sym("μ"),
                Expr::func("∇²", vec![Expr::sym("v")])
            )
        )
    )
}

/// Euler equations (inviscid flow): ∂v/∂t + (v·∇)v = -∇p/ρ + f
pub fn euler_equation_symbolic() -> Expr {
    let lhs = Expr::add(
        Expr::func("∂/∂t", vec![Expr::sym("v")]),
        Expr::func("(v·∇)", vec![Expr::sym("v")])
    );

    let rhs = Expr::add(
        Expr::mul(
            Expr::num(-1),
            Expr::mul(
                Expr::func("∇", vec![Expr::sym("p")]),
                Expr::pow(Expr::sym("ρ"), Expr::num(-1))
            )
        ),
        Expr::sym("f")
    );

    Expr::add(lhs, Expr::mul(Expr::num(-1), rhs))
}

/// Poiseuille flow velocity profile (pipe flow)
///
/// v(r) = (ΔP/4μL)(R² - r²)
pub fn poiseuille_flow_velocity(
    radius_position: &str,
    pipe_radius: &str,
    pressure_drop: &str,
    viscosity: &str,
    length: &str,
) -> Expr {
    let r = Expr::sym(radius_position);
    let big_r = Expr::sym(pipe_radius);
    let delta_p = Expr::sym(pressure_drop);
    let mu = Expr::sym(viscosity);
    let l = Expr::sym(length);

    // (ΔP/4μL)
    let prefactor = Expr::mul(
        delta_p,
        Expr::pow(
            Expr::mul(Expr::mul(Expr::num(4), mu), l),
            Expr::num(-1)
        )
    );

    // (R² - r²)
    let radius_term = Expr::add(
        Expr::pow(big_r, Expr::num(2)),
        Expr::mul(Expr::num(-1), Expr::pow(r, Expr::num(2)))
    );

    simplify(&Expr::mul(prefactor, radius_term))
}

/// Drag force on a sphere (Stokes drag)
///
/// F_d = 6πμRv
pub fn stokes_drag_force(
    viscosity: &str,
    radius: &str,
    velocity: &str,
) -> Expr {
    let mu = Expr::sym(viscosity);
    let r = Expr::sym(radius);
    let v = Expr::sym(velocity);

    // 6πμRv
    simplify(&Expr::mul(
        Expr::mul(
            Expr::mul(Expr::num(6), Expr::sym("π")),
            Expr::mul(mu, r)
        ),
        v
    ))
}

/// Drag coefficient: C_d = F_d / (0.5 ρ v² A)
pub fn drag_coefficient(
    drag_force: &str,
    density: &str,
    velocity: &str,
    area: &str,
) -> Expr {
    let f_d = Expr::sym(drag_force);
    let rho = Expr::sym(density);
    let v = Expr::sym(velocity);
    let a = Expr::sym(area);

    // C_d = F_d / (0.5 ρ v² A)
    let denominator = Expr::mul(
        Expr::mul(Expr::rational_unchecked(1, 2), rho),
        Expr::mul(Expr::pow(v, Expr::num(2)), a)
    );

    simplify(&Expr::mul(f_d, Expr::pow(denominator, Expr::num(-1))))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vorticity_2d() {
        let omega = vorticity_2d("u", "v");
        println!("Vorticity ω_z = {}", omega);
    }

    #[test]
    fn test_stream_function() {
        let (u, v) = stream_function_velocity("ψ");
        println!("u = {}", u);
        println!("v = {}", v);
    }

    #[test]
    fn test_bernoulli() {
        let b = bernoulli_equation("p", "v", "h", "ρ", "g");
        println!("Bernoulli: {}", b);
    }

    #[test]
    fn test_reynolds_number() {
        let re = reynolds_number("v", "L", "ν");
        println!("Reynolds number: {}", re);
    }

    #[test]
    fn test_poiseuille_flow() {
        let v = poiseuille_flow_velocity("r", "R", "ΔP", "μ", "L");
        println!("Poiseuille velocity profile: {}", v);
    }

    #[test]
    fn test_stokes_drag() {
        let f = stokes_drag_force("μ", "R", "v");
        println!("Stokes drag: {}", f);
    }

    #[test]
    fn test_navier_stokes_symbolic() {
        let ns = navier_stokes_momentum_symbolic();
        println!("Navier-Stokes: {}", ns);
    }
}
