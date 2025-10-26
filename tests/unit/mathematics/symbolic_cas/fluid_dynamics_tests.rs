// Unit tests for mathematics::symbolic_cas::fluid_dynamics
use computational_engine::mathematics::symbolic_cas::fluid_dynamics::*;

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
