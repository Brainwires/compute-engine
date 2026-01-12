// Unit tests for materials_science::diffusion
use computational_engine::materials_science::diffusion::*;

use super::*;

    #[test]
    fn test_diffusion_distance() {
        // Typical metal diffusion: D ~ 10⁻¹⁴ m²/s at moderate temp
        // After 1 hour (3600 s)
        let d = 1e-14;
        let t = 3600.0;
        let x = diffusion_distance(d, t);

        assert!(x > 1e-6 && x < 1e-5); // Should be a few micrometers
    }

    #[test]
    fn test_arrhenius_diffusion() {
        // Typical values for substitutional diffusion in metals
        let d0 = 1e-4; // m²/s
        let q = 200e3; // J/mol (200 kJ/mol)
        let temp = 1000.0; // K

        let d = arrhenius_diffusion_coefficient(d0, q, temp);
        assert!(d > 0.0 && d < d0); // Should be less than D₀
    }

    #[test]
    fn test_diffusion_time() {
        // How long to diffuse 1 μm with D = 10⁻¹⁴ m²/s?
        let d = 1e-14;
        let x = 1e-6;
        let t = diffusion_time(d, x);

        assert!(t > 40.0 && t < 60.0); // ~50 seconds
    }

    #[test]
    fn test_membrane_flux() {
        // D = 10⁻¹⁰ m²/s, ΔC = 100 mol/m³, L = 100 μm
        let j = membrane_flux_steady_state(1e-10, 100.0, 0.0, 100e-6);
        assert!((j - 1e-4).abs() < 1e-6); // 10⁻⁴ mol/(m²·s)
    }

    #[test]
    fn test_interdiffusion() {
        // 50-50 mixture with different diffusion coefficients
        let xa = 0.5;
        let xb = 0.5;
        let da = 1e-14;
        let db = 2e-14;

        let d_int = interdiffusion_coefficient(xa, xb, da, db);
        assert!((d_int - 1.5e-14).abs() < 1e-16); // Average
    }

    #[test]
    fn test_einstein_diffusion() {
        // Water at 298K: η ≈ 0.001 Pa·s, r = 1 nm
        let d = einstein_diffusion_coefficient(298.0, 0.001, 1e-9);
        assert!(d > 1e-10 && d < 1e-9); // Should be ~10⁻⁹ m²/s
    }

    #[test]
    fn test_concentration_semi_infinite() {
        // At t=0, any position x>0 should have C ≈ 0
        // At large time, concentration should approach C₀ everywhere
        let c0 = 100.0;
        let d = 1e-14;

        // Short time, far distance -> low concentration
        let c1 = concentration_semi_infinite(c0, 1e-4, d, 1.0);
        assert!(c1 < c0 * 0.1);

        // Long time, same distance -> higher concentration
        let c2 = concentration_semi_infinite(c0, 1e-4, d, 1e8);
        assert!(c2 > c0 * 0.5);
    }

    #[test]
    fn test_penetration_depth() {
        // For D = 10⁻¹⁴ m²/s, t = 1 hour
        let depth = penetration_depth(1e-14, 3600.0, 0.5);
        assert!(depth > 1e-5 && depth < 1.5e-5); // ~14 micrometers
    }
