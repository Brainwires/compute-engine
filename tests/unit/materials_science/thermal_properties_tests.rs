// Unit tests for materials_science::thermal_properties
use computational_engine::materials_science::thermal_properties::*;

use super::*;

    #[test]
    fn test_dulong_petit() {
        // For 1 mole of atoms (N_A atoms)
        let cv = dulong_petit_heat_capacity(N_A);
        let expected = 3.0 * R; // 3R ≈ 24.9 J/(mol·K)
        assert!((cv - expected).abs() < 0.1);
    }

    #[test]
    fn test_thermal_expansion() {
        // Steel: α ≈ 12×10⁻⁶ /K
        let alpha = 12e-6;
        let original_length = 1.0; // 1 meter
        let delta_t = 100.0; // 100K temperature rise

        let delta_l = thermal_expansion_length(original_length, alpha, delta_t);
        assert!((delta_l - 0.0012).abs() < 1e-6); // 1.2 mm expansion

        let new_length = length_after_expansion(original_length, alpha, delta_t);
        assert!((new_length - 1.0012).abs() < 1e-6);
    }

    #[test]
    fn test_thermal_stress() {
        // Steel constrained expansion
        let e = 200e9; // 200 GPa
        let alpha = 12e-6; // /K
        let delta_t = 100.0; // K

        let stress = thermal_stress_constrained(e, alpha, delta_t);
        assert!((stress - 240e6).abs() < 1e6); // 240 MPa
    }

    #[test]
    fn test_thermal_diffusivity() {
        // Copper: κ ≈ 400 W/(m·K), ρ ≈ 8960 kg/m³, C_p ≈ 385 J/(kg·K)
        let kappa = 400.0;
        let rho = 8960.0;
        let cp = 385.0;

        let alpha = thermal_diffusivity(kappa, rho, cp);
        assert!(alpha > 1e-4 && alpha < 2e-4); // ~1.16×10⁻⁴ m²/s
    }

    #[test]
    fn test_wiedemann_franz() {
        // Copper at 300K: σ ≈ 6×10⁷ S/m
        let sigma = 6e7;
        let temp = 300.0;

        let kappa = thermal_conductivity_wiedemann_franz(sigma, temp);
        assert!(kappa > 400.0 && kappa < 500.0); // Should be ~440 W/(m·K)
    }

    #[test]
    fn test_heat_transfer_slab() {
        // 1 W/(m·K) material, 1 m² area, 1 m thick, 100K difference
        let q = heat_transfer_slab(1.0, 1.0, 400.0, 300.0, 1.0);
        assert!((q - 100.0).abs() < 0.01); // 100 W
    }

    #[test]
    fn test_biot_number() {
        // h = 100 W/(m²·K), L = 0.01 m, κ = 50 W/(m·K)
        let bi = biot_number(100.0, 0.01, 50.0);
        assert!((bi - 0.02).abs() < 0.001); // Bi << 1 (lumped capacitance valid)
    }

    #[test]
    fn test_stefan_boltzmann() {
        // Black body (ε=1) at 1000K, 1 m², surrounding at 300K
        let q = stefan_boltzmann_radiation(1.0, 1.0, 1000.0, 300.0);
        assert!(q > 50000.0 && q < 60000.0); // ~56 kW
    }
