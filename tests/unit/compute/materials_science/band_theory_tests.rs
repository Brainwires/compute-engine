// Unit tests for materials_science::band_theory
use computational_engine::materials_science::band_theory::*;

use super::*;

    #[test]
    fn test_fermi_energy() {
        // Typical metal: n ~ 10²⁸ m⁻³
        let ef = fermi_energy_free_electron(1e28);
        assert!(ef > 0.0 && ef < 20.0); // Typical range 2-10 eV
    }

    #[test]
    fn test_material_classification() {
        assert_eq!(classify_material(0.0), MaterialClass::Conductor);
        assert_eq!(classify_material(1.1), MaterialClass::Semiconductor);
        assert_eq!(classify_material(5.0), MaterialClass::Insulator);
    }

    #[test]
    fn test_intrinsic_carriers() {
        // Silicon at 300K: E_g = 1.12 eV
        let ni = intrinsic_carrier_concentration(1.12, 300.0);
        assert!(ni > 1e15 && ni < 1e16); // ~10¹⁶ m⁻³ for Si
    }

    #[test]
    fn test_plasma_frequency() {
        let omega_p = plasma_frequency(1e28);
        let freq_hz = omega_p / (2.0 * PI);
        assert!(freq_hz > 8e14 && freq_hz < 1e15); // ~9×10^14 Hz
    }

    #[test]
    fn test_thermionic_emission() {
        // Tungsten: φ ≈ 4.5 eV at 2500K
        let j = thermionic_emission_current(2500.0, 4.5);
        assert!(j > 6e3 && j < 7e3); // ~6370 A/m^2
    }
