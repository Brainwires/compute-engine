// Unit tests for materials_science::mechanical_properties
use computational_engine::materials_science::mechanical_properties::*;

use super::*;

    #[test]
    fn test_stress_strain() {
        let stress = engineering_stress(1000.0, 0.001); // 1000N on 1mm²
        assert_eq!(stress, 1e6); // 1 MPa

        let strain = engineering_strain(100.0, 101.0);
        assert_eq!(strain, 0.01); // 1% strain
    }

    #[test]
    fn test_elastic_constants() {
        let e = youngs_from_shear_poisson(80e9, 0.3);
        assert!((e - 208e9).abs() < 1e9);

        let g = shear_from_youngs_poisson(200e9, 0.3);
        assert!((g - 76.9e9).abs() < 1e9);
    }

    #[test]
    fn test_von_mises() {
        // Uniaxial tension: σ₁=σ, σ₂=σ₃=0
        let sigma_vm = von_mises_stress(100e6, 0.0, 0.0);
        assert!((sigma_vm - 100e6).abs() < 1e3);
    }

    #[test]
    fn test_stress_concentration() {
        let kt_circle = stress_concentration_circular();
        assert_eq!(kt_circle, 3.0);

        let kt_ellipse = stress_concentration_ellipse(10.0, 1.0);
        assert_eq!(kt_ellipse, 21.0); // 1 + 2*10
    }

    #[test]
    fn test_hardness_conversion() {
        let bhn = rockwell_c_to_brinell(50.0);
        assert!((bhn - 950.0).abs() < 50.0); // Approximate
    }
