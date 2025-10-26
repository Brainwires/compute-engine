// Unit tests for materials_science::crystal_structures
use computational_engine::materials_science::crystal_structures::*;

use super::*;

    #[test]
    fn test_cubic_volume() {
        let cell = UnitCell::cubic(4.0);
        let volume = unit_cell_volume(&cell);
        assert!((volume - 64.0).abs() < 0.001);
    }

    #[test]
    fn test_interplanar_spacing() {
        let d = interplanar_spacing_cubic(4.0, 1, 0, 0);
        assert!((d - 4.0).abs() < 0.001);

        let d = interplanar_spacing_cubic(4.0, 1, 1, 0);
        assert!((d - 2.828).abs() < 0.01);

        let d = interplanar_spacing_cubic(4.0, 1, 1, 1);
        assert!((d - 2.309).abs() < 0.01);
    }

    #[test]
    fn test_bragg_angle() {
        // Cu Kα wavelength = 1.5418 Å
        let theta = bragg_angle(2.0, 1.5418, 1);
        assert!(theta.is_some());
        let angle = theta.unwrap();
        assert!(angle > 0.0 && angle < 90.0);
    }

    #[test]
    fn test_packing_factors() {
        let apf_sc = atomic_packing_factor(BravaisLattice::SimpleCubic);
        assert!((apf_sc - 0.524).abs() < 0.001);

        let apf_bcc = atomic_packing_factor(BravaisLattice::BodyCenteredCubic);
        assert!((apf_bcc - 0.680).abs() < 0.001);

        let apf_fcc = atomic_packing_factor(BravaisLattice::FaceCenteredCubic);
        assert!((apf_fcc - 0.740).abs() < 0.001);
    }

    #[test]
    fn test_systematic_absences() {
        // BCC: only h+k+l even allowed
        assert!(is_reflection_allowed_bcc(1, 1, 0)); // 2 - even
        assert!(!is_reflection_allowed_bcc(1, 0, 0)); // 1 - odd

        // FCC: all even or all odd
        assert!(is_reflection_allowed_fcc(2, 0, 0)); // all even
        assert!(is_reflection_allowed_fcc(1, 1, 1)); // all odd
        assert!(!is_reflection_allowed_fcc(1, 0, 0)); // mixed
    }

    #[test]
    fn test_crystal_density() {
        // Silicon: FCC (diamond cubic), a=5.43 Å, M=28.09 g/mol, 8 atoms/cell for diamond structure
        let vol = 5.43_f64.powi(3);
        let density = crystal_density(8.0, 28.09, vol);
        assert!((density - 2.33).abs() < 0.1); // Silicon density ~2.33 g/cm³
    }
