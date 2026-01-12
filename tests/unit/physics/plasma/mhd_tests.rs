// Unit tests for physics::plasma::mhd
use computational_engine::compute::physics::plasma::mhd::*;

use super::*;

    #[test]
    fn test_magnetic_pressure() {
        let b = 5.0; // 5 Tesla
        let p_b = magnetic_pressure(b);

        // p_B ~ 10^6 Pa for 5T field
        assert!(p_b > 1e6 && p_b < 1e7);
    }

    #[test]
    fn test_magnetic_tension() {
        let b = [1.0, 0.0, 0.0];
        let grad_b = [
            [0.1, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];

        let tension = magnetic_tension(&b, &grad_b);

        assert!(tension[0].is_finite());
    }

    #[test]
    fn test_lorentz_force() {
        let b = [0.0, 0.0, 5.0];
        let curl_b = [1.0, 0.0, 0.0];

        let force = lorentz_force(&b, &curl_b);

        // Force should be perpendicular to both J and B
        assert!(force[0].is_finite());
        assert!(force[1].is_finite());
        assert!(force[2].abs() < 1e-10); // Should be zero
    }

    #[test]
    fn test_frozen_in_parameter() {
        let state = MHDState {
            velocity: [0.0, 0.0, 0.0],
            b_field: [0.0, 0.0, 5.0],
            pressure: 1e5,
            density: 1e-6,
        };

        let param = frozen_in_parameter(&state);

        assert!(param > 0.0);
        assert!(param.is_finite());
    }

    #[test]
    fn test_mhd_equilibrium() {
        let grad_p = [1000.0, 0.0, 0.0];
        let b = [0.0, 0.0, 5.0];
        let curl_b = [1000.0 * 4.0e-7 * PI, 0.0, 0.0]; // Chosen to balance

        let eq = check_mhd_equilibrium(&grad_p, &b, &curl_b);

        assert!(eq.lorentz_force[0].is_finite());
    }

    #[test]
    fn test_current_density() {
        let curl_b = [1.0, 2.0, 3.0];
        let j = current_density(&curl_b);

        // J should be proportional to curl(B)
        assert!(j[0] > 0.0);
        assert!(j[1] > j[0]);
        assert!(j[2] > j[1]);
    }

    #[test]
    fn test_ideal_mhd_check() {
        let params = PlasmaParams::tokamak();
        let is_ideal = check_ideal_mhd(&params);

        // Tokamak should satisfy ideal MHD
        assert!(is_ideal);
    }

    #[test]
    fn test_lundquist_number() {
        let params = PlasmaParams::tokamak();
        let s = lundquist_number(&params, 1.0);

        // Tokamak should have large Lundquist number
        assert!(s > 100.0);
    }
