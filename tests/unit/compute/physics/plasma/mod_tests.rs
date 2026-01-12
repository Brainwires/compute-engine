// Unit tests for physics::plasma::mod
use computational_engine::compute::physics::plasma::mod::*;

use super::*;

    #[test]
    fn test_debye_length_tokamak() {
        let params = PlasmaParams::tokamak();
        let lambda_d = debye_length(&params);

        // Tokamak: λ_D ~ 10^-5 m
        assert!(lambda_d > 1e-6 && lambda_d < 1e-4);
    }

    #[test]
    fn test_plasma_frequency() {
        let n_e = 1e20; // Tokamak density
        let omega_pe = plasma_frequency(n_e);

        // ω_pe ~ 10^12 rad/s
        assert!(omega_pe > 1e11 && omega_pe < 1e13);
    }

    #[test]
    fn test_cyclotron_frequencies() {
        let b = 5.0; // 5 Tesla

        let omega_ce = electron_cyclotron_frequency(b);
        let omega_ci = ion_cyclotron_frequency(b, 1.0, M_P);

        // Electron cyclotron frequency much higher than ion
        assert!(omega_ce > omega_ci * 1000.0);
    }

    #[test]
    fn test_thermal_velocity() {
        let t_e = 1000.0; // 1 keV
        let v_th_e = thermal_velocity(t_e, M_E);

        // Electron thermal velocity ~ 10^7 m/s
        assert!(v_th_e > 1e6 && v_th_e < 1e8);
    }

    #[test]
    fn test_larmor_radius() {
        let v_perp = 1e6; // 1 Mm/s
        let omega_c = 1e8; // rad/s

        let r_l = larmor_radius(v_perp, omega_c);

        assert!(r_l > 0.0);
        assert!(r_l.is_finite());
    }

    #[test]
    fn test_plasma_beta() {
        let params = PlasmaParams::tokamak();
        let beta = plasma_beta(&params);

        // Tokamak typically has β < 0.1
        assert!(beta > 0.0);
        assert!(beta < 1.0);
    }

    #[test]
    fn test_alfven_velocity() {
        let b = 5.0; // Tesla
        let rho = 1e-6; // kg/m³ (low density plasma)

        let v_a = alfven_velocity(b, rho);

        // Alfvén velocity ~ 10^6 m/s for typical plasma
        assert!(v_a > 0.0);
        assert!(v_a.is_finite());
    }

    #[test]
    fn test_debye_number() {
        let params = PlasmaParams::tokamak();
        let n_d = debye_number(&params);

        // Should be >> 1 for good plasma
        assert!(n_d > 1000.0);
    }

    #[test]
    fn test_is_plasma() {
        let tokamak = PlasmaParams::tokamak();
        assert!(is_plasma(&tokamak));

        let corona = PlasmaParams::solar_corona();
        assert!(is_plasma(&corona));
    }

    #[test]
    fn test_temperature_conversion() {
        let params = PlasmaParams::tokamak();
        let t_k = params.electron_temp_kelvin();

        // 10 keV ~ 100 million K
        assert!(t_k > 1e8 && t_k < 1e9);
    }
