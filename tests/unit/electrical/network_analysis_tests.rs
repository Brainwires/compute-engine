// Unit tests for electrical::network_analysis
use computational_engine::electrical::network_analysis::*;

use super::*;
    use nalgebra::{dmatrix, dvector};

    #[test]
    fn test_thevenin_norton_conversion() {
        let th = TheveninEquivalent {
            v_th: 12.0,
            r_th: 6.0,
        };

        let norton = thevenin_to_norton(&th);
        assert_eq!(norton.i_n, 2.0);
        assert_eq!(norton.r_n, 6.0);

        let th2 = norton_to_thevenin(&norton);
        assert_eq!(th2.v_th, 12.0);
        assert_eq!(th2.r_th, 6.0);
    }

    #[test]
    fn test_maximum_power_transfer() {
        let p_max = maximum_power_transfer(12.0, 6.0);
        assert_eq!(p_max, 6.0); // 144 / 24 = 6W

        let r_opt = optimal_load_resistance(6.0);
        assert_eq!(r_opt, 6.0);
    }

    #[test]
    fn test_mesh_analysis() {
        // Simple two-mesh circuit
        let r = dmatrix![
            5.0, -1.0;
            -1.0, 4.0
        ];
        let v = dvector![10.0, 5.0];
        let i = mesh_analysis(&r, &v);

        // Check that currents are reasonable
        assert!(i[0] > 0.0 && i[0] < 5.0);
        assert!(i[1] > 0.0 && i[1] < 5.0);
    }

    #[test]
    fn test_delta_wye_conversion() {
        let (r_a, r_b, r_c) = delta_to_wye(30.0, 30.0, 30.0);
        assert_eq!(r_a, 10.0);
        assert_eq!(r_b, 10.0);
        assert_eq!(r_c, 10.0);

        let (r_ab, r_bc, r_ca) = wye_to_delta(10.0, 10.0, 10.0);
        assert_eq!(r_ab, 30.0);
        assert_eq!(r_bc, 30.0);
        assert_eq!(r_ca, 30.0);
    }
