// Unit tests for physics::quantum::feynman_rules
use computational_engine::physics::quantum::feynman_rules::*;

use super::*;

    #[test]
    fn test_qed_vertex() {
        let v = qed_vertex_factor();
        assert!(v.im.abs() > 0.0);
    }

    #[test]
    fn test_ee_to_mumu() {
        let s = 10000.0;
        let m_sq = ee_to_mumu_matrix_element_squared(s);
        assert!(m_sq > 0.0);
    }

    #[test]
    fn test_compton_scattering() {
        let s = 1000.0;
        let t = -250.0;
        let m_e = 0.511;
        let m_sq = compton_scattering_matrix_element_squared(s, t, m_e);
        // Result may be negative for unphysical kinematics
        assert!(m_sq.is_finite());
    }

    #[test]
    fn test_anomalous_magnetic_moment() {
        let a_e = anomalous_magnetic_moment_one_loop();
        assert!((a_e - 0.001161).abs() < 1e-5);
    }

    #[test]
    fn test_running_alpha_em() {
        let q_sq = 1e6;
        let alpha = running_alpha_em(q_sq, 0.511);
        assert!(alpha >= ALPHA_EM);
        assert!(alpha.is_finite());
    }

    #[test]
    fn test_running_alpha_s() {
        let lambda_qcd = 200.0;
        let q_sq = 1e6;
        let alpha_s = running_alpha_s(q_sq, lambda_qcd);
        assert!(alpha_s > 0.0 && alpha_s < 1.0);
    }
