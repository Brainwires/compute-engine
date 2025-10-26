// Unit tests for physics::quantum::particle_decays
use computational_engine::physics::quantum::particle_decays::*;

use super::*;

    #[test]
    fn test_muon_decay_width() {
        let m_mu = 105.66;
        let width = muon_decay_width(m_mu);
        assert!(width > 0.0);
        assert!(width.is_finite());
    }

    #[test]
    fn test_muon_lifetime() {
        let m_mu = 105.66;
        let width = muon_decay_width(m_mu);
        let tau = lifetime_from_width(width);
        assert!(tau > 1e-6 && tau < 3e-6);
    }

    #[test]
    fn test_w_boson_width() {
        let m_w = 80379.0;
        let width = w_boson_decay_width(m_w);
        assert!(width > 0.0);
        assert!(width.is_finite());
    }

    #[test]
    fn test_z_boson_width() {
        let m_z = 91187.6;
        let width = z_boson_decay_width(m_z);
        assert!(width > 0.0);
        assert!(width.is_finite());
    }

    #[test]
    fn test_higgs_to_bottom() {
        let m_h = 125090.0;
        let m_b = 4180.0;
        let n_c = 3.0;
        let width = higgs_to_fermion_width(m_h, m_b, n_c);
        assert!(width > 1000.0);
    }

    #[test]
    fn test_top_quark_width() {
        let m_t = 173070.0;
        let m_w = 80379.0;
        let width = top_quark_decay_width(m_t, m_w);
        assert!(width > 0.0);
        assert!(width.is_finite());
    }

    #[test]
    fn test_branching_ratios() {
        let channels = vec![("a", 1.0), ("b", 2.0), ("c", 1.0)];
        let brs = calculate_branching_ratios(channels);
        let total: f64 = brs.iter().map(|br| br.branching_ratio).sum();
        assert!((total - 1.0).abs() < 1e-10);
    }
