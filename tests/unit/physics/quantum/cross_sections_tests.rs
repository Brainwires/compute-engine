// Unit tests for physics::quantum::cross_sections
use computational_engine::physics::quantum::cross_sections::*;

use super::*;

    #[test]
    fn test_ee_to_mumu() {
        let s = 10000.0; // MeV²
        let m_mu = 105.66; // MeV
        let sigma = ee_to_mumu_cross_section(s, m_mu);
        assert!(sigma >= 0.0);
        assert!(sigma.is_finite());
    }

    #[test]
    fn test_compton_thomson_limit() {
        let omega = 0.001; // MeV
        let m_e = 0.511; // MeV
        let sigma = compton_scattering_cross_section(omega, m_e);
        assert!(sigma > 0.0);
        assert!(sigma.is_finite());
    }

    #[test]
    fn test_rutherford_scattering() {
        let theta = PI / 4.0;
        let dsigma = rutherford_scattering_differential(theta, 1.0, 79.0, 5.0);
        assert!(dsigma > 0.0);
    }

    #[test]
    fn test_breit_wigner() {
        let m_z = 91187.6;
        let gamma_z = 2495.0;
        let peak_sigma = 1.0;

        let s_on = m_z * m_z;
        let sigma_on = breit_wigner_cross_section(s_on, m_z, gamma_z, peak_sigma);

        let s_off = (m_z + 10000.0).powi(2);
        let sigma_off = breit_wigner_cross_section(s_off, m_z, gamma_z, peak_sigma);

        assert!(sigma_on > sigma_off);
        assert!(sigma_on.is_finite());
    }

    #[test]
    fn test_event_rate() {
        let sigma = 1.0;
        let lumi = 1e34;
        let rate = event_rate(sigma, lumi);
        assert!(rate > 1e9 && rate < 1e11);
    }
