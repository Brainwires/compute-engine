// Unit tests for physics::quantum::propagators
use computational_engine::physics::quantum::propagators::*;

use super::*;

    #[test]
    fn test_scalar_propagator() {
        let p_sq = 100.0;
        let mass = 5.0;
        let prop = scalar_propagator(p_sq, mass);
        assert!(prop.norm() > 0.0);
    }

    #[test]
    fn test_photon_propagator() {
        let q_sq = 1.0;
        let prop = photon_propagator(q_sq);
        assert!(prop.norm() > 0.0);
    }

    #[test]
    fn test_yukawa_potential() {
        let r = 1.0;
        let mass = 100.0;
        let coupling = 0.1;
        let v = yukawa_potential(r, mass, coupling);
        assert!(v < 0.0);
    }

    #[test]
    fn test_coulomb_potential() {
        let r = 1.0;
        let alpha = ALPHA_EM;
        let v = coulomb_potential(r, alpha);
        assert!(v < 0.0);
    }
