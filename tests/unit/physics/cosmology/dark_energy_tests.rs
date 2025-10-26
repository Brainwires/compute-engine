// Unit tests for physics::cosmology::dark_energy
use computational_engine::physics::cosmology::dark_energy::*;

use super::*;

    #[test]
    fn test_lambda_cdm_eos() {
        let model = DarkEnergyModel::LambdaCDM;

        assert_eq!(model.w(0.0), -1.0);
        assert_eq!(model.w(1.0), -1.0);
        assert_eq!(model.w(1000.0), -1.0);
    }

    #[test]
    fn test_quintessence_eos() {
        let model = DarkEnergyModel::Quintessence { w0: -0.8 };

        let w = model.w(0.0);
        assert!(w > -1.0 && w < -1.0/3.0);
    }

    #[test]
    fn test_phantom_eos() {
        let model = DarkEnergyModel::Phantom { w0: -1.2 };

        let w = model.w(0.0);
        assert!(w < -1.0);
    }

    #[test]
    fn test_cpl_evolution() {
        let model = DarkEnergyModel::CPL { w0: -0.9, wa: -0.1 };

        let w_today = model.w(0.0);
        let w_past = model.w(1.0);

        // w should vary with redshift
        assert!((w_today - (-0.9)).abs() < 0.01);
        assert!(w_past != w_today);
    }

    #[test]
    fn test_lambda_density_constant() {
        let model = DarkEnergyModel::LambdaCDM;

        let rho_today = model.density_evolution(Redshift(0.0));
        let rho_past = model.density_evolution(Redshift(1.0));

        // Cosmological constant has constant density
        assert_eq!(rho_today, 1.0);
        assert_eq!(rho_past, 1.0);
    }

    #[test]
    fn test_quintessence_density_evolution() {
        let model = DarkEnergyModel::Quintessence { w0: -0.8 };

        let rho_today = model.density_evolution(Redshift(0.0));
        let rho_past = model.density_evolution(Redshift(1.0));

        // Quintessence density should increase going back in time
        assert!(rho_past > rho_today);
    }

    #[test]
    fn test_big_rip_phantom() {
        let params = CosmologyParams::planck_2018();

        // Phantom energy (w < -1) leads to Big Rip
        let t_rip = big_rip_time(&params, -1.2);
        assert!(t_rip.is_some());
        assert!(t_rip.unwrap() > 0.0);

        // Lambda CDM (w = -1) has no Big Rip
        let t_rip_lambda = big_rip_time(&params, -1.0);
        assert!(t_rip_lambda.is_none());
    }

    #[test]
    fn test_hubble_parameter_de() {
        let params = CosmologyParams::planck_2018();
        let model = DarkEnergyModel::LambdaCDM;

        let h_today = hubble_parameter_de(Redshift(0.0), &params, &model);
        let h0 = params.hubble_constant_si();

        // Should equal H0 today
        assert!((h_today - h0).abs() / h0 < 0.01);
    }

    #[test]
    fn test_future_evolution() {
        let params = CosmologyParams::planck_2018();
        let model = DarkEnergyModel::LambdaCDM;

        // Simulate 1 billion years into future
        let evolution = future_evolution(3.15e16, 100, &params, &model);

        assert_eq!(evolution.times.len(), 101);
        assert_eq!(evolution.scale_factors.len(), 101);

        // Scale factor should increase
        assert!(evolution.scale_factors.last().unwrap() > &1.0);

        // All values should be finite
        for a in &evolution.scale_factors {
            assert!(a.is_finite());
            assert!(*a > 0.0);
        }
    }

    #[test]
    fn test_accelerating_expansion() {
        let params = CosmologyParams::planck_2018();
        let model = DarkEnergyModel::LambdaCDM;

        let evolution = future_evolution(1e17, 100, &params, &model);

        // Check that expansion is accelerating (da/dt increasing)
        let early_rate = (evolution.scale_factors[10] - evolution.scale_factors[0]) / 10.0;
        let late_rate = (evolution.scale_factors[100] - evolution.scale_factors[90]) / 10.0;

        assert!(late_rate > early_rate);
    }
