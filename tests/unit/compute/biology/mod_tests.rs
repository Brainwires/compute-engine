// Unit tests for biology::mod
use computational_engine::biology::mod::*;

use super::*;

    // Michaelis-Menten Tests (3 tests)
    #[test]
    fn test_michaelis_menten_at_km() {
        let params = BiologyParams {
            vmax: Some(100.0),
            km: Some(10.0),
            substrate_concentration: Some(10.0),
            ..Default::default()
        };
        let result = calculate_michaelis_menten(&params).unwrap();
        assert!((result.value - 50.0).abs() < 0.1); // At Km, v = Vmax/2
    }

    #[test]
    fn test_michaelis_menten_saturated() {
        let params = BiologyParams {
            vmax: Some(100.0),
            km: Some(1.0),
            substrate_concentration: Some(100.0), // >> Km
            ..Default::default()
        };
        let result = calculate_michaelis_menten(&params).unwrap();
        assert!(result.value > 99.0); // Nearly Vmax
        assert!(result.interpretation.contains("saturated"));
    }

    #[test]
    fn test_michaelis_menten_low_substrate() {
        let params = BiologyParams {
            vmax: Some(100.0),
            km: Some(100.0),
            substrate_concentration: Some(10.0), // << Km
            ..Default::default()
        };
        let result = calculate_michaelis_menten(&params).unwrap();
        assert!(result.value < 10.0); // Much less than Vmax
    }

    // Lineweaver-Burk Tests (2 tests)
    #[test]
    fn test_lineweaver_burk_calculation() {
        let params = BiologyParams {
            vmax: Some(100.0),
            km: Some(10.0),
            substrate_concentration: Some(20.0),
            ..Default::default()
        };
        let result = calculate_lineweaver_burk(&params).unwrap();
        let additional = result.additional_data.unwrap();
        assert!(additional["slope_km_over_vmax"] > 0.0);
        assert!((additional["y_intercept_1_over_vmax"] - 0.01).abs() < 0.001);
    }

    #[test]
    fn test_lineweaver_burk_high_substrate() {
        let params = BiologyParams {
            vmax: Some(50.0),
            km: Some(5.0),
            substrate_concentration: Some(100.0),
            ..Default::default()
        };
        let result = calculate_lineweaver_burk(&params).unwrap();
        assert!(result.value < 0.1); // 1/v is small when v is high
    }

    // Pharmacokinetics Tests (4 tests)
    #[test]
    fn test_pharmacokinetics_initial() {
        let params = BiologyParams {
            dose: Some(500.0),
            volume_distribution: Some(50.0),
            elimination_rate: Some(0.1),
            time: Some(0.0), // Initial time
            ..Default::default()
        };
        let result = calculate_pharmacokinetics(&params).unwrap();
        assert!((result.value - 10.0).abs() < 0.1); // C(0) = Dose/V = 500/50 = 10
    }

    #[test]
    fn test_pharmacokinetics_decay() {
        let params = BiologyParams {
            dose: Some(100.0),
            volume_distribution: Some(10.0),
            elimination_rate: Some(0.693), // k = 0.693 gives t½ = 1h
            time: Some(1.0),               // One half-life
            ..Default::default()
        };
        let result = calculate_pharmacokinetics(&params).unwrap();
        assert!((result.value - 5.0).abs() < 0.5); // Should be ~half of 10.0
    }

    #[test]
    fn test_pharmacokinetics_bioavailability() {
        let params = BiologyParams {
            dose: Some(100.0),
            volume_distribution: Some(10.0),
            bioavailability: Some(0.5), // 50% bioavailable
            elimination_rate: Some(0.1),
            time: Some(0.0),
            ..Default::default()
        };
        let result = calculate_pharmacokinetics(&params).unwrap();
        assert!((result.value - 5.0).abs() < 0.1); // C = 0.5*100/10 = 5
    }

    #[test]
    fn test_pharmacokinetics_half_life() {
        let params = BiologyParams {
            dose: Some(100.0),
            volume_distribution: Some(10.0),
            elimination_rate: Some(0.693),
            time: Some(0.0),
            ..Default::default()
        };
        let result = calculate_pharmacokinetics(&params).unwrap();
        let additional = result.additional_data.unwrap();
        assert!((additional["half_life_hours"] - 1.0).abs() < 0.1);
    }

    // Hardy-Weinberg Tests (4 tests)
    #[test]
    fn test_hardy_weinberg_sum_equals_one() {
        let params = BiologyParams {
            allele_frequency_p: Some(0.6),
            ..Default::default()
        };
        let result = calculate_hardy_weinberg(&params).unwrap();
        let additional = result.additional_data.unwrap();
        let sum = result.value
            + additional["freq_Aa_heterozygous"]
            + additional["freq_aa_homozygous_recessive"];
        assert!((sum - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_hardy_weinberg_equal_alleles() {
        let params = BiologyParams {
            allele_frequency_p: Some(0.5),
            ..Default::default()
        };
        let result = calculate_hardy_weinberg(&params).unwrap();
        let additional = result.additional_data.unwrap();
        // p=q=0.5: AA=0.25, Aa=0.50, aa=0.25
        assert!((result.value - 0.25).abs() < 0.01);
        assert!((additional["freq_Aa_heterozygous"] - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_hardy_weinberg_dominant_allele() {
        let params = BiologyParams {
            allele_frequency_p: Some(0.9),
            ..Default::default()
        };
        let result = calculate_hardy_weinberg(&params).unwrap();
        assert!((result.value - 0.81).abs() < 0.01); // AA = 0.9^2 = 0.81
    }

    #[test]
    fn test_hardy_weinberg_recessive_allele() {
        let params = BiologyParams {
            allele_frequency_p: Some(0.1),
            ..Default::default()
        };
        let result = calculate_hardy_weinberg(&params).unwrap();
        let additional = result.additional_data.unwrap();
        assert!((additional["freq_aa_homozygous_recessive"] - 0.81).abs() < 0.01); // aa = 0.9^2
    }

    // Goldman Equation Tests (3 tests)
    #[test]
    fn test_goldman_equation_basic() {
        let params = BiologyParams {
            ion_concentrations_inside: Some(vec![140.0, 12.0, 4.0]), // K+, Na+, Cl-
            ion_concentrations_outside: Some(vec![5.0, 145.0, 116.0]),
            permeabilities: Some(vec![1.0, 0.04, 0.45]), // Relative permeabilities
            temperature: Some(310.0),
            ..Default::default()
        };
        let result = calculate_goldman(&params).unwrap();
        // Resting potential should be negative (around -70 mV)
        assert!(result.value < -50.0 && result.value > -90.0);
    }

    #[test]
    fn test_goldman_equation_potassium_dominated() {
        let params = BiologyParams {
            ion_concentrations_inside: Some(vec![140.0, 10.0]),
            ion_concentrations_outside: Some(vec![5.0, 145.0]),
            permeabilities: Some(vec![1.0, 0.0]), // Only K+ permeable
            temperature: Some(310.0),
            ..Default::default()
        };
        let result = calculate_goldman(&params).unwrap();
        // Should be close to Nernst potential for K+ (~ -88 mV)
        assert!(result.value < -80.0);
    }

    #[test]
    fn test_goldman_equation_sodium_dominated() {
        let params = BiologyParams {
            ion_concentrations_inside: Some(vec![140.0, 10.0]),
            ion_concentrations_outside: Some(vec![5.0, 145.0]),
            permeabilities: Some(vec![0.0, 1.0]), // Only Na+ permeable
            temperature: Some(310.0),
            ..Default::default()
        };
        let result = calculate_goldman(&params).unwrap();
        // Should be close to Nernst potential for Na+ (~ +67 mV)
        assert!(result.value > 50.0);
    }

    // Allometric Scaling Tests (3 tests)
    #[test]
    fn test_allometric_metabolic_scaling() {
        let params = BiologyParams {
            body_mass: Some(70.0), // Human
            scaling_type: Some("metabolic".to_string()),
            ..Default::default()
        };
        let result = calculate_allometric(&params).unwrap();
        // BMR for 70kg human should be around 1500-2000 kcal/day
        assert!(result.value > 1000.0 && result.value < 3000.0);
        assert_eq!(result.unit, "kcal/day");
    }

    #[test]
    fn test_allometric_surface_area() {
        let params = BiologyParams {
            body_mass: Some(70.0),
            scaling_type: Some("surface_area".to_string()),
            ..Default::default()
        };
        let result = calculate_allometric(&params).unwrap();
        // BSA for 70kg human should be around 1.7-1.9 m²
        assert!(result.value > 1.5 && result.value < 2.5);
        assert_eq!(result.unit, "m²");
    }

    #[test]
    fn test_allometric_lifespan() {
        let params = BiologyParams {
            body_mass: Some(5000.0), // Elephant
            scaling_type: Some("lifespan".to_string()),
            ..Default::default()
        };
        let result = calculate_allometric(&params).unwrap();
        // Larger animals live longer
        assert!(result.value > 40.0);
    }
