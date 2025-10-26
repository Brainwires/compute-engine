use super::{calculate_biology, BiologyInput, BiologyOperation, BiologyParams};

#[test]
fn test_michaelis_menten_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::MichaelisMenten,
        parameters: BiologyParams {
            vmax: Some(100.0),
            km: Some(10.0),
            substrate_concentration: Some(5.0),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value > 0.0);
    assert_eq!(res.unit, "μmol/(min·mg enzyme)");
}

#[test]
fn test_lineweaver_burk_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::LineweaverBurk,
        parameters: BiologyParams {
            vmax: Some(100.0),
            km: Some(10.0),
            substrate_concentration: Some(20.0),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value > 0.0);
    assert_eq!(res.unit, "1/v");
}

#[test]
fn test_pharmacokinetics_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::Pharmacokinetics,
        parameters: BiologyParams {
            dose: Some(500.0),
            volume_distribution: Some(50.0),
            bioavailability: Some(1.0),
            elimination_rate: Some(0.1),
            time: Some(5.0),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value > 0.0);
    assert_eq!(res.unit, "mg/L");
    assert!(res.additional_data.is_some());
}

#[test]
fn test_hardy_weinberg_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::HardyWeinberg,
        parameters: BiologyParams {
            allele_frequency_p: Some(0.6),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value >= 0.0 && res.value <= 1.0);
    assert_eq!(res.unit, "frequency");
    assert!(res.additional_data.is_some());
}

#[test]
fn test_goldman_equation_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::GoldmanEquation,
        parameters: BiologyParams {
            ion_concentrations_inside: Some(vec![140.0, 12.0, 4.0]),
            ion_concentrations_outside: Some(vec![5.0, 145.0, 116.0]),
            permeabilities: Some(vec![1.0, 0.04, 0.45]),
            temperature: Some(310.0),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert_eq!(res.unit, "mV");
}

#[test]
fn test_allometric_metabolic_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::AllometricScaling,
        parameters: BiologyParams {
            body_mass: Some(70.0),
            scaling_type: Some("metabolic".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value > 0.0);
    assert_eq!(res.unit, "kcal/day");
}

#[test]
fn test_allometric_surface_area_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::AllometricScaling,
        parameters: BiologyParams {
            body_mass: Some(70.0),
            scaling_type: Some("surface_area".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value > 0.0);
    assert_eq!(res.unit, "m²");
}

#[test]
fn test_allometric_lifespan_basic() {
    let input = BiologyInput {
        operation: BiologyOperation::AllometricScaling,
        parameters: BiologyParams {
            body_mass: Some(5.0),
            scaling_type: Some("lifespan".to_string()),
            ..Default::default()
        },
    };

    let result = calculate_biology(input);
    assert!(result.is_ok());
    let res = result.unwrap();
    assert!(res.value > 0.0);
    assert_eq!(res.unit, "years");
}

#[test]
fn test_biology_params_default() {
    let params = BiologyParams::default();
    assert!(params.vmax.is_none());
    assert!(params.km.is_none());
    assert!(params.substrate_concentration.is_none());
    assert!(params.dose.is_none());
    assert!(params.allele_frequency_p.is_none());
}
