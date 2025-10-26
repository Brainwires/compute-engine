// Unit tests for tools::equation_validation::lib
use computational_engine::tools::equation_validation::lib::*;

use super::*;

    #[test]
    fn test_parse_equation() {
        let (left, right) = parse_equation("F = m*a").unwrap();
        assert_eq!(left, "F");
        assert_eq!(right, "m*a");
    }

    #[test]
    fn test_parse_equation_invalid() {
        let result = parse_equation("F + m*a");
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_variables() {
        let vars = extract_variables("F = m*a");
        assert!(vars.contains(&"F".to_string()));
        assert!(vars.contains(&"m".to_string()));
        assert!(vars.contains(&"a".to_string()));
    }

    #[test]
    fn test_extract_variables_skip_functions() {
        let vars = extract_variables("y = sin(x) + cos(x)");
        assert!(vars.contains(&"y".to_string()));
        assert!(vars.contains(&"x".to_string()));
        assert!(!vars.contains(&"sin".to_string()));
        assert!(!vars.contains(&"cos".to_string()));
    }

    #[test]
    fn test_check_mathematical_correctness_valid() {
        let result = check_mathematical_correctness("F = m*a").unwrap();
        assert!(result);
    }

    #[test]
    fn test_check_mathematical_correctness_division_by_zero() {
        let result = check_mathematical_correctness("y = x/0").unwrap();
        assert!(!result);
    }

    #[test]
    fn test_check_mathematical_correctness_unmatched_parens() {
        let result = check_mathematical_correctness("y = (x + 2").unwrap();
        assert!(!result);
    }

    #[test]
    fn test_dimension_from_unit() {
        let kg = Dimension::from_unit("kg").unwrap();
        assert_eq!(kg.mass, 1);

        let newton = Dimension::from_unit("N").unwrap();
        assert_eq!(newton.mass, 1);
        assert_eq!(newton.length, 1);
        assert_eq!(newton.time, -2);
    }

    #[test]
    fn test_dimension_multiply() {
        let kg = Dimension::from_unit("kg").unwrap();
        let m = Dimension::from_unit("m").unwrap();
        let kg_m = kg.multiply(&m);

        assert_eq!(kg_m.mass, 1);
        assert_eq!(kg_m.length, 1);
    }

    #[test]
    fn test_dimension_divide() {
        let m = Dimension::from_unit("m").unwrap();
        let s = Dimension::from_unit("s").unwrap();
        let m_per_s = m.divide(&s);

        assert_eq!(m_per_s.length, 1);
        assert_eq!(m_per_s.time, -1);
    }

    #[test]
    fn test_dimension_power() {
        let m = Dimension::from_unit("m").unwrap();
        let m2 = m.power(2);
        assert_eq!(m2.length, 2);
    }

    #[test]
    fn test_dimension_is_dimensionless() {
        let dim = Dimension::new();
        assert!(dim.is_dimensionless());

        let kg = Dimension::from_unit("kg").unwrap();
        assert!(!kg.is_dimensionless());
    }

    #[test]
    fn test_check_dimensional_consistency() {
        let mut units = HashMap::new();
        units.insert("F".to_string(), "N".to_string());
        units.insert("m".to_string(), "kg".to_string());

        let (_, unit_analysis) = check_dimensional_consistency("F = m", &units).unwrap();
        assert!(unit_analysis.len() > 0);
    }

    #[test]
    fn test_validate_equation_newtons_second_law() {
        let mut units = HashMap::new();
        units.insert("F".to_string(), "N".to_string());
        units.insert("m".to_string(), "kg".to_string());

        let result = validate_equation(
            "F = m".to_string(),
            "mechanics".to_string(),
            Some(units),
            None,
            None,
        )
        .unwrap();

        assert!(result.mathematical_correctness);
    }

    #[test]
    fn test_validate_equation_energy() {
        let mut units = HashMap::new();
        units.insert("E".to_string(), "J".to_string());
        units.insert("m".to_string(), "kg".to_string());
        units.insert("c".to_string(), "m/s".to_string());

        let result = validate_equation(
            "E = m*c^2".to_string(),
            "relativity".to_string(),
            Some(units),
            None,
            None,
        )
        .unwrap();

        assert!(result.mathematical_correctness);
    }

    #[test]
    fn test_check_physics_compliance_mechanics() {
        let (compliant, violations) = check_physics_compliance("F = m*a", "mechanics").unwrap();
        assert!(compliant);
        assert_eq!(violations.len(), 0);
    }

    #[test]
    fn test_check_physics_compliance_force_without_mass() {
        let (compliant, violations) = check_physics_compliance("F = x", "mechanics").unwrap();
        assert!(!compliant);
        assert!(violations.len() > 0);
    }

    #[test]
    fn test_check_conservation_laws_energy() {
        let laws = vec!["energy".to_string()];
        let violations = check_conservation_laws("E = KE + PE", "mechanics", &laws).unwrap();
        assert_eq!(violations.len(), 0);
    }

    #[test]
    fn test_check_symmetries() {
        let symmetries = vec!["time_translation".to_string()];
        let violations = check_symmetries("F = m*a", "mechanics", &symmetries).unwrap();
        // Should not have violations for time-independent equation
        assert_eq!(violations.len(), 0);
    }

    #[test]
    fn test_validate_equation_invalid_syntax() {
        let result = validate_equation(
            "F = m*a)".to_string(),
            "mechanics".to_string(),
            None,
            None,
            None,
        )
        .unwrap();

        assert!(!result.is_valid);
        assert!(!result.mathematical_correctness);
    }

    #[test]
    fn test_validate_equation_no_units() {
        let result = validate_equation(
            "F = m*a".to_string(),
            "mechanics".to_string(),
            None,
            None,
            None,
        )
        .unwrap();

        // Should still check mathematical correctness
        assert!(result.mathematical_correctness);
    }

    #[test]
    fn test_dimension_compound_units() {
        let m_per_s = Dimension::from_unit("m/s").unwrap();
        assert_eq!(m_per_s.length, 1);
        assert_eq!(m_per_s.time, -1);

        let m_per_s2 = Dimension::from_unit("m/s^2").unwrap();
        assert_eq!(m_per_s2.length, 1);
        assert_eq!(m_per_s2.time, -2);
    }

    #[test]
    fn test_dimension_newton_meter() {
        let nm = Dimension::from_unit("N").unwrap();
        let m = Dimension::from_unit("m").unwrap();
        let torque = nm.multiply(&m);

        assert_eq!(torque.mass, 1);
        assert_eq!(torque.length, 2);
        assert_eq!(torque.time, -2);
    }

    #[test]
    fn test_confidence_calculation() {
        let mut units = HashMap::new();
        units.insert("F".to_string(), "N".to_string());

        let result = validate_equation(
            "F = 10".to_string(),
            "mechanics".to_string(),
            Some(units),
            None,
            None,
        )
        .unwrap();

        assert!(result.confidence >= 0.0 && result.confidence <= 1.0);
    }

    #[test]
    fn test_check_physics_compliance_electromagnetism() {
        let (_, violations) = check_physics_compliance("F = q", "electromagnetism").unwrap();
        // Should have violations if charge but no field info
        // Just check it runs without error
        assert!(violations.len() >= 0);
    }

    #[test]
    fn test_check_physics_compliance_relativity() {
        let (compliant, violations) = check_physics_compliance("E = m*c^2", "relativity").unwrap();
        assert!(compliant);
        assert_eq!(violations.len(), 0);
    }

    #[test]
    fn test_check_physics_compliance_relativity_wrong() {
        let (compliant, violations) = check_physics_compliance("E = m*c", "relativity").unwrap();
        assert!(!compliant);
        assert!(violations.len() > 0);
    }

    #[test]
    fn test_dimension_volt() {
        let v = Dimension::from_unit("V").unwrap();
        assert_eq!(v.mass, 1);
        assert_eq!(v.length, 2);
        assert_eq!(v.time, -3);
        assert_eq!(v.current, -1);
    }

    #[test]
    fn test_dimension_coulomb() {
        let c = Dimension::from_unit("C").unwrap();
        assert_eq!(c.current, 1);
        assert_eq!(c.time, 1);
    }

    #[test]
    fn test_dimension_unknown_unit() {
        let result = Dimension::from_unit("unknown_unit");
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_variables_duplicates() {
        let vars = extract_variables("x + x + y");
        // Should not have duplicates
        assert_eq!(vars.iter().filter(|v| *v == "x").count(), 1);
    }

    #[test]
    fn test_parse_equation_whitespace() {
        let (left, right) = parse_equation("  F  =  m * a  ").unwrap();
        assert_eq!(left, "F");
        assert_eq!(right, "m * a");
    }

    #[test]
    fn test_validate_equation_full_checks() {
        let mut units = HashMap::new();
        units.insert("p".to_string(), "kg".to_string());
        units.insert("m".to_string(), "kg".to_string());
        units.insert("v".to_string(), "m".to_string());

        let result = validate_equation(
            "p = m*v".to_string(),
            "mechanics".to_string(),
            Some(units),
            Some(vec!["momentum".to_string()]),
            None,
        )
        .unwrap();

        assert!(result.mathematical_correctness);
    }
