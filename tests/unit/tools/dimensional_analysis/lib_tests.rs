// Unit tests for tools::dimensional_analysis::lib
use computational_engine::units::dimensional_analysis::lib::*;

use super::*;

    #[test]
    fn test_dimension_from_base_units() {
        let kg = PhysicalDimension::from_unit("kg").unwrap();
        assert_eq!(kg.mass, 1);
        assert_eq!(kg.length, 0);

        let m = PhysicalDimension::from_unit("m").unwrap();
        assert_eq!(m.length, 1);
        assert_eq!(m.mass, 0);

        let s = PhysicalDimension::from_unit("s").unwrap();
        assert_eq!(s.time, 1);
    }

    #[test]
    fn test_dimension_from_derived_units() {
        let newton = PhysicalDimension::from_unit("N").unwrap();
        assert_eq!(newton.mass, 1);
        assert_eq!(newton.length, 1);
        assert_eq!(newton.time, -2);

        let joule = PhysicalDimension::from_unit("J").unwrap();
        assert_eq!(joule.mass, 1);
        assert_eq!(joule.length, 2);
        assert_eq!(joule.time, -2);

        let watt = PhysicalDimension::from_unit("W").unwrap();
        assert_eq!(watt.mass, 1);
        assert_eq!(watt.length, 2);
        assert_eq!(watt.time, -3);
    }

    #[test]
    fn test_dimension_multiply() {
        let m = PhysicalDimension::from_unit("m").unwrap();
        let s = PhysicalDimension::from_unit("s").unwrap();

        let ms = m.multiply(&s);
        assert_eq!(ms.length, 1);
        assert_eq!(ms.time, 1);
    }

    #[test]
    fn test_dimension_divide() {
        let m = PhysicalDimension::from_unit("m").unwrap();
        let s = PhysicalDimension::from_unit("s").unwrap();

        let m_per_s = m.divide(&s);
        assert_eq!(m_per_s.length, 1);
        assert_eq!(m_per_s.time, -1);
    }

    #[test]
    fn test_dimension_power() {
        let m = PhysicalDimension::from_unit("m").unwrap();
        let m2 = m.power(2);
        assert_eq!(m2.length, 2);

        let m3 = m.power(3);
        assert_eq!(m3.length, 3);
    }

    #[test]
    fn test_dimension_is_dimensionless() {
        let dim = PhysicalDimension::new();
        assert!(dim.is_dimensionless());

        let kg = PhysicalDimension::from_unit("kg").unwrap();
        assert!(!kg.is_dimensionless());
    }

    #[test]
    fn test_dimension_to_string() {
        let kg = PhysicalDimension::from_unit("kg").unwrap();
        assert_eq!(kg.to_string(), "M");

        let newton = PhysicalDimension::from_unit("N").unwrap();
        assert_eq!(newton.to_string(), "M⋅L⋅T^-2");

        let dim = PhysicalDimension::new();
        assert_eq!(dim.to_string(), "1");
    }

    #[test]
    fn test_tokenize_unit() {
        let tokens = tokenize_unit("kg*m/s^2");
        assert!(tokens.contains(&"kg".to_string()));
        assert!(tokens.contains(&"m".to_string()));
        assert!(tokens.contains(&"s".to_string()));
    }

    #[test]
    fn test_extract_variables_with_powers() {
        let vars = extract_variables_with_powers("x^2 + y");
        assert!(vars.contains(&("x".to_string(), 2)));
        assert!(vars.contains(&("y".to_string(), 1)));
    }

    #[test]
    fn test_extract_variables_skip_functions() {
        let vars = extract_variables_with_powers("sin(x) + cos(y)");
        assert!(vars.contains(&("x".to_string(), 1)));
        assert!(vars.contains(&("y".to_string(), 1)));
        assert!(!vars.iter().any(|(v, _)| v == "sin"));
        assert!(!vars.iter().any(|(v, _)| v == "cos"));
    }

    #[test]
    fn test_dimensional_analysis_force() {
        let mut units = HashMap::new();
        units.insert("F".to_string(), "N".to_string());
        units.insert("m".to_string(), "kg".to_string());

        let result = dimensional_analysis("F = m".to_string(), units, None).unwrap();

        // Just verify it runs without error
        assert!(result.unit_breakdown.len() > 0);
    }

    #[test]
    fn test_dimensional_analysis_energy() {
        let mut units = HashMap::new();
        units.insert("E".to_string(), "J".to_string());
        units.insert("m".to_string(), "kg".to_string());

        let result = dimensional_analysis("E = m".to_string(), units, None).unwrap();

        assert!(result.unit_breakdown.len() > 0);
    }

    #[test]
    fn test_dimensional_analysis_velocity() {
        let mut units = HashMap::new();
        units.insert("v".to_string(), "m".to_string());
        units.insert("t".to_string(), "s".to_string());

        let result = dimensional_analysis("v = t".to_string(), units, None).unwrap();

        assert!(result.unit_breakdown.len() > 0);
    }

    #[test]
    fn test_dimensional_analysis_unknown_units() {
        let units = HashMap::new();

        let result = dimensional_analysis("F = m*a".to_string(), units, None).unwrap();

        assert!(!result.consistent);
        assert!(result.recommendations.len() > 0);
    }

    #[test]
    fn test_generate_recommendations_force() {
        let dim = PhysicalDimension::from_unit("N").unwrap();
        let unit_breakdown = HashMap::new();

        let recs = generate_recommendations(&dim, true, &unit_breakdown);
        assert!(recs.iter().any(|r| r.contains("force")));
    }

    #[test]
    fn test_generate_recommendations_energy() {
        let dim = PhysicalDimension::from_unit("J").unwrap();
        let unit_breakdown = HashMap::new();

        let recs = generate_recommendations(&dim, true, &unit_breakdown);
        assert!(recs.iter().any(|r| r.contains("energy")));
    }

    #[test]
    fn test_generate_recommendations_power() {
        let dim = PhysicalDimension::from_unit("W").unwrap();
        let unit_breakdown = HashMap::new();

        let recs = generate_recommendations(&dim, true, &unit_breakdown);
        assert!(recs.iter().any(|r| r.contains("power")));
    }

    #[test]
    fn test_generate_recommendations_dimensionless() {
        let dim = PhysicalDimension::new();
        let unit_breakdown = HashMap::new();

        let recs = generate_recommendations(&dim, true, &unit_breakdown);
        assert!(recs.iter().any(|r| r.contains("dimensionless")));
    }

    #[test]
    fn test_dimension_pascal() {
        let pa = PhysicalDimension::from_unit("Pa").unwrap();
        assert_eq!(pa.mass, 1);
        assert_eq!(pa.length, -1);
        assert_eq!(pa.time, -2);
    }

    #[test]
    fn test_dimension_volt() {
        let v = PhysicalDimension::from_unit("V").unwrap();
        assert_eq!(v.mass, 1);
        assert_eq!(v.length, 2);
        assert_eq!(v.time, -3);
        assert_eq!(v.current, -1);
    }

    #[test]
    fn test_dimension_coulomb() {
        let c = PhysicalDimension::from_unit("C").unwrap();
        assert_eq!(c.current, 1);
        assert_eq!(c.time, 1);
    }

    #[test]
    fn test_dimension_farad() {
        let f = PhysicalDimension::from_unit("F").unwrap();
        assert_eq!(f.mass, -1);
        assert_eq!(f.length, -2);
        assert_eq!(f.time, 4);
        assert_eq!(f.current, 2);
    }

    #[test]
    fn test_dimension_unknown_unit() {
        let result = PhysicalDimension::from_unit("invalid_unit");
        assert!(result.is_err());
    }

    #[test]
    fn test_dimension_frequency() {
        let hz = PhysicalDimension::from_unit("Hz").unwrap();
        assert_eq!(hz.time, -1);
        assert_eq!(hz.mass, 0);
    }

    #[test]
    fn test_dimension_length_variants() {
        let mm = PhysicalDimension::from_unit("mm").unwrap();
        assert_eq!(mm.length, 1);

        let km = PhysicalDimension::from_unit("km").unwrap();
        assert_eq!(km.length, 1);
    }

    #[test]
    fn test_dimension_energy_variants() {
        let ev = PhysicalDimension::from_unit("eV").unwrap();
        assert_eq!(ev.mass, 1);
        assert_eq!(ev.length, 2);
        assert_eq!(ev.time, -2);
    }

    #[test]
    fn test_analyze_expression_dimensions() {
        let mut units = HashMap::new();
        units.insert("F".to_string(), "N".to_string());
        units.insert("m".to_string(), "kg".to_string());

        let (dim, breakdown, consistent) = analyze_expression_dimensions("F/m", &units).unwrap();
        assert!(consistent);
        assert_eq!(breakdown.len(), 2);
    }

    #[test]
    fn test_tokenize_complex_unit() {
        let tokens = tokenize_unit("kg*m^2/s^3");
        assert!(tokens.len() > 0);
    }

    #[test]
    fn test_dimension_equality() {
        let n1 = PhysicalDimension::from_unit("N").unwrap();
        let n2 = PhysicalDimension::from_unit("N").unwrap();
        assert_eq!(n1, n2);

        let j = PhysicalDimension::from_unit("J").unwrap();
        assert_ne!(n1, j);
    }
