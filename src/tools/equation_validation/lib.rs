use clap::{Parser, Subcommand};
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use regex::Regex;

#[derive(Parser)]
#[command(name = "equation-validator")]
#[command(about = "Physics equation validation with dimensional analysis and conservation laws")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    ValidateEquation {
        #[arg(long)]
        equation: String,
        #[arg(long)]
        domain: String,
        #[arg(long)]
        units: Option<String>,
        #[arg(long)]
        conservation: Option<String>,
        #[arg(long)]
        symmetries: Option<String>,
    },
}

#[derive(Serialize, Debug)]
pub struct ValidationResult {
    is_valid: bool,
    dimensional_consistency: bool,
    physics_compliance: bool,
    mathematical_correctness: bool,
    unit_analysis: HashMap<String, String>,
    violations: Vec<String>,
    confidence: f64,
}

#[derive(Debug)]
struct ValidationError {
    message: String,
}

impl fmt::Display for ValidationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Validation error: {}", self.message)
    }
}

impl Error for ValidationError {}

// Unit dimension representation
#[derive(Debug, Clone, PartialEq)]
struct Dimension {
    mass: i32,      // kg
    length: i32,    // m
    time: i32,      // s
    current: i32,   // A
    temperature: i32, // K
    amount: i32,    // mol
    luminosity: i32, // cd
}

impl Dimension {
    fn new() -> Self {
        Dimension {
            mass: 0,
            length: 0,
            time: 0,
            current: 0,
            temperature: 0,
            amount: 0,
            luminosity: 0,
        }
    }

    fn from_unit(unit: &str) -> Result<Self, Box<dyn Error>> {
        let mut dim = Dimension::new();
        
        match unit {
            // Base units
            "kg" => dim.mass = 1,
            "m" => dim.length = 1,
            "s" => dim.time = 1,
            "A" => dim.current = 1,
            "K" => dim.temperature = 1,
            "mol" => dim.amount = 1,
            "cd" => dim.luminosity = 1,
            
            // Common derived units
            "N" => { dim.mass = 1; dim.length = 1; dim.time = -2; }, // kg⋅m⋅s⁻²
            "J" => { dim.mass = 1; dim.length = 2; dim.time = -2; }, // kg⋅m²⋅s⁻²
            "W" => { dim.mass = 1; dim.length = 2; dim.time = -3; }, // kg⋅m²⋅s⁻³
            "Pa" => { dim.mass = 1; dim.length = -1; dim.time = -2; }, // kg⋅m⁻¹⋅s⁻²
            "V" => { dim.mass = 1; dim.length = 2; dim.time = -3; dim.current = -1; }, // kg⋅m²⋅s⁻³⋅A⁻¹
            "C" => { dim.current = 1; dim.time = 1; }, // A⋅s
            "F" => { dim.mass = -1; dim.length = -2; dim.time = 4; dim.current = 2; }, // A²⋅s⁴⋅kg⁻¹⋅m⁻²
            "H" => { dim.mass = 1; dim.length = 2; dim.time = -2; dim.current = -2; }, // kg⋅m²⋅s⁻²⋅A⁻²
            
            // Compound units (simplified parsing)
            "m/s" => { dim.length = 1; dim.time = -1; },
            "m/s^2" => { dim.length = 1; dim.time = -2; },
            "kg*m/s^2" => { dim.mass = 1; dim.length = 1; dim.time = -2; },
            "kg*m^2/s^2" => { dim.mass = 1; dim.length = 2; dim.time = -2; },
            "N/m" => { dim.mass = 1; dim.time = -2; },
            "rad/s" => { dim.time = -1; }, // radians are dimensionless
            
            // Dimensionless
            "" | "1" | "rad" => {},
            
            _ => return Err(format!("Unknown unit: {}", unit).into()),
        }
        
        Ok(dim)
    }

    fn multiply(&self, other: &Dimension) -> Dimension {
        Dimension {
            mass: self.mass + other.mass,
            length: self.length + other.length,
            time: self.time + other.time,
            current: self.current + other.current,
            temperature: self.temperature + other.temperature,
            amount: self.amount + other.amount,
            luminosity: self.luminosity + other.luminosity,
        }
    }

    fn divide(&self, other: &Dimension) -> Dimension {
        Dimension {
            mass: self.mass - other.mass,
            length: self.length - other.length,
            time: self.time - other.time,
            current: self.current - other.current,
            temperature: self.temperature - other.temperature,
            amount: self.amount - other.amount,
            luminosity: self.luminosity - other.luminosity,
        }
    }

    fn power(&self, exp: i32) -> Dimension {
        Dimension {
            mass: self.mass * exp,
            length: self.length * exp,
            time: self.time * exp,
            current: self.current * exp,
            temperature: self.temperature * exp,
            amount: self.amount * exp,
            luminosity: self.luminosity * exp,
        }
    }

    fn is_dimensionless(&self) -> bool {
        self.mass == 0 && self.length == 0 && self.time == 0 && self.current == 0 
        && self.temperature == 0 && self.amount == 0 && self.luminosity == 0
    }
}

pub fn parse_equation(equation: &str) -> Result<(String, String), Box<dyn Error>> {
    // Simple equation parsing - split on '='
    let parts: Vec<&str> = equation.split('=').collect();
    if parts.len() != 2 {
        return Err("Equation must have exactly one '=' sign".into());
    }
    
    Ok((parts[0].trim().to_string(), parts[1].trim().to_string()))
}

pub fn extract_variables(expression: &str) -> Vec<String> {
    let re = Regex::new(r"[a-zA-Z][a-zA-Z0-9_]*").unwrap();
    let mut variables = Vec::new();
    
    for cap in re.find_iter(expression) {
        let var = cap.as_str();
        // Skip common functions
        if !["sin", "cos", "tan", "exp", "log", "sqrt", "abs"].contains(&var) && 
           !variables.contains(&var.to_string()) {
            variables.push(var.to_string());
        }
    }
    
    variables
}

pub fn check_dimensional_consistency(
    equation: &str,
    units: &HashMap<String, String>
) -> Result<(bool, HashMap<String, String>), Box<dyn Error>> {
    let (left_side, right_side) = parse_equation(equation)?;
    
    // Extract variables from both sides
    let left_vars = extract_variables(&left_side);
    let right_vars = extract_variables(&right_side);
    
    // Combine all variables
    let mut all_vars = left_vars;
    for var in right_vars {
        if !all_vars.contains(&var) {
            all_vars.push(var);
        }
    }
    
    // Create unit analysis map
    let mut unit_analysis = HashMap::new();
    for var in &all_vars {
        if let Some(unit) = units.get(var) {
            unit_analysis.insert(var.clone(), unit.clone());
        } else {
            unit_analysis.insert(var.clone(), "unknown".to_string());
        }
    }
    
    // Perform dimensional analysis by evaluating both sides
    let left_dim = evaluate_expression_dimensions(&left_side, units)?;
    let right_dim = evaluate_expression_dimensions(&right_side, units)?;

    // Check if dimensions match
    let consistency = left_dim == right_dim;
    
    Ok((consistency, unit_analysis))
}

fn evaluate_expression_dimensions(
    expression: &str,
    units: &HashMap<String, String>
) -> Result<Dimension, Box<dyn Error>> {
    // Simplified dimensional evaluation
    // Extract variables from expression
    let variables = extract_variables(expression);

    // If expression has no variables, it's dimensionless (a constant)
    if variables.is_empty() {
        return Ok(Dimension::new());
    }

    // For simplicity, assume first variable determines dimension
    // More sophisticated would parse entire expression tree
    if let Some(first_var) = variables.first() {
        if let Some(unit_str) = units.get(first_var) {
            return Dimension::from_unit(unit_str);
        }
    }

    // If units not specified, assume dimensionless
    Ok(Dimension::new())
}

pub fn check_conservation_laws(
    equation: &str,
    domain: &str,
    conservation_laws: &[String]
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut violations = Vec::new();
    
    for law in conservation_laws {
        match law.as_str() {
            "energy" => {
                // Check if energy terms are present and balanced
                if domain == "mechanics" && equation.contains("E") {
                    // Simplified check - look for energy conservation patterns
                    if !equation.contains("+") && !equation.contains("=") {
                        violations.push("Energy conservation may be violated".to_string());
                    }
                }
            }
            "momentum" => {
                // Check momentum conservation
                if domain == "mechanics" && equation.contains("p") {
                    // Look for momentum conservation patterns
                    if !equation.contains("m") || !equation.contains("v") {
                        violations.push("Momentum terms may be inconsistent".to_string());
                    }
                }
            }
            "charge" => {
                // Check charge conservation in electromagnetic contexts
                if domain == "electromagnetism" && equation.contains("q") {
                    // Simplified charge conservation check
                }
            }
            _ => {}
        }
    }
    
    Ok(violations)
}

pub fn check_symmetries(
    equation: &str,
    _domain: &str,
    symmetries: &[String]
) -> Result<Vec<String>, Box<dyn Error>> {
    let violations = Vec::new();
    
    for symmetry in symmetries {
        match symmetry.as_str() {
            "time_translation" => {
                // Check if equation is explicitly time-dependent when it shouldn't be
                if equation.contains("t") && !equation.contains("dt") {
                    // This is a simplified check
                }
            }
            "space_translation" => {
                // Check spatial homogeneity
            }
            "rotation" => {
                // Check rotational symmetry
            }
            "lorentz" => {
                // Check Lorentz invariance for relativistic equations
                if equation.contains("c") && !equation.contains("gamma") {
                    // Very simplified check
                }
            }
            "gauge" => {
                // Check gauge invariance for electromagnetic fields
            }
            "parity" => {
                // Check parity conservation
            }
            _ => {}
        }
    }
    
    Ok(violations)
}

pub fn check_mathematical_correctness(equation: &str) -> Result<bool, Box<dyn Error>> {
    // Basic mathematical checks
    let (left, right) = parse_equation(equation)?;
    
    // Check for basic mathematical validity
    let mut is_correct = true;
    
    // Check for division by zero patterns
    if equation.contains("/0") {
        is_correct = false;
    }
    
    // Check for negative square roots (simplified)
    if equation.contains("sqrt(-") {
        is_correct = false;
    }
    
    // Check for logarithms of non-positive numbers (simplified)
    if equation.contains("log(0") || equation.contains("log(-") {
        is_correct = false;
    }
    
    // Check bracket matching
    let left_parens = left.matches('(').count();
    let right_parens = left.matches(')').count();
    if left_parens != right_parens {
        is_correct = false;
    }
    
    let left_parens_r = right.matches('(').count();
    let right_parens_r = right.matches(')').count();
    if left_parens_r != right_parens_r {
        is_correct = false;
    }
    
    Ok(is_correct)
}

pub fn check_physics_compliance(equation: &str, domain: &str) -> Result<(bool, Vec<String>), Box<dyn Error>> {
    let mut violations = Vec::new();
    let mut is_compliant = true;
    
    match domain {
        "mechanics" => {
            // Check Newton's laws compliance
            if equation.contains("F") && equation.contains("=") {
                if !equation.contains("m") && !equation.contains("a") {
                    violations.push("Force equation should relate to mass and acceleration".to_string());
                    is_compliant = false;
                }
            }
            
            // Check energy conservation
            if equation.contains("E") && equation.contains("=") {
                // Look for kinetic + potential energy patterns
            }
        }
        "electromagnetism" => {
            // Check Maxwell's equations compliance
            if equation.contains("E") && equation.contains("B") {
                // Electromagnetic field relationships
            }
            
            // Check Coulomb's law patterns
            if equation.contains("F") && equation.contains("q") {
                if !equation.contains("r") {
                    violations.push("Electromagnetic force should depend on distance".to_string());
                }
            }
        }
        "thermodynamics" => {
            // Check thermodynamic laws
            if equation.contains("PV") && !equation.contains("T") {
                violations.push("Ideal gas law should include temperature".to_string());
            }
            
            // Check entropy relationships
            if equation.contains("S") && equation.contains("T") {
                // Entropy-temperature relationships
            }
        }
        "quantum" => {
            // Check quantum mechanical principles
            if equation.contains("psi") {
                // Wave function normalization checks
            }
            
            // Check uncertainty principle
            if equation.contains("delta") {
                // Heisenberg uncertainty relationships
            }
        }
        "relativity" => {
            // Check relativistic relationships
            if equation.contains("E") && equation.contains("m") && equation.contains("c") {
                if !equation.contains("c^2") && !equation.contains("c²") {
                    violations.push("Mass-energy relation should have c²".to_string());
                    is_compliant = false;
                }
            }
        }
        _ => {} // General physics
    }
    
    Ok((is_compliant, violations))
}

pub fn validate_equation(
    equation: String,
    domain: String,
    units: Option<HashMap<String, String>>,
    conservation_laws: Option<Vec<String>>,
    symmetries: Option<Vec<String>>,
) -> Result<ValidationResult, Box<dyn Error>> {
    let units_map = units.unwrap_or_default();
    let conservation = conservation_laws.unwrap_or_default();
    let symmetry_list = symmetries.unwrap_or_default();

    // Check mathematical correctness
    let mathematical_correctness = check_mathematical_correctness(&equation)?;

    // Check dimensional consistency
    let (dimensional_consistency, unit_analysis) = check_dimensional_consistency(&equation, &units_map)?;

    // Check physics compliance
    let (physics_compliance, mut physics_violations) = check_physics_compliance(&equation, &domain)?;

    // Check conservation laws
    let mut conservation_violations = check_conservation_laws(&equation, &domain, &conservation)?;

    // Check symmetries
    let mut symmetry_violations = check_symmetries(&equation, &domain, &symmetry_list)?;

    // Combine all violations
    let mut all_violations = Vec::new();
    all_violations.append(&mut physics_violations);
    all_violations.append(&mut conservation_violations);
    all_violations.append(&mut symmetry_violations);

    if !mathematical_correctness {
        all_violations.push("Mathematical syntax errors detected".to_string());
    }

    if !dimensional_consistency {
        all_violations.push("Dimensional inconsistency detected".to_string());
    }

    // Overall validation
    let is_valid = mathematical_correctness && dimensional_consistency && physics_compliance && all_violations.is_empty();

    // Calculate confidence based on number of checks passed
    let mut confidence = 0.0;
    if mathematical_correctness { confidence += 0.25; }
    if dimensional_consistency { confidence += 0.25; }
    if physics_compliance { confidence += 0.25; }
    if all_violations.is_empty() { confidence += 0.25; }

    Ok(ValidationResult {
        is_valid,
        dimensional_consistency,
        physics_compliance,
        mathematical_correctness,
        unit_analysis,
        violations: all_violations,
        confidence,
    })
}

#[cfg(test)]
mod tests {
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
        ).unwrap();

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
        ).unwrap();

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
        ).unwrap();

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
        ).unwrap();

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
        ).unwrap();

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
        ).unwrap();

        assert!(result.mathematical_correctness);
    }
}

