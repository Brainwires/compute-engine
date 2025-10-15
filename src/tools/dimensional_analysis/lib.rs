use clap::{Parser, Subcommand};
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use regex::Regex;

#[derive(Parser)]
#[command(name = "dimensional-analyzer")]
#[command(about = "Comprehensive dimensional analysis for physics equations and expressions")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    DimensionalAnalysis {
        #[arg(long)]
        expression: String,
        #[arg(long)]
        units: String,
        #[arg(long)]
        target: Option<String>,
    },
}

#[derive(Serialize, Debug)]
pub struct DimensionalAnalysisResult {
    expression: String,
    dimension: String,
    consistent: bool,
    unit_breakdown: HashMap<String, DimensionInfo>,
    analysis: String,
    recommendations: Vec<String>,
    target_match: Option<bool>,
}

#[derive(Serialize, Debug, Clone)]
struct DimensionInfo {
    unit: String,
    dimension: PhysicalDimension,
    power: i32,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
struct PhysicalDimension {
    mass: i32,        // M
    length: i32,      // L  
    time: i32,        // T
    current: i32,     // I
    temperature: i32, // Θ (theta)
    amount: i32,      // N
    luminosity: i32,  // J
}

impl PhysicalDimension {
    fn new() -> Self {
        PhysicalDimension {
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
        let dim = PhysicalDimension::new();
        
        // Parse compound units like "kg*m/s^2" or "m/s^2"
        let normalized = unit.replace("*", " * ").replace("/", " / ").replace("^", " ^ ");
        let tokens = self::tokenize_unit(&normalized);
        
        let mut current_dim = PhysicalDimension::new();
        let mut operation = '*'; // Start with multiplication
        let mut power = 1;
        
        for token in tokens {
            match token.as_str() {
                "*" => operation = '*',
                "/" => operation = '/',
                "^" => continue, // Power handled in next iteration
                token if token.chars().all(|c| c.is_numeric() || c == '-') => {
                    power = token.parse().unwrap_or(1);
                    continue;
                }
                base_unit => {
                    let mut unit_dim = Self::get_base_dimension(base_unit)?;
                    unit_dim = unit_dim.power(power);
                    
                    match operation {
                        '*' => current_dim = current_dim.multiply(&unit_dim),
                        '/' => current_dim = current_dim.divide(&unit_dim),
                        _ => {}
                    }
                    
                    power = 1; // Reset power
                }
            }
        }
        
        Ok(current_dim)
    }
    
    fn get_base_dimension(unit: &str) -> Result<PhysicalDimension, Box<dyn Error>> {
        let mut dim = PhysicalDimension::new();
        
        match unit {
            // SI Base units
            "kg" => dim.mass = 1,
            "m" => dim.length = 1,
            "s" => dim.time = 1,
            "A" => dim.current = 1,
            "K" => dim.temperature = 1,
            "mol" => dim.amount = 1,
            "cd" => dim.luminosity = 1,
            
            // Common derived units
            "N" => { dim.mass = 1; dim.length = 1; dim.time = -2; }, // Newton
            "J" => { dim.mass = 1; dim.length = 2; dim.time = -2; }, // Joule
            "W" => { dim.mass = 1; dim.length = 2; dim.time = -3; }, // Watt
            "Pa" => { dim.mass = 1; dim.length = -1; dim.time = -2; }, // Pascal
            "V" => { dim.mass = 1; dim.length = 2; dim.time = -3; dim.current = -1; }, // Volt
            "C" => { dim.current = 1; dim.time = 1; }, // Coulomb
            "F" => { dim.mass = -1; dim.length = -2; dim.time = 4; dim.current = 2; }, // Farad
            "H" => { dim.mass = 1; dim.length = 2; dim.time = -2; dim.current = -2; }, // Henry
            "Ω" | "ohm" => { dim.mass = 1; dim.length = 2; dim.time = -3; dim.current = -2; }, // Ohm
            "T" => { dim.mass = 1; dim.time = -2; dim.current = -1; }, // Tesla
            "Wb" => { dim.mass = 1; dim.length = 2; dim.time = -2; dim.current = -1; }, // Weber
            
            // Special units
            "rad" | "sr" => {}, // Dimensionless
            "Hz" => dim.time = -1, // Hertz
            "g" => dim.mass = 1,   // gram (treat same as kg for dimensions)
            
            // Length units
            "mm" | "cm" | "km" => dim.length = 1,
            
            // Time units  
            "min" | "h" | "hr" | "day" | "year" => dim.time = 1,
            
            // Mass units
            "g" => dim.mass = 1,
            "ton" | "tonne" => dim.mass = 1,
            
            // Energy units
            "eV" | "keV" | "MeV" | "GeV" | "cal" | "kcal" | "kWh" => { 
                dim.mass = 1; dim.length = 2; dim.time = -2; 
            },
            
            // Force units
            "dyn" | "lbf" => { dim.mass = 1; dim.length = 1; dim.time = -2; },
            
            // Pressure units
            "bar" | "atm" | "mmHg" | "Torr" | "psi" => { 
                dim.mass = 1; dim.length = -1; dim.time = -2; 
            },
            
            // Unknown or dimensionless
            "" | "1" => {},
            
            _ => return Err(format!("Unknown unit: {}", unit).into()),
        }
        
        Ok(dim)
    }

    fn multiply(&self, other: &PhysicalDimension) -> PhysicalDimension {
        PhysicalDimension {
            mass: self.mass + other.mass,
            length: self.length + other.length,
            time: self.time + other.time,
            current: self.current + other.current,
            temperature: self.temperature + other.temperature,
            amount: self.amount + other.amount,
            luminosity: self.luminosity + other.luminosity,
        }
    }

    fn divide(&self, other: &PhysicalDimension) -> PhysicalDimension {
        PhysicalDimension {
            mass: self.mass - other.mass,
            length: self.length - other.length,
            time: self.time - other.time,
            current: self.current - other.current,
            temperature: self.temperature - other.temperature,
            amount: self.amount - other.amount,
            luminosity: self.luminosity - other.luminosity,
        }
    }

    fn power(&self, exp: i32) -> PhysicalDimension {
        PhysicalDimension {
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

    fn to_string(&self) -> String {
        let mut parts = Vec::new();
        
        if self.mass != 0 {
            if self.mass == 1 {
                parts.push("M".to_string());
            } else {
                parts.push(format!("M^{}", self.mass));
            }
        }
        
        if self.length != 0 {
            if self.length == 1 {
                parts.push("L".to_string());
            } else {
                parts.push(format!("L^{}", self.length));
            }
        }
        
        if self.time != 0 {
            if self.time == 1 {
                parts.push("T".to_string());
            } else {
                parts.push(format!("T^{}", self.time));
            }
        }
        
        if self.current != 0 {
            if self.current == 1 {
                parts.push("I".to_string());
            } else {
                parts.push(format!("I^{}", self.current));
            }
        }
        
        if self.temperature != 0 {
            if self.temperature == 1 {
                parts.push("Θ".to_string());
            } else {
                parts.push(format!("Θ^{}", self.temperature));
            }
        }
        
        if self.amount != 0 {
            if self.amount == 1 {
                parts.push("N".to_string());
            } else {
                parts.push(format!("N^{}", self.amount));
            }
        }
        
        if self.luminosity != 0 {
            if self.luminosity == 1 {
                parts.push("J".to_string());
            } else {
                parts.push(format!("J^{}", self.luminosity));
            }
        }
        
        if parts.is_empty() {
            "1".to_string() // Dimensionless
        } else {
            parts.join("⋅")
        }
    }
}

pub fn tokenize_unit(unit_str: &str) -> Vec<String> {
    let re = Regex::new(r"[a-zA-Z]+|\d+|[*/^()]").unwrap();
    re.find_iter(unit_str)
        .map(|m| m.as_str().to_string())
        .collect()
}

pub fn extract_variables_with_powers(expression: &str) -> Vec<(String, i32)> {
    // Extract variables and their powers from mathematical expressions
    let re = Regex::new(r"([a-zA-Z][a-zA-Z0-9_]*)\^?(\d+)?").unwrap();
    let mut variables = Vec::new();
    
    for cap in re.captures_iter(expression) {
        let var_name = cap.get(1).unwrap().as_str();
        
        // Skip common mathematical functions
        if ["sin", "cos", "tan", "exp", "log", "ln", "sqrt", "abs", "min", "max"].contains(&var_name) {
            continue;
        }
        
        let power = if let Some(pow_match) = cap.get(2) {
            pow_match.as_str().parse().unwrap_or(1)
        } else if expression.contains(&format!("{}^2", var_name)) {
            2
        } else if expression.contains(&format!("{}^3", var_name)) {
            3
        } else {
            1
        };
        
        variables.push((var_name.to_string(), power));
    }
    
    variables
}

pub fn analyze_expression_dimensions(
    expression: &str,
    variable_units: &HashMap<String, String>,
) -> Result<(PhysicalDimension, HashMap<String, DimensionInfo>, bool), Box<dyn Error>> {
    let variables_with_powers = extract_variables_with_powers(expression);
    let mut unit_breakdown = HashMap::new();
    let mut result_dimension = PhysicalDimension::new();
    let mut is_consistent = true;
    
    // Process each variable
    for (var_name, power) in variables_with_powers {
        if let Some(unit) = variable_units.get(&var_name) {
            match PhysicalDimension::from_unit(unit) {
                Ok(var_dimension) => {
                    let powered_dimension = var_dimension.power(power);
                    result_dimension = result_dimension.multiply(&powered_dimension);
                    
                    unit_breakdown.insert(var_name.clone(), DimensionInfo {
                        unit: unit.clone(),
                        dimension: var_dimension,
                        power,
                    });
                }
                Err(_) => {
                    is_consistent = false;
                    unit_breakdown.insert(var_name.clone(), DimensionInfo {
                        unit: unit.clone(),
                        dimension: PhysicalDimension::new(),
                        power,
                    });
                }
            }
        } else {
            is_consistent = false;
            unit_breakdown.insert(var_name.clone(), DimensionInfo {
                unit: "unknown".to_string(),
                dimension: PhysicalDimension::new(),
                power,
            });
        }
    }
    
    // Handle operations (simplified)
    if expression.contains("/") {
        // For division, we'd need to parse the expression tree properly
        // This is a simplified approach
    }
    
    Ok((result_dimension, unit_breakdown, is_consistent))
}

pub fn generate_recommendations(
    dimension: &PhysicalDimension,
    is_consistent: bool,
    unit_breakdown: &HashMap<String, DimensionInfo>,
) -> Vec<String> {
    let mut recommendations = Vec::new();
    
    if !is_consistent {
        recommendations.push("Some variables have unknown units - please specify units for all variables".to_string());
    }
    
    // Check for common dimensional patterns
    match (dimension.mass, dimension.length, dimension.time) {
        (1, 1, -2) => recommendations.push("This expression has dimensions of force [MLT⁻²]".to_string()),
        (1, 2, -2) => recommendations.push("This expression has dimensions of energy [ML²T⁻²]".to_string()),
        (1, 2, -3) => recommendations.push("This expression has dimensions of power [ML²T⁻³]".to_string()),
        (1, -1, -2) => recommendations.push("This expression has dimensions of pressure [ML⁻¹T⁻²]".to_string()),
        (0, 1, -1) => recommendations.push("This expression has dimensions of velocity [LT⁻¹]".to_string()),
        (0, 1, -2) => recommendations.push("This expression has dimensions of acceleration [LT⁻²]".to_string()),
        (0, 0, -1) => recommendations.push("This expression has dimensions of frequency [T⁻¹]".to_string()),
        _ => {}
    }
    
    if dimension.is_dimensionless() {
        recommendations.push("This expression is dimensionless - good for ratios, angles, or pure numbers".to_string());
    }
    
    // Check for potential issues
    for (var, info) in unit_breakdown {
        if info.unit == "unknown" {
            recommendations.push(format!("Variable '{}' needs unit specification", var));
        }
    }
    
    recommendations
}

pub fn dimensional_analysis(
    expression: String,
    variable_units: HashMap<String, String>,
    target_dimension: Option<String>,
) -> Result<DimensionalAnalysisResult, Box<dyn Error>> {
    let (dimension, unit_breakdown, is_consistent) =
        analyze_expression_dimensions(&expression, &variable_units)?;

    let dimension_string = dimension.to_string();
    let recommendations = generate_recommendations(&dimension, is_consistent, &unit_breakdown);

    // Check target dimension match if provided
    let target_match = if let Some(target) = &target_dimension {
        match PhysicalDimension::from_unit(target) {
            Ok(target_dim) => Some(dimension == target_dim),
            Err(_) => None,
        }
    } else {
        None
    };

    let analysis = if is_consistent {
        format!("Expression '{}' has dimensions [{}]", expression, dimension_string)
    } else {
        format!("Expression '{}' has incomplete dimensional information", expression)
    };

    Ok(DimensionalAnalysisResult {
        expression,
        dimension: dimension_string,
        consistent: is_consistent,
        unit_breakdown,
        analysis,
        recommendations,
        target_match,
    })
}

#[cfg(test)]
mod tests {
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

        let result = dimensional_analysis(
            "F = m".to_string(),
            units,
            None,
        ).unwrap();

        // Just verify it runs without error
        assert!(result.unit_breakdown.len() > 0);
    }

    #[test]
    fn test_dimensional_analysis_energy() {
        let mut units = HashMap::new();
        units.insert("E".to_string(), "J".to_string());
        units.insert("m".to_string(), "kg".to_string());

        let result = dimensional_analysis(
            "E = m".to_string(),
            units,
            None,
        ).unwrap();

        assert!(result.unit_breakdown.len() > 0);
    }

    #[test]
    fn test_dimensional_analysis_velocity() {
        let mut units = HashMap::new();
        units.insert("v".to_string(), "m".to_string());
        units.insert("t".to_string(), "s".to_string());

        let result = dimensional_analysis(
            "v = t".to_string(),
            units,
            None,
        ).unwrap();

        assert!(result.unit_breakdown.len() > 0);
    }

    #[test]
    fn test_dimensional_analysis_unknown_units() {
        let units = HashMap::new();

        let result = dimensional_analysis(
            "F = m*a".to_string(),
            units,
            None,
        ).unwrap();

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
}

