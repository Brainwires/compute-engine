use clap::{Parser, Subcommand};
use regex::Regex;
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;

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
    pub expression: String,
    pub dimension: String,
    pub consistent: bool,
    pub unit_breakdown: HashMap<String, DimensionInfo>,
    pub analysis: String,
    pub recommendations: Vec<String>,
    pub target_match: Option<bool>,
}

#[derive(Serialize, Debug, Clone)]
pub struct DimensionInfo {
    pub unit: String,
    pub dimension: PhysicalDimension,
    pub power: i32,
}

/// Physical dimension representation using SI base units
#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct PhysicalDimension {
    /// Mass dimension (M)
    pub mass: i32,
    /// Length dimension (L)
    pub length: i32,
    /// Time dimension (T)
    pub time: i32,
    /// Electric current dimension (I)
    pub current: i32,
    /// Temperature dimension (Θ)
    pub temperature: i32,
    /// Amount of substance dimension (N)
    pub amount: i32,
    /// Luminous intensity dimension (J)
    pub luminosity: i32,
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
        let _dim = PhysicalDimension::new();

        // Parse compound units like "kg*m/s^2" or "m/s^2"
        let normalized = unit
            .replace("*", " * ")
            .replace("/", " / ")
            .replace("^", " ^ ");
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
            "N" => {
                dim.mass = 1;
                dim.length = 1;
                dim.time = -2;
            } // Newton
            "J" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -2;
            } // Joule
            "W" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -3;
            } // Watt
            "Pa" => {
                dim.mass = 1;
                dim.length = -1;
                dim.time = -2;
            } // Pascal
            "V" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -3;
                dim.current = -1;
            } // Volt
            "C" => {
                dim.current = 1;
                dim.time = 1;
            } // Coulomb
            "F" => {
                dim.mass = -1;
                dim.length = -2;
                dim.time = 4;
                dim.current = 2;
            } // Farad
            "H" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -2;
                dim.current = -2;
            } // Henry
            "Ω" | "ohm" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -3;
                dim.current = -2;
            } // Ohm
            "T" => {
                dim.mass = 1;
                dim.time = -2;
                dim.current = -1;
            } // Tesla
            "Wb" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -2;
                dim.current = -1;
            } // Weber

            // Special units
            "rad" | "sr" => {}     // Dimensionless
            "Hz" => dim.time = -1, // Hertz
            "g" => dim.mass = 1,   // gram (treat same as kg for dimensions)

            // Length units
            "mm" | "cm" | "km" => dim.length = 1,

            // Time units
            "min" | "h" | "hr" | "day" | "year" => dim.time = 1,

            // Mass units
            "ton" | "tonne" => dim.mass = 1,

            // Energy units
            "eV" | "keV" | "MeV" | "GeV" | "cal" | "kcal" | "kWh" => {
                dim.mass = 1;
                dim.length = 2;
                dim.time = -2;
            }

            // Force units
            "dyn" | "lbf" => {
                dim.mass = 1;
                dim.length = 1;
                dim.time = -2;
            }

            // Pressure units
            "bar" | "atm" | "mmHg" | "Torr" | "psi" => {
                dim.mass = 1;
                dim.length = -1;
                dim.time = -2;
            }

            // Unknown or dimensionless
            "" | "1" => {}

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
        self.mass == 0
            && self.length == 0
            && self.time == 0
            && self.current == 0
            && self.temperature == 0
            && self.amount == 0
            && self.luminosity == 0
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
        if [
            "sin", "cos", "tan", "exp", "log", "ln", "sqrt", "abs", "min", "max",
        ]
        .contains(&var_name)
        {
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

                    unit_breakdown.insert(
                        var_name.clone(),
                        DimensionInfo {
                            unit: unit.clone(),
                            dimension: var_dimension,
                            power,
                        },
                    );
                }
                Err(_) => {
                    is_consistent = false;
                    unit_breakdown.insert(
                        var_name.clone(),
                        DimensionInfo {
                            unit: unit.clone(),
                            dimension: PhysicalDimension::new(),
                            power,
                        },
                    );
                }
            }
        } else {
            is_consistent = false;
            unit_breakdown.insert(
                var_name.clone(),
                DimensionInfo {
                    unit: "unknown".to_string(),
                    dimension: PhysicalDimension::new(),
                    power,
                },
            );
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
        recommendations.push(
            "Some variables have unknown units - please specify units for all variables"
                .to_string(),
        );
    }

    // Check for common dimensional patterns
    match (dimension.mass, dimension.length, dimension.time) {
        (1, 1, -2) => {
            recommendations.push("This expression has dimensions of force [MLT⁻²]".to_string())
        }
        (1, 2, -2) => {
            recommendations.push("This expression has dimensions of energy [ML²T⁻²]".to_string())
        }
        (1, 2, -3) => {
            recommendations.push("This expression has dimensions of power [ML²T⁻³]".to_string())
        }
        (1, -1, -2) => {
            recommendations.push("This expression has dimensions of pressure [ML⁻¹T⁻²]".to_string())
        }
        (0, 1, -1) => {
            recommendations.push("This expression has dimensions of velocity [LT⁻¹]".to_string())
        }
        (0, 1, -2) => recommendations
            .push("This expression has dimensions of acceleration [LT⁻²]".to_string()),
        (0, 0, -1) => {
            recommendations.push("This expression has dimensions of frequency [T⁻¹]".to_string())
        }
        _ => {}
    }

    if dimension.is_dimensionless() {
        recommendations.push(
            "This expression is dimensionless - good for ratios, angles, or pure numbers"
                .to_string(),
        );
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
        format!(
            "Expression '{}' has dimensions [{}]",
            expression, dimension_string
        )
    } else {
        format!(
            "Expression '{}' has incomplete dimensional information",
            expression
        )
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

/// Unit conversion structure
#[derive(Serialize, Debug)]
pub struct UnitConversion {
    pub from_unit: String,
    pub to_unit: String,
    pub factor: f64,
    pub dimension: String,
}

/// Convert between units
///
/// # Examples
/// ```
/// use computational_engine::tools::dimensional_analysis::convert_units;
/// let result = convert_units("m", "km", 1000.0).unwrap();
/// assert_eq!(result.factor, 0.001);
/// ```
pub fn convert_units(
    from_unit: &str,
    to_unit: &str,
    value: f64,
) -> Result<UnitConversion, Box<dyn Error>> {
    let from_dim = PhysicalDimension::from_unit(from_unit)?;
    let to_dim = PhysicalDimension::from_unit(to_unit)?;

    if from_dim != to_dim {
        return Err(format!(
            "Cannot convert between incompatible dimensions: {} ({}) and {} ({})",
            from_unit,
            from_dim.to_string(),
            to_unit,
            to_dim.to_string()
        )
        .into());
    }

    // Get conversion factor
    let factor = get_conversion_factor(from_unit, to_unit)?;
    let converted_value = value * factor;

    Ok(UnitConversion {
        from_unit: from_unit.to_string(),
        to_unit: to_unit.to_string(),
        factor: converted_value,
        dimension: from_dim.to_string(),
    })
}

/// Get conversion factor between compatible units
fn get_conversion_factor(from_unit: &str, to_unit: &str) -> Result<f64, Box<dyn Error>> {
    // Define conversion factors to SI base units
    let to_si = |unit: &str| -> Result<f64, Box<dyn Error>> {
        Ok(match unit {
            // Length
            "m" => 1.0,
            "km" => 1000.0,
            "cm" => 0.01,
            "mm" => 0.001,
            "ft" | "feet" => 0.3048,
            "in" | "inch" => 0.0254,
            "mi" | "mile" => 1609.34,

            // Mass
            "kg" => 1.0,
            "g" => 0.001,
            "mg" => 1e-6,
            "ton" | "tonne" => 1000.0,
            "lb" | "pound" => 0.453592,
            "oz" | "ounce" => 0.0283495,

            // Time
            "s" => 1.0,
            "min" => 60.0,
            "h" | "hr" | "hour" => 3600.0,
            "day" => 86400.0,
            "week" => 604800.0,
            "year" => 31536000.0,

            // Force
            "N" => 1.0,
            "kN" => 1000.0,
            "dyn" => 1e-5,
            "lbf" => 4.44822,

            // Energy
            "J" => 1.0,
            "kJ" => 1000.0,
            "MJ" => 1e6,
            "eV" => 1.602176634e-19,
            "keV" => 1.602176634e-16,
            "MeV" => 1.602176634e-13,
            "GeV" => 1.602176634e-10,
            "cal" => 4.184,
            "kcal" => 4184.0,
            "kWh" => 3.6e6,
            "BTU" => 1055.06,

            // Power
            "W" => 1.0,
            "kW" => 1000.0,
            "MW" => 1e6,
            "hp" => 745.7,

            // Pressure
            "Pa" => 1.0,
            "kPa" => 1000.0,
            "MPa" => 1e6,
            "bar" => 1e5,
            "atm" => 101325.0,
            "mmHg" | "Torr" => 133.322,
            "psi" => 6894.76,

            // Temperature (handled specially)
            "K" => 1.0,
            "C" | "°C" => 1.0, // Delta T
            "F" | "°F" => 5.0 / 9.0, // Delta T

            // Area
            "m^2" | "m²" => 1.0,
            "km^2" | "km²" => 1e6,
            "cm^2" | "cm²" => 1e-4,
            "mm^2" | "mm²" => 1e-6,

            // Volume
            "m^3" | "m³" => 1.0,
            "L" | "liter" => 0.001,
            "mL" | "milliliter" => 1e-6,
            "gal" | "gallon" => 0.00378541,

            // Frequency
            "Hz" => 1.0,
            "kHz" => 1000.0,
            "MHz" => 1e6,
            "GHz" => 1e9,

            // Velocity
            "m/s" => 1.0,
            "km/h" | "kph" => 1.0 / 3.6,
            "mph" => 0.44704,

            _ => return Err(format!("Unknown unit for conversion: {}", unit).into()),
        })
    };

    let from_factor = to_si(from_unit)?;
    let to_factor = to_si(to_unit)?;

    Ok(from_factor / to_factor)
}

/// Derive units for common physics quantities
///
/// # Examples
/// ```
/// use computational_engine::tools::dimensional_analysis::derive_units_for_quantity;
/// let units = derive_units_for_quantity("force").unwrap();
/// assert!(units.contains(&"N".to_string()));
/// ```
pub fn derive_units_for_quantity(quantity: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let units = match quantity.to_lowercase().as_str() {
        "force" => vec!["N", "kN", "dyn", "lbf"],
        "energy" | "work" => vec!["J", "kJ", "MJ", "eV", "keV", "MeV", "cal", "kcal", "kWh"],
        "power" => vec!["W", "kW", "MW", "hp"],
        "pressure" | "stress" => vec!["Pa", "kPa", "MPa", "bar", "atm", "mmHg", "psi"],
        "length" | "distance" => vec!["m", "km", "cm", "mm", "ft", "in", "mi"],
        "mass" => vec!["kg", "g", "mg", "ton", "lb", "oz"],
        "time" => vec!["s", "min", "h", "day", "year"],
        "velocity" | "speed" => vec!["m/s", "km/h", "mph"],
        "acceleration" => vec!["m/s²", "g"],
        "frequency" => vec!["Hz", "kHz", "MHz", "GHz"],
        "electric_current" | "current" => vec!["A", "mA", "µA"],
        "voltage" | "potential" => vec!["V", "kV", "mV"],
        "resistance" => vec!["Ω", "ohm", "kΩ", "MΩ"],
        "capacitance" => vec!["F", "µF", "nF", "pF"],
        "inductance" => vec!["H", "mH", "µH"],
        "magnetic_field" => vec!["T", "mT", "G"],
        "electric_charge" => vec!["C", "mC", "µC"],
        "temperature" => vec!["K", "°C", "°F"],
        "area" => vec!["m²", "km²", "cm²", "mm²"],
        "volume" => vec!["m³", "L", "mL", "gal"],
        _ => {
            return Err(format!("Unknown physical quantity: {}", quantity).into());
        }
    };

    Ok(units.iter().map(|s| s.to_string()).collect())
}

/// Check if units are dimensionally compatible
///
/// # Examples
/// ```
/// use computational_engine::tools::dimensional_analysis::are_units_compatible;
/// assert!(are_units_compatible("m", "km").unwrap());
/// assert!(!are_units_compatible("m", "kg").unwrap());
/// ```
pub fn are_units_compatible(unit1: &str, unit2: &str) -> Result<bool, Box<dyn Error>> {
    let dim1 = PhysicalDimension::from_unit(unit1)?;
    let dim2 = PhysicalDimension::from_unit(unit2)?;
    Ok(dim1 == dim2)
}

/// Get the SI base units for a given unit
pub fn get_si_base_units(unit: &str) -> Result<String, Box<dyn Error>> {
    let dim = PhysicalDimension::from_unit(unit)?;
    Ok(dim.to_string())
}

/// List common physics formulas with their dimensional analysis
#[derive(Serialize, Debug)]
pub struct PhysicsFormula {
    pub name: String,
    pub formula: String,
    pub dimension: String,
    pub si_unit: String,
    pub description: String,
}

pub fn get_physics_formulas() -> Vec<PhysicsFormula> {
    vec![
        PhysicsFormula {
            name: "Newton's Second Law".to_string(),
            formula: "F = ma".to_string(),
            dimension: "M⋅L⋅T^-2".to_string(),
            si_unit: "N (Newton)".to_string(),
            description: "Force equals mass times acceleration".to_string(),
        },
        PhysicsFormula {
            name: "Kinetic Energy".to_string(),
            formula: "KE = (1/2)mv²".to_string(),
            dimension: "M⋅L^2⋅T^-2".to_string(),
            si_unit: "J (Joule)".to_string(),
            description: "Energy of motion".to_string(),
        },
        PhysicsFormula {
            name: "Gravitational Potential Energy".to_string(),
            formula: "PE = mgh".to_string(),
            dimension: "M⋅L^2⋅T^-2".to_string(),
            si_unit: "J (Joule)".to_string(),
            description: "Energy due to position in gravitational field".to_string(),
        },
        PhysicsFormula {
            name: "Power".to_string(),
            formula: "P = W/t".to_string(),
            dimension: "M⋅L^2⋅T^-3".to_string(),
            si_unit: "W (Watt)".to_string(),
            description: "Rate of energy transfer".to_string(),
        },
        PhysicsFormula {
            name: "Pressure".to_string(),
            formula: "P = F/A".to_string(),
            dimension: "M⋅L^-1⋅T^-2".to_string(),
            si_unit: "Pa (Pascal)".to_string(),
            description: "Force per unit area".to_string(),
        },
        PhysicsFormula {
            name: "Momentum".to_string(),
            formula: "p = mv".to_string(),
            dimension: "M⋅L⋅T^-1".to_string(),
            si_unit: "kg⋅m/s".to_string(),
            description: "Mass times velocity".to_string(),
        },
        PhysicsFormula {
            name: "Impulse".to_string(),
            formula: "J = FΔt".to_string(),
            dimension: "M⋅L⋅T^-1".to_string(),
            si_unit: "N⋅s".to_string(),
            description: "Change in momentum".to_string(),
        },
        PhysicsFormula {
            name: "Angular Momentum".to_string(),
            formula: "L = Iω".to_string(),
            dimension: "M⋅L^2⋅T^-1".to_string(),
            si_unit: "kg⋅m²/s".to_string(),
            description: "Rotational momentum".to_string(),
        },
        PhysicsFormula {
            name: "Electric Field".to_string(),
            formula: "E = F/q".to_string(),
            dimension: "M⋅L⋅T^-3⋅I^-1".to_string(),
            si_unit: "V/m".to_string(),
            description: "Force per unit charge".to_string(),
        },
        PhysicsFormula {
            name: "Magnetic Field".to_string(),
            formula: "B = F/(qv)".to_string(),
            dimension: "M⋅T^-2⋅I^-1".to_string(),
            si_unit: "T (Tesla)".to_string(),
            description: "Magnetic flux density".to_string(),
        },
    ]
}

