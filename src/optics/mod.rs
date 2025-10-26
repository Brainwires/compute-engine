/**
 * Optics Module
 *
 * Implements fundamental optics formulas:
 * - Thin lens equation
 * - Snell's law (refraction)
 * - Diffraction grating
 * - Fresnel equations (reflection/transmission)
 */
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpticsInput {
    pub operation: OpticsOperation,
    pub parameters: OpticsParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OpticsOperation {
    ThinLens,
    SnellsLaw,
    DiffractionGrating,
    FresnelEquations,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpticsParams {
    // Thin lens
    pub focal_length: Option<f64>,    // f (m)
    pub object_distance: Option<f64>, // do (m)
    pub image_distance: Option<f64>,  // di (m)

    // Snell's law
    pub n1: Option<f64>,     // refractive index 1
    pub n2: Option<f64>,     // refractive index 2
    pub theta1: Option<f64>, // incident angle (degrees)
    pub theta2: Option<f64>, // refracted angle (degrees)

    // Diffraction
    pub grating_spacing: Option<f64>, // d (m)
    pub wavelength: Option<f64>,      // λ (m)
    pub order: Option<i32>,           // m
    pub angle: Option<f64>,           // θ (degrees)

    // Fresnel
    pub polarization: Option<String>, // "s" or "p"
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpticsResult {
    pub primary_value: f64,
    pub unit: String,
    pub formula_used: String,
    pub secondary_values: Option<std::collections::HashMap<String, f64>>,
    pub interpretation: String,
}

pub fn calculate_optics(input: OpticsInput) -> Result<OpticsResult, String> {
    match input.operation {
        OpticsOperation::ThinLens => calculate_thin_lens(&input.parameters),
        OpticsOperation::SnellsLaw => calculate_snells_law(&input.parameters),
        OpticsOperation::DiffractionGrating => calculate_diffraction(&input.parameters),
        OpticsOperation::FresnelEquations => calculate_fresnel(&input.parameters),
    }
}

/// Thin lens equation: 1/f = 1/do + 1/di
fn calculate_thin_lens(params: &OpticsParams) -> Result<OpticsResult, String> {
    let f = params.focal_length;
    let d_o = params.object_distance;
    let d_i = params.image_distance;

    // Need at least 2 of the 3 values
    let count = [f.is_some(), d_o.is_some(), d_i.is_some()]
        .iter()
        .filter(|&&x| x)
        .count();
    if count < 2 {
        return Err(
            "Need at least 2 values: focal_length, object_distance, or image_distance".to_string(),
        );
    }

    let mut secondary = std::collections::HashMap::new();

    let (result, unit, interpretation) = if let Some(f_val) = f {
        if let Some(do_val) = d_o {
            // Calculate di: 1/di = 1/f - 1/do
            if do_val <= 0.0 {
                return Err("Object distance must be positive".to_string());
            }
            let di = 1.0 / (1.0 / f_val - 1.0 / do_val);
            let magnification = -di / do_val;

            secondary.insert("image_distance".to_string(), di);
            secondary.insert("magnification".to_string(), magnification);

            let img_type = if di > 0.0 { "real" } else { "virtual" };
            let orientation = if magnification > 0.0 {
                "upright"
            } else {
                "inverted"
            };

            (
                di.abs(),
                "m".to_string(),
                format!(
                    "{} {} image, {}, magnification: {:.2}x",
                    img_type,
                    orientation,
                    if magnification.abs() > 1.0 {
                        "enlarged"
                    } else {
                        "reduced"
                    },
                    magnification.abs()
                ),
            )
        } else {
            // d_i is provided, calculate d_o
            let di_val = d_i.unwrap();
            let do_calc = 1.0 / (1.0 / f_val - 1.0 / di_val);
            secondary.insert("object_distance".to_string(), do_calc);
            (
                do_calc.abs(),
                "m".to_string(),
                format!("Object distance: {:.3} m", do_calc),
            )
        }
    } else {
        // Both do and di provided, calculate f
        let do_val = d_o.unwrap();
        let di_val = d_i.unwrap();
        let f_calc = 1.0 / (1.0 / do_val + 1.0 / di_val);

        let lens_type = if f_calc > 0.0 {
            "converging"
        } else {
            "diverging"
        };
        secondary.insert("focal_length".to_string(), f_calc);

        (f_calc.abs(), "m".to_string(), format!("{} lens", lens_type))
    };

    Ok(OpticsResult {
        primary_value: result,
        unit,
        formula_used: "Thin lens: 1/f = 1/do + 1/di".to_string(),
        secondary_values: Some(secondary),
        interpretation,
    })
}

/// Snell's law: n1·sin(θ1) = n2·sin(θ2)
fn calculate_snells_law(params: &OpticsParams) -> Result<OpticsResult, String> {
    let n1 = params.n1.ok_or("Refractive index n1 required")?;
    let n2 = params.n2.ok_or("Refractive index n2 required")?;

    if n1 <= 0.0 || n2 <= 0.0 {
        return Err("Refractive indices must be positive".to_string());
    }

    let mut secondary = std::collections::HashMap::new();

    let (angle_result, is_incident) = if let Some(theta1_deg) = params.theta1 {
        // Calculate refracted angle
        let theta1 = theta1_deg.to_radians();
        let sin_theta2 = (n1 / n2) * theta1.sin();

        if sin_theta2.abs() > 1.0 {
            // Total internal reflection
            let critical_angle = ((n2 / n1).asin()).to_degrees();
            secondary.insert("critical_angle".to_string(), critical_angle);
            return Ok(OpticsResult {
                primary_value: critical_angle,
                unit: "degrees".to_string(),
                formula_used: "Snell's law: n1·sin(θ1) = n2·sin(θ2)".to_string(),
                secondary_values: Some(secondary),
                interpretation: format!(
                    "Total internal reflection! Critical angle: {:.2}°",
                    critical_angle
                ),
            });
        }

        let theta2 = sin_theta2.asin().to_degrees();
        secondary.insert("refracted_angle".to_string(), theta2);
        (theta2, false)
    } else if let Some(theta2_deg) = params.theta2 {
        // Calculate incident angle
        let theta2 = theta2_deg.to_radians();
        let sin_theta1 = (n2 / n1) * theta2.sin();
        let theta1 = sin_theta1.asin().to_degrees();
        secondary.insert("incident_angle".to_string(), theta1);
        (theta1, true)
    } else {
        return Err("Either theta1 or theta2 required".to_string());
    };

    // Calculate critical angle
    let critical = if n1 > n2 {
        Some((n2 / n1).asin().to_degrees())
    } else {
        None
    };

    if let Some(crit) = critical {
        secondary.insert("critical_angle".to_string(), crit);
    }

    let angle_type = if is_incident { "incident" } else { "refracted" };
    let bending = if n2 > n1 { "toward" } else { "away from" };

    Ok(OpticsResult {
        primary_value: angle_result,
        unit: "degrees".to_string(),
        formula_used: "Snell's law: n1·sin(θ1) = n2·sin(θ2)".to_string(),
        secondary_values: Some(secondary),
        interpretation: format!(
            "{} angle: {:.2}° (light bends {} normal)",
            angle_type, angle_result, bending
        ),
    })
}

/// Diffraction grating: d·sin(θ) = m·λ
fn calculate_diffraction(params: &OpticsParams) -> Result<OpticsResult, String> {
    let d = params.grating_spacing.ok_or("Grating spacing d required")?;
    let wavelength = params.wavelength.ok_or("Wavelength λ required")?;
    let m = params.order.unwrap_or(1);

    if d <= 0.0 || wavelength <= 0.0 {
        return Err("Spacing and wavelength must be positive".to_string());
    }

    let mut secondary = std::collections::HashMap::new();

    // Calculate angle: sin(θ) = m·λ/d
    let sin_theta = (m as f64) * wavelength / d;

    if sin_theta.abs() > 1.0 {
        return Err(format!("Order {} not observable (sin(θ) > 1)", m));
    }

    let theta = sin_theta.asin().to_degrees();

    // Calculate maximum order
    let m_max = (d / wavelength).floor() as i32;
    secondary.insert("max_order".to_string(), m_max as f64);
    secondary.insert("grating_lines_per_mm".to_string(), 1e-3 / d);

    Ok(OpticsResult {
        primary_value: theta,
        unit: "degrees".to_string(),
        formula_used: "Grating: d·sin(θ) = m·λ".to_string(),
        secondary_values: Some(secondary),
        interpretation: format!(
            "Order {} diffraction at {:.2}° (max order: {})",
            m, theta, m_max
        ),
    })
}

/// Fresnel equations for reflection/transmission coefficients
fn calculate_fresnel(params: &OpticsParams) -> Result<OpticsResult, String> {
    let n1 = params.n1.ok_or("Refractive index n1 required")?;
    let n2 = params.n2.ok_or("Refractive index n2 required")?;
    let theta1_deg = params.theta1.ok_or("Incident angle required")?;
    let polarization = params
        .polarization
        .as_ref()
        .map(|s| s.as_str())
        .unwrap_or("s");

    let theta1 = theta1_deg.to_radians();
    let sin_theta2 = (n1 / n2) * theta1.sin();

    if sin_theta2.abs() > 1.0 {
        return Err("Total internal reflection - no transmission".to_string());
    }

    let theta2 = sin_theta2.asin();

    let mut secondary = std::collections::HashMap::new();

    let (r, t) = match polarization {
        "s" => {
            // s-polarization (perpendicular to plane of incidence)
            let r_s =
                (n1 * theta1.cos() - n2 * theta2.cos()) / (n1 * theta1.cos() + n2 * theta2.cos());
            let t_s = (2.0 * n1 * theta1.cos()) / (n1 * theta1.cos() + n2 * theta2.cos());
            (r_s, t_s)
        }
        "p" => {
            // p-polarization (parallel to plane of incidence)
            let r_p =
                (n2 * theta1.cos() - n1 * theta2.cos()) / (n2 * theta1.cos() + n1 * theta2.cos());
            let t_p = (2.0 * n1 * theta1.cos()) / (n2 * theta1.cos() + n1 * theta2.cos());
            (r_p, t_p)
        }
        _ => return Err("Polarization must be 's' or 'p'".to_string()),
    };

    let reflectance = r * r; // Intensity reflection coefficient
    let transmittance = (n2 * theta2.cos()) / (n1 * theta1.cos()) * t * t;

    secondary.insert("reflection_coefficient".to_string(), reflectance);
    secondary.insert("transmission_coefficient".to_string(), transmittance);
    secondary.insert("refracted_angle".to_string(), theta2.to_degrees());

    // Calculate Brewster's angle for p-polarization
    if polarization == "p" {
        let brewster = (n2 / n1).atan().to_degrees();
        secondary.insert("brewsters_angle".to_string(), brewster);
    }

    Ok(OpticsResult {
        primary_value: reflectance * 100.0,
        unit: "% reflected".to_string(),
        formula_used: format!("Fresnel equations ({}-polarization)", polarization),
        secondary_values: Some(secondary),
        interpretation: format!(
            "{:.1}% reflected, {:.1}% transmitted",
            reflectance * 100.0,
            transmittance * 100.0
        ),
    })
}


impl Default for OpticsParams {
    fn default() -> Self {
        Self {
            focal_length: None,
            object_distance: None,
            image_distance: None,
            n1: None,
            n2: None,
            theta1: None,
            theta2: None,
            grating_spacing: None,
            wavelength: None,
            order: None,
            angle: None,
            polarization: None,
        }
    }
}
