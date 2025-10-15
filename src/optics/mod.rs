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
    pub focal_length: Option<f64>,      // f (m)
    pub object_distance: Option<f64>,    // do (m)
    pub image_distance: Option<f64>,     // di (m)

    // Snell's law
    pub n1: Option<f64>,                 // refractive index 1
    pub n2: Option<f64>,                 // refractive index 2
    pub theta1: Option<f64>,             // incident angle (degrees)
    pub theta2: Option<f64>,             // refracted angle (degrees)

    // Diffraction
    pub grating_spacing: Option<f64>,    // d (m)
    pub wavelength: Option<f64>,         // λ (m)
    pub order: Option<i32>,              // m
    pub angle: Option<f64>,              // θ (degrees)

    // Fresnel
    pub polarization: Option<String>,    // "s" or "p"
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
    let count = [f.is_some(), d_o.is_some(), d_i.is_some()].iter().filter(|&&x| x).count();
    if count < 2 {
        return Err("Need at least 2 values: focal_length, object_distance, or image_distance".to_string());
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
            let orientation = if magnification > 0.0 { "upright" } else { "inverted" };

            (di.abs(), "m".to_string(),
             format!("{} {} image, {}, magnification: {:.2}x",
                     img_type, orientation, if magnification.abs() > 1.0 { "enlarged" } else { "reduced" }, magnification.abs()))
        } else {
            // d_i is provided, calculate d_o
            let di_val = d_i.unwrap();
            let do_calc = 1.0 / (1.0 / f_val - 1.0 / di_val);
            secondary.insert("object_distance".to_string(), do_calc);
            (do_calc.abs(), "m".to_string(), format!("Object distance: {:.3} m", do_calc))
        }
    } else {
        // Both do and di provided, calculate f
        let do_val = d_o.unwrap();
        let di_val = d_i.unwrap();
        let f_calc = 1.0 / (1.0 / do_val + 1.0 / di_val);

        let lens_type = if f_calc > 0.0 { "converging" } else { "diverging" };
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
                interpretation: format!("Total internal reflection! Critical angle: {:.2}°", critical_angle),
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
        interpretation: format!("{} angle: {:.2}° (light bends {} normal)", angle_type, angle_result, bending),
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
        interpretation: format!("Order {} diffraction at {:.2}° (max order: {})", m, theta, m_max),
    })
}

/// Fresnel equations for reflection/transmission coefficients
fn calculate_fresnel(params: &OpticsParams) -> Result<OpticsResult, String> {
    let n1 = params.n1.ok_or("Refractive index n1 required")?;
    let n2 = params.n2.ok_or("Refractive index n2 required")?;
    let theta1_deg = params.theta1.ok_or("Incident angle required")?;
    let polarization = params.polarization.as_ref()
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
            let r_s = (n1 * theta1.cos() - n2 * theta2.cos()) /
                      (n1 * theta1.cos() + n2 * theta2.cos());
            let t_s = (2.0 * n1 * theta1.cos()) /
                      (n1 * theta1.cos() + n2 * theta2.cos());
            (r_s, t_s)
        }
        "p" => {
            // p-polarization (parallel to plane of incidence)
            let r_p = (n2 * theta1.cos() - n1 * theta2.cos()) /
                      (n2 * theta1.cos() + n1 * theta2.cos());
            let t_p = (2.0 * n1 * theta1.cos()) /
                      (n2 * theta1.cos() + n1 * theta2.cos());
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
        interpretation: format!("{:.1}% reflected, {:.1}% transmitted",
                               reflectance * 100.0, transmittance * 100.0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Thin Lens Tests (4 tests)
    #[test]
    fn test_thin_lens_converging() {
        let params = OpticsParams {
            focal_length: Some(0.1), // 10 cm converging lens
            object_distance: Some(0.2), // 20 cm
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // 1/di = 1/f - 1/do = 1/0.1 - 1/0.2 = 10 - 5 = 5 → di = 0.2 m
        assert!((secondary["image_distance"] - 0.2).abs() < 0.01);
        assert!(result.interpretation.contains("real"));
    }

    #[test]
    fn test_thin_lens_diverging() {
        let params = OpticsParams {
            focal_length: Some(-0.15), // 15 cm diverging lens
            object_distance: Some(0.3), // 30 cm
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // 1/di = 1/f - 1/do = 1/-0.15 - 1/0.3 = -6.67 - 3.33 = -10 → di = -0.1 m (virtual)
        assert!(secondary["image_distance"] < 0.0);
        assert!(result.interpretation.contains("virtual"));
    }

    #[test]
    fn test_thin_lens_magnification() {
        let params = OpticsParams {
            focal_length: Some(0.05), // 5 cm
            object_distance: Some(0.1), // 10 cm (2f)
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // At 2f, image forms at 2f with magnification = -1
        assert!((secondary["image_distance"] - 0.1).abs() < 0.01);
        assert!((secondary["magnification"].abs() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_thin_lens_virtual_image() {
        let params = OpticsParams {
            focal_length: Some(0.2), // 20 cm
            object_distance: Some(0.1), // 10 cm (inside focal length)
            ..Default::default()
        };

        let result = calculate_thin_lens(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // Inside f gives virtual, upright, enlarged image
        assert!(secondary["image_distance"] < 0.0);
        assert!(secondary["magnification"] > 0.0); // Upright
        assert!(secondary["magnification"].abs() > 1.0); // Enlarged
    }

    // Snell's Law Tests (4 tests)
    #[test]
    fn test_snells_law_air_to_glass() {
        let params = OpticsParams {
            n1: Some(1.0),   // air
            n2: Some(1.5),   // glass
            theta1: Some(30.0),
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        // Should bend toward normal (smaller angle)
        assert!(result.primary_value < 30.0);
        // n1·sin(30°) = 1.5·sin(θ2) → sin(θ2) = 0.5/1.5 = 0.333 → θ2 ≈ 19.47°
        assert!((result.primary_value - 19.47).abs() < 0.5);
    }

    #[test]
    fn test_snells_law_total_internal_reflection() {
        let params = OpticsParams {
            n1: Some(1.5),  // glass
            n2: Some(1.0),  // air
            theta1: Some(50.0), // Greater than critical angle
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        // Should report total internal reflection
        assert!(result.interpretation.contains("Total internal reflection"));
        // Critical angle = arcsin(1.0/1.5) ≈ 41.8°
        assert!((result.primary_value - 41.8).abs() < 1.0);
    }

    #[test]
    fn test_snells_law_reverse_calculation() {
        let params = OpticsParams {
            n1: Some(1.0),
            n2: Some(1.33), // water
            theta2: Some(22.0), // Refracted angle given
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // Calculate incident angle from refracted angle
        // sin(θ1) = 1.33·sin(22°) / 1.0 ≈ 0.498 → θ1 ≈ 29.9°
        assert!(secondary["incident_angle"] > 22.0);
        assert!((secondary["incident_angle"] - 29.9).abs() < 1.0);
    }

    #[test]
    fn test_snells_law_critical_angle() {
        let params = OpticsParams {
            n1: Some(1.5),
            n2: Some(1.0),
            theta1: Some(30.0), // Below critical angle
            ..Default::default()
        };

        let result = calculate_snells_law(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // Critical angle should be calculated
        assert!(secondary.contains_key("critical_angle"));
        // Critical = arcsin(1.0/1.5) ≈ 41.8°
        assert!((secondary["critical_angle"] - 41.8).abs() < 1.0);
    }

    // Diffraction Grating Tests (3 tests)
    #[test]
    fn test_diffraction_first_order() {
        let params = OpticsParams {
            grating_spacing: Some(2e-6), // 2 μm (500 lines/mm)
            wavelength: Some(600e-9),    // 600 nm (red light)
            order: Some(1),
            ..Default::default()
        };

        let result = calculate_diffraction(&params).unwrap();
        // sin(θ) = m·λ/d = 1·600e-9/2e-6 = 0.3 → θ ≈ 17.46°
        assert!((result.primary_value - 17.46).abs() < 1.0);
        assert_eq!(result.unit, "degrees");
    }

    #[test]
    fn test_diffraction_higher_order() {
        let params = OpticsParams {
            grating_spacing: Some(1e-6), // 1 μm
            wavelength: Some(500e-9),    // 500 nm
            order: Some(2),              // Second order
            ..Default::default()
        };

        let result = calculate_diffraction(&params).unwrap();
        // sin(θ) = 2·500e-9/1e-6 = 1.0 → θ = 90° (grazing angle)
        assert!((result.primary_value - 90.0).abs() < 1.0);
        assert!(result.interpretation.contains("Order 2"));
    }

    #[test]
    fn test_diffraction_max_order() {
        let params = OpticsParams {
            grating_spacing: Some(3e-6), // 3 μm
            wavelength: Some(600e-9),    // 600 nm
            order: Some(1),
            ..Default::default()
        };

        let result = calculate_diffraction(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // m_max = d/λ = 3e-6/600e-9 = 5
        assert_eq!(secondary["max_order"] as i32, 5);
    }

    // Fresnel Equations Tests (3 tests)
    #[test]
    fn test_fresnel_s_polarization() {
        let params = OpticsParams {
            n1: Some(1.0),   // air
            n2: Some(1.5),   // glass
            theta1: Some(45.0),
            polarization: Some("s".to_string()),
            ..Default::default()
        };

        let result = calculate_fresnel(&params).unwrap();
        // Should have reflection and transmission coefficients
        let secondary = result.secondary_values.unwrap();
        assert!(secondary.contains_key("reflection_coefficient"));
        assert!(secondary.contains_key("transmission_coefficient"));
        // Sum should be close to 1 (energy conservation)
        let sum = secondary["reflection_coefficient"] + secondary["transmission_coefficient"];
        assert!((sum - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_fresnel_p_polarization() {
        let params = OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(30.0),
            polarization: Some("p".to_string()),
            ..Default::default()
        };

        let result = calculate_fresnel(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // p-polarization should include Brewster's angle
        assert!(secondary.contains_key("brewsters_angle"));
        // Brewster = arctan(n2/n1) = arctan(1.5) ≈ 56.3°
        assert!((secondary["brewsters_angle"] - 56.3).abs() < 1.0);
    }

    #[test]
    fn test_fresnel_brewster_angle() {
        let params = OpticsParams {
            n1: Some(1.0),
            n2: Some(1.5),
            theta1: Some(56.3), // At Brewster's angle
            polarization: Some("p".to_string()),
            ..Default::default()
        };

        let result = calculate_fresnel(&params).unwrap();
        let secondary = result.secondary_values.unwrap();
        // At Brewster's angle, p-polarized reflection should be minimal
        assert!(secondary["reflection_coefficient"] < 0.01);
        // Transmission should be nearly 100%
        assert!(secondary["transmission_coefficient"] > 0.9);
    }
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
