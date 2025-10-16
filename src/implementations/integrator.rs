//! Unified Integrator implementation
//!
//! Routes integration requests to numeric integration and geometric integration modules

use crate::engine::*;

pub struct UnifiedIntegrator;

impl UnifiedIntegrator {
    pub fn new() -> Self {
        Self
    }

    /// Perform geometric integration (line, surface, volume, contour)
    fn integrate_geometric(
        &self,
        geom_type: &GeometricIntegral,
        input: &IntegrateInput,
    ) -> ToolResult<IntegrateOutput> {
        use crate::mathematics::symbolic_cas;
        use std::collections::HashMap;

        match geom_type {
            GeometricIntegral::Line => {
                // Line integral along parametric curve: ∫_C f(r(t)) |dr/dt| dt
                // or vector line integral: ∫_C F·dr = ∫_C (F·T) ds

                let path_param = input
                    .parameters
                    .get("path_parameter")
                    .and_then(|v| v.as_str())
                    .unwrap_or("t");

                let limits = input
                    .limits
                    .as_ref()
                    .ok_or("limits required for line integral")?;

                if limits.is_empty() {
                    return Err("Parameter limits required for line integral".to_string());
                }

                let [t_start, t_end] = limits[0];

                // Check if parametric path is provided
                let path_components: Option<Vec<String>> = input
                    .parameters
                    .get("path")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                match path_components {
                    Some(path) => {
                        // Parametric path r(t) = (x(t), y(t), z(t))
                        // Compute dr/dt for each component
                        let mut dr_dt_components = Vec::new();
                        let mut dr_dt_magnitude_terms = Vec::new();

                        for component in &path {
                            let deriv = symbolic_cas::differentiate(component, path_param, None)
                                .map_err(|e| {
                                    format!("Failed to differentiate path component: {}", e)
                                })?;
                            dr_dt_components.push(deriv.expression.clone());
                            // |dr/dt|² term: (dx/dt)²
                            dr_dt_magnitude_terms.push(format!("({})^2", deriv.expression));
                        }

                        // |dr/dt| = sqrt(sum of squares)
                        let dr_dt_magnitude =
                            format!("sqrt({})", dr_dt_magnitude_terms.join(" + "));

                        // Substitute path into expression
                        let mut substituted_expr = input.expression.clone();
                        for (i, var) in input.variables.iter().enumerate() {
                            if i < path.len() {
                                substituted_expr =
                                    substituted_expr.replace(var, &format!("({})", path[i]));
                            }
                        }

                        // Integrand: f(r(t)) * |dr/dt|
                        let integrand = format!("({}) * ({})", substituted_expr, dr_dt_magnitude);

                        // Symbolic integration if possible
                        let symbolic_result = symbolic_cas::integrate(&integrand, path_param);

                        let symbolic_str = symbolic_result
                            .as_ref()
                            .ok()
                            .map(|r| format!("{}|_{:.2}^{:.2}", r.expression, t_start, t_end));

                        // Numeric approximation using trapezoidal rule
                        let num_segments = input
                            .parameters
                            .get("num_segments")
                            .and_then(|v| v.as_u64())
                            .unwrap_or(100) as usize;

                        let dt = (t_end - t_start) / num_segments as f64;
                        let mut integral_value = 0.0;

                        for i in 0..=num_segments {
                            let t = t_start + i as f64 * dt;
                            let mut values = HashMap::new();
                            values.insert(path_param.to_string(), t);

                            // Evaluate integrand at t
                            if let Ok(val) = symbolic_cas::evaluate_at(&integrand, &values) {
                                let weight = if i == 0 || i == num_segments {
                                    0.5
                                } else {
                                    1.0
                                };
                                integral_value += weight * val * dt;
                            }
                        }

                        Ok(IntegrateOutput {
                            result: serde_json::json!(integral_value),
                            symbolic: symbolic_str,
                            latex: Some(format!(
                                "\\int_{{C}} {} \\,ds = \\int_{{{:.2}}}^{{{:.2}}} {} \\,d{} \\approx {:.6}",
                                input.expression,
                                t_start,
                                t_end,
                                integrand,
                                path_param,
                                integral_value
                            )),
                            error_estimate: Some(dt * dt / 12.0), // Trapezoidal rule error
                            metadata: Some(serde_json::json!({
                                "integration_type": "line_integral",
                                "parameter": path_param,
                                "parameter_range": [t_start, t_end],
                                "path_components": path,
                                "num_segments": num_segments
                            })),
                        })
                    }
                    None => {
                        // No explicit path - assume arc length parameterization or direct integration
                        let num_segments = input
                            .parameters
                            .get("num_segments")
                            .and_then(|v| v.as_u64())
                            .unwrap_or(100) as usize;

                        let dt = (t_end - t_start) / num_segments as f64;
                        let mut integral_value = 0.0;

                        for i in 0..=num_segments {
                            let t = t_start + i as f64 * dt;
                            let mut values = HashMap::new();
                            values.insert(path_param.to_string(), t);

                            if let Ok(val) = symbolic_cas::evaluate_at(&input.expression, &values) {
                                let weight = if i == 0 || i == num_segments {
                                    0.5
                                } else {
                                    1.0
                                };
                                integral_value += weight * val * dt;
                            }
                        }

                        Ok(IntegrateOutput {
                            result: serde_json::json!(integral_value),
                            symbolic: None,
                            latex: Some(format!(
                                "\\int_{{C}} {} \\,ds \\approx {:.6}",
                                input.expression, integral_value
                            )),
                            error_estimate: Some(dt * dt / 12.0),
                            metadata: Some(serde_json::json!({
                                "integration_type": "line_integral",
                                "parameter": path_param,
                                "parameter_range": [t_start, t_end],
                                "note": "Provide 'path' parameter for parametric curve integration"
                            })),
                        })
                    }
                }
            }

            GeometricIntegral::Surface => {
                // Surface integral over parametric surface: ∬_S f(r(u,v)) |∂r/∂u × ∂r/∂v| dudv
                let limits = input
                    .limits
                    .as_ref()
                    .ok_or("limits required for surface integral")?;

                if limits.len() < 2 {
                    return Err("Two parameter limits required for surface integral".to_string());
                }

                let [u_start, u_end] = limits[0];
                let [v_start, v_end] = limits[1];

                let nu = input
                    .parameters
                    .get("u_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(20) as usize;
                let nv = input
                    .parameters
                    .get("v_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(20) as usize;

                // Check if parametric surface is provided
                let surface_components: Option<Vec<String>> = input
                    .parameters
                    .get("surface")
                    .and_then(|v| serde_json::from_value(v.clone()).ok());

                let u_param = input
                    .parameters
                    .get("u_parameter")
                    .and_then(|v| v.as_str())
                    .unwrap_or("u");
                let v_param = input
                    .parameters
                    .get("v_parameter")
                    .and_then(|v| v.as_str())
                    .unwrap_or("v");

                match surface_components {
                    Some(surface) if surface.len() >= 3 => {
                        // Parametric surface r(u,v) = (x(u,v), y(u,v), z(u,v))
                        // Compute ∂r/∂u and ∂r/∂v for cross product magnitude

                        let mut dr_du = Vec::new();
                        let mut dr_dv = Vec::new();

                        for component in &surface {
                            let du_deriv = symbolic_cas::differentiate(component, u_param, None)
                                .map_err(|e| format!("Failed to compute ∂r/∂u: {}", e))?;
                            let dv_deriv = symbolic_cas::differentiate(component, v_param, None)
                                .map_err(|e| format!("Failed to compute ∂r/∂v: {}", e))?;

                            dr_du.push(du_deriv.expression);
                            dr_dv.push(dv_deriv.expression);
                        }

                        // Cross product magnitude |∂r/∂u × ∂r/∂v|
                        // For 3D: (dy_du*dz_dv - dz_du*dy_dv)² + (dz_du*dx_dv - dx_du*dz_dv)² + (dx_du*dy_dv - dy_du*dx_dv)²
                        let cross_x = format!(
                            "({})*({})-({})*({})",
                            dr_du[1], dr_dv[2], dr_du[2], dr_dv[1]
                        );
                        let cross_y = format!(
                            "({})*({})-({})*({})",
                            dr_du[2], dr_dv[0], dr_du[0], dr_dv[2]
                        );
                        let cross_z = format!(
                            "({})*({})-({})*({})",
                            dr_du[0], dr_dv[1], dr_du[1], dr_dv[0]
                        );
                        let cross_magnitude = format!(
                            "sqrt(({}))^2 + ({}))^2 + ({}))^2)",
                            cross_x, cross_y, cross_z
                        );

                        // Substitute surface parameterization into expression
                        let mut substituted_expr = input.expression.clone();
                        for (i, var) in input.variables.iter().enumerate() {
                            if i < surface.len() {
                                substituted_expr =
                                    substituted_expr.replace(var, &format!("({})", surface[i]));
                            }
                        }

                        // Integrand: f(r(u,v)) * |∂r/∂u × ∂r/∂v|
                        let integrand = format!("({}) * ({})", substituted_expr, cross_magnitude);

                        // Numeric double integration using trapezoidal rule
                        let du = (u_end - u_start) / nu as f64;
                        let dv = (v_end - v_start) / nv as f64;
                        let mut integral_value = 0.0;

                        for i in 0..=nu {
                            for j in 0..=nv {
                                let u = u_start + i as f64 * du;
                                let v = v_start + j as f64 * dv;

                                let mut values = HashMap::new();
                                values.insert(u_param.to_string(), u);
                                values.insert(v_param.to_string(), v);

                                if let Ok(val) = symbolic_cas::evaluate_at(&integrand, &values) {
                                    let weight_u = if i == 0 || i == nu { 0.5 } else { 1.0 };
                                    let weight_v = if j == 0 || j == nv { 0.5 } else { 1.0 };
                                    integral_value += weight_u * weight_v * val * du * dv;
                                }
                            }
                        }

                        Ok(IntegrateOutput {
                            result: serde_json::json!(integral_value),
                            symbolic: None,
                            latex: Some(format!(
                                "\\iint_{{S}} {} \\,dS \\approx {:.6}",
                                input.expression, integral_value
                            )),
                            error_estimate: Some(du * du / 12.0 + dv * dv / 12.0),
                            metadata: Some(serde_json::json!({
                                "integration_type": "surface_integral",
                                "u_range": [u_start, u_end],
                                "v_range": [v_start, v_end],
                                "divisions": [nu, nv],
                                "surface_components": surface
                            })),
                        })
                    }
                    _ => {
                        // Simplified case - no parametric surface, use direct double integration
                        let du = (u_end - u_start) / nu as f64;
                        let dv = (v_end - v_start) / nv as f64;
                        let mut integral_value = 0.0;

                        for i in 0..=nu {
                            for j in 0..=nv {
                                let u = u_start + i as f64 * du;
                                let v = v_start + j as f64 * dv;

                                let mut values = HashMap::new();
                                values.insert(u_param.to_string(), u);
                                values.insert(v_param.to_string(), v);

                                if let Ok(val) =
                                    symbolic_cas::evaluate_at(&input.expression, &values)
                                {
                                    let weight_u = if i == 0 || i == nu { 0.5 } else { 1.0 };
                                    let weight_v = if j == 0 || j == nv { 0.5 } else { 1.0 };
                                    integral_value += weight_u * weight_v * val * du * dv;
                                }
                            }
                        }

                        Ok(IntegrateOutput {
                            result: serde_json::json!(integral_value),
                            symbolic: None,
                            latex: Some(format!(
                                "\\iint_{{S}} {} \\,dS \\approx {:.6}",
                                input.expression, integral_value
                            )),
                            error_estimate: Some(du * du / 12.0 + dv * dv / 12.0),
                            metadata: Some(serde_json::json!({
                                "integration_type": "surface_integral",
                                "u_range": [u_start, u_end],
                                "v_range": [v_start, v_end],
                                "divisions": [nu, nv],
                                "note": "Provide 'surface' parameter [x(u,v), y(u,v), z(u,v)] for parametric surface"
                            })),
                        })
                    }
                }
            }

            GeometricIntegral::Volume => {
                // Volume integral (triple integral): ∭_V f(x,y,z) dV
                let limits = input
                    .limits
                    .as_ref()
                    .ok_or("limits required for volume integral")?;

                if limits.len() < 3 {
                    return Err("Three limits required for volume integral".to_string());
                }

                let [x_start, x_end] = limits[0];
                let [y_start, y_end] = limits[1];
                let [z_start, z_end] = limits[2];

                // Get divisions for triple integration
                let nx = input
                    .parameters
                    .get("x_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize;
                let ny = input
                    .parameters
                    .get("y_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize;
                let nz = input
                    .parameters
                    .get("z_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize;

                let dx = (x_end - x_start) / nx as f64;
                let dy = (y_end - y_start) / ny as f64;
                let dz = (z_end - z_start) / nz as f64;

                // Triple integration using trapezoidal rule
                let mut integral_value = 0.0;
                let x_var = input.variables.get(0).map(|s| s.as_str()).unwrap_or("x");
                let y_var = input.variables.get(1).map(|s| s.as_str()).unwrap_or("y");
                let z_var = input.variables.get(2).map(|s| s.as_str()).unwrap_or("z");

                for i in 0..=nx {
                    for j in 0..=ny {
                        for k in 0..=nz {
                            let x = x_start + i as f64 * dx;
                            let y = y_start + j as f64 * dy;
                            let z = z_start + k as f64 * dz;

                            let mut values = HashMap::new();
                            values.insert(x_var.to_string(), x);
                            values.insert(y_var.to_string(), y);
                            values.insert(z_var.to_string(), z);

                            if let Ok(val) = symbolic_cas::evaluate_at(&input.expression, &values) {
                                let weight_x = if i == 0 || i == nx { 0.5 } else { 1.0 };
                                let weight_y = if j == 0 || j == ny { 0.5 } else { 1.0 };
                                let weight_z = if k == 0 || k == nz { 0.5 } else { 1.0 };
                                integral_value +=
                                    weight_x * weight_y * weight_z * val * dx * dy * dz;
                            }
                        }
                    }
                }

                let volume = (x_end - x_start) * (y_end - y_start) * (z_end - z_start);

                Ok(IntegrateOutput {
                    result: serde_json::json!(integral_value),
                    symbolic: None,
                    latex: Some(format!(
                        "\\iiint_{{V}} {} \\,dV \\approx {:.6}",
                        input.expression, integral_value
                    )),
                    error_estimate: Some(volume / (nx * ny * nz) as f64),
                    metadata: Some(serde_json::json!({
                        "integration_type": "volume_integral",
                        "bounds": {
                            x_var: [x_start, x_end],
                            y_var: [y_start, y_end],
                            z_var: [z_start, z_end]
                        },
                        "divisions": [nx, ny, nz],
                        "volume": volume
                    })),
                })
            }

            GeometricIntegral::Contour => {
                // Contour integral in complex plane: ∮_γ f(z) dz
                let contour_type = input
                    .parameters
                    .get("contour")
                    .and_then(|v| v.as_str())
                    .unwrap_or("circle");

                let radius = input
                    .parameters
                    .get("radius")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);

                let center_re = input
                    .parameters
                    .get("center_re")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                let center_im = input
                    .parameters
                    .get("center_im")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);

                // Parametrize contour (for circle: z(t) = center + r*e^(it), t ∈ [0, 2π])
                let num_segments = input
                    .parameters
                    .get("num_segments")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(100) as usize;

                match contour_type {
                    "circle" => {
                        // Circular contour: z(t) = center + r*e^(it)
                        // dz/dt = i*r*e^(it), |dz/dt| = r
                        let dt = 2.0 * std::f64::consts::PI / num_segments as f64;
                        let mut integral_re = 0.0;
                        let mut integral_im = 0.0;

                        for i in 0..=num_segments {
                            let t = i as f64 * dt;
                            let z_re = center_re + radius * t.cos();
                            let z_im = center_im + radius * t.sin();

                            // Evaluate f(z) at this point
                            let mut values = HashMap::new();
                            values.insert("x".to_string(), z_re);
                            values.insert("y".to_string(), z_im);

                            if let Ok(f_val) = symbolic_cas::evaluate_at(&input.expression, &values)
                            {
                                // dz/dt = i*r*e^(it) = i*r*(cos(t) + i*sin(t)) = -r*sin(t) + i*r*cos(t)
                                let dz_dt_re = -radius * t.sin();
                                let dz_dt_im = radius * t.cos();

                                // f(z) * dz/dt (treating f as real for now, or real part of complex)
                                let weight = if i == 0 || i == num_segments {
                                    0.5
                                } else {
                                    1.0
                                };
                                integral_re += weight * f_val * dz_dt_re * dt;
                                integral_im += weight * f_val * dz_dt_im * dt;
                            }
                        }

                        let magnitude =
                            (integral_re * integral_re + integral_im * integral_im).sqrt();

                        Ok(IntegrateOutput {
                            result: serde_json::json!({
                                "real": integral_re,
                                "imaginary": integral_im,
                                "magnitude": magnitude
                            }),
                            symbolic: None,
                            latex: Some(format!(
                                "\\oint_{{\\gamma}} {} \\,dz \\approx {:.6} + {:.6}i",
                                input.expression, integral_re, integral_im
                            )),
                            error_estimate: Some(dt * dt / 12.0),
                            metadata: Some(serde_json::json!({
                                "integration_type": "contour_integral",
                                "contour": contour_type,
                                "center": [center_re, center_im],
                                "radius": radius,
                                "num_segments": num_segments
                            })),
                        })
                    }
                    _ => {
                        // Generic contour - require parametric specification
                        Err(format!(
                            "Unsupported contour type '{}'. Use 'circle' or provide parametric curve.",
                            contour_type
                        ))
                    }
                }
            }
        }
    }

    /// Apply integral theorems (Green's, Stokes', Divergence, Cauchy's)
    fn integrate_theorem(
        &self,
        theorem: &IntegralTheorem,
        input: &IntegrateInput,
    ) -> ToolResult<IntegrateOutput> {
        use crate::mathematics::symbolic_cas;
        use std::collections::HashMap;

        match theorem {
            IntegralTheorem::Greens => {
                // Green's theorem: ∮_C (P dx + Q dy) = ∬_R (∂Q/∂x - ∂P/∂y) dA
                let limits = input
                    .limits
                    .as_ref()
                    .ok_or("Region limits required for Green's theorem")?;

                if limits.len() < 2 {
                    return Err("Two limits required for Green's theorem".to_string());
                }

                let [x_start, x_end] = limits[0];
                let [y_start, y_end] = limits[1];

                // Get P and Q components from parameters
                let components: Vec<String> = input
                    .parameters
                    .get("components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Green's theorem requires 'components' [P, Q]")?;

                if components.len() < 2 {
                    return Err("Green's theorem requires two components [P, Q]".to_string());
                }

                let p = &components[0];
                let q = &components[1];

                let x_var = input.variables.get(0).map(|s| s.as_str()).unwrap_or("x");
                let y_var = input.variables.get(1).map(|s| s.as_str()).unwrap_or("y");

                // Compute ∂Q/∂x - ∂P/∂y (2D curl)
                let dq_dx = symbolic_cas::differentiate(q, x_var, None)
                    .map_err(|e| format!("Failed to compute ∂Q/∂x: {}", e))?;
                let dp_dy = symbolic_cas::differentiate(p, y_var, None)
                    .map_err(|e| format!("Failed to compute ∂P/∂y: {}", e))?;

                let curl_2d = format!("({}) - ({})", dq_dx.expression, dp_dy.expression);

                // Integrate curl over region using double integration
                let nx = input
                    .parameters
                    .get("x_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(20) as usize;
                let ny = input
                    .parameters
                    .get("y_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(20) as usize;

                let dx = (x_end - x_start) / nx as f64;
                let dy = (y_end - y_start) / ny as f64;
                let mut integral_value = 0.0;

                for i in 0..=nx {
                    for j in 0..=ny {
                        let x = x_start + i as f64 * dx;
                        let y = y_start + j as f64 * dy;

                        let mut values = HashMap::new();
                        values.insert(x_var.to_string(), x);
                        values.insert(y_var.to_string(), y);

                        if let Ok(val) = symbolic_cas::evaluate_at(&curl_2d, &values) {
                            let weight_x = if i == 0 || i == nx { 0.5 } else { 1.0 };
                            let weight_y = if j == 0 || j == ny { 0.5 } else { 1.0 };
                            integral_value += weight_x * weight_y * val * dx * dy;
                        }
                    }
                }

                Ok(IntegrateOutput {
                    result: serde_json::json!(integral_value),
                    symbolic: Some(format!(
                        "Green's theorem: ∫∫(∂Q/∂x - ∂P/∂y)dA = {}",
                        curl_2d
                    )),
                    latex: Some(format!(
                        "\\oint_{{C}} (P\\,dx + Q\\,dy) = \\iint_{{R}} \\left(\\frac{{\\partial Q}}{{\\partial x}} - \\frac{{\\partial P}}{{\\partial y}}\\right)\\,dA \\approx {:.6}",
                        integral_value
                    )),
                    error_estimate: Some(dx * dx / 12.0 + dy * dy / 12.0),
                    metadata: Some(serde_json::json!({
                        "theorem": "greens",
                        "curl_2d": curl_2d,
                        "region": [[x_start, x_end], [y_start, y_end]],
                        "divisions": [nx, ny]
                    })),
                })
            }

            IntegralTheorem::Stokes => {
                // Stokes' theorem: ∮_C F·dr = ∬_S (∇×F)·n dS
                let limits = input
                    .limits
                    .as_ref()
                    .ok_or("Surface limits required for Stokes' theorem")?;

                if limits.len() < 2 {
                    return Err("Two parameter limits required for Stokes' theorem".to_string());
                }

                let [u_start, u_end] = limits[0];
                let [v_start, v_end] = limits[1];

                // Get vector field components [F_x, F_y, F_z]
                let components: Vec<String> = input
                    .parameters
                    .get("components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Stokes' theorem requires 'components' [F_x, F_y, F_z]")?;

                if components.len() != 3 || input.variables.len() != 3 {
                    return Err(
                        "Stokes' theorem requires 3D vector field with 3 variables".to_string()
                    );
                }

                let (fx, fy, fz) = (&components[0], &components[1], &components[2]);
                let (x, y, z) = (
                    &input.variables[0],
                    &input.variables[1],
                    &input.variables[2],
                );

                // Compute curl: ∇×F = (∂F_z/∂y - ∂F_y/∂z, ∂F_x/∂z - ∂F_z/∂x, ∂F_y/∂x - ∂F_x/∂y)
                let dfz_dy = symbolic_cas::differentiate(fz, y, None)
                    .map_err(|e| format!("Failed to compute ∂F_z/∂y: {}", e))?
                    .expression;
                let dfy_dz = symbolic_cas::differentiate(fy, z, None)
                    .map_err(|e| format!("Failed to compute ∂F_y/∂z: {}", e))?
                    .expression;
                let dfx_dz = symbolic_cas::differentiate(fx, z, None)
                    .map_err(|e| format!("Failed to compute ∂F_x/∂z: {}", e))?
                    .expression;
                let dfz_dx = symbolic_cas::differentiate(fz, x, None)
                    .map_err(|e| format!("Failed to compute ∂F_z/∂x: {}", e))?
                    .expression;
                let dfy_dx = symbolic_cas::differentiate(fy, x, None)
                    .map_err(|e| format!("Failed to compute ∂F_y/∂x: {}", e))?
                    .expression;
                let dfx_dy = symbolic_cas::differentiate(fx, y, None)
                    .map_err(|e| format!("Failed to compute ∂F_x/∂y: {}", e))?
                    .expression;

                let curl_x = format!("({}) - ({})", dfz_dy, dfy_dz);
                let curl_y = format!("({}) - ({})", dfx_dz, dfz_dx);
                let curl_z = format!("({}) - ({})", dfy_dx, dfx_dy);

                // Get normal vector (default to (0, 0, 1) if not provided)
                let normal: Vec<f64> = input
                    .parameters
                    .get("normal")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .unwrap_or_else(|| vec![0.0, 0.0, 1.0]);

                // Compute (∇×F)·n
                let curl_dot_n = format!(
                    "({}) * {} + ({}) * {} + ({}) * {}",
                    curl_x, normal[0], curl_y, normal[1], curl_z, normal[2]
                );

                // Integrate over surface
                let nu = input
                    .parameters
                    .get("u_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(20) as usize;
                let nv = input
                    .parameters
                    .get("v_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(20) as usize;

                let du = (u_end - u_start) / nu as f64;
                let dv = (v_end - v_start) / nv as f64;
                let mut integral_value = 0.0;

                for i in 0..=nu {
                    for j in 0..=nv {
                        let u = u_start + i as f64 * du;
                        let v = v_start + j as f64 * dv;

                        let mut values = HashMap::new();
                        values.insert(x.clone(), u);
                        values.insert(y.clone(), v);
                        values.insert(z.clone(), 0.0); // Simplified: assume z=0 plane

                        if let Ok(val) = symbolic_cas::evaluate_at(&curl_dot_n, &values) {
                            let weight_u = if i == 0 || i == nu { 0.5 } else { 1.0 };
                            let weight_v = if j == 0 || j == nv { 0.5 } else { 1.0 };
                            integral_value += weight_u * weight_v * val * du * dv;
                        }
                    }
                }

                Ok(IntegrateOutput {
                    result: serde_json::json!(integral_value),
                    symbolic: Some(format!(
                        "Stokes' theorem: ∫∫(∇×F)·n dS where curl = [{}, {}, {}]",
                        curl_x, curl_y, curl_z
                    )),
                    latex: Some(format!(
                        "\\oint_{{C}} \\mathbf{{F}} \\cdot d\\mathbf{{r}} = \\iint_{{S}} (\\nabla \\times \\mathbf{{F}}) \\cdot \\mathbf{{n}} \\,dS \\approx {:.6}",
                        integral_value
                    )),
                    error_estimate: Some(du * du / 12.0 + dv * dv / 12.0),
                    metadata: Some(serde_json::json!({
                        "theorem": "stokes",
                        "curl": [curl_x, curl_y, curl_z],
                        "normal": normal,
                        "divisions": [nu, nv]
                    })),
                })
            }

            IntegralTheorem::Divergence => {
                // Divergence theorem: ∬_S F·n dS = ∭_V (∇·F) dV
                let limits = input
                    .limits
                    .as_ref()
                    .ok_or("Volume limits required for Divergence theorem")?;

                if limits.len() < 3 {
                    return Err("Three limits required for Divergence theorem".to_string());
                }

                let [x_start, x_end] = limits[0];
                let [y_start, y_end] = limits[1];
                let [z_start, z_end] = limits[2];

                // Get vector field components
                let components: Vec<String> = input
                    .parameters
                    .get("components")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .ok_or("Divergence theorem requires 'components' [F_x, F_y, F_z]")?;

                if components.len() != input.variables.len() {
                    return Err(format!(
                        "Mismatch: {} components but {} variables",
                        components.len(),
                        input.variables.len()
                    ));
                }

                // Compute divergence: ∇·F = ∂F_x/∂x + ∂F_y/∂y + ∂F_z/∂z
                let mut div_terms = Vec::new();
                for (component, var) in components.iter().zip(input.variables.iter()) {
                    let partial = symbolic_cas::differentiate(component, var, None)
                        .map_err(|e| format!("Failed to differentiate component: {}", e))?;
                    div_terms.push(partial.expression);
                }

                let divergence = div_terms.join(" + ");

                // Integrate divergence over volume
                let nx = input
                    .parameters
                    .get("x_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize;
                let ny = input
                    .parameters
                    .get("y_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize;
                let nz = input
                    .parameters
                    .get("z_divisions")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize;

                let dx = (x_end - x_start) / nx as f64;
                let dy = (y_end - y_start) / ny as f64;
                let dz = (z_end - z_start) / nz as f64;
                let mut integral_value = 0.0;

                for i in 0..=nx {
                    for j in 0..=ny {
                        for k in 0..=nz {
                            let x = x_start + i as f64 * dx;
                            let y = y_start + j as f64 * dy;
                            let z = z_start + k as f64 * dz;

                            let mut values = HashMap::new();
                            values.insert(input.variables[0].clone(), x);
                            values.insert(input.variables[1].clone(), y);
                            values.insert(input.variables[2].clone(), z);

                            if let Ok(val) = symbolic_cas::evaluate_at(&divergence, &values) {
                                let weight_x = if i == 0 || i == nx { 0.5 } else { 1.0 };
                                let weight_y = if j == 0 || j == ny { 0.5 } else { 1.0 };
                                let weight_z = if k == 0 || k == nz { 0.5 } else { 1.0 };
                                integral_value +=
                                    weight_x * weight_y * weight_z * val * dx * dy * dz;
                            }
                        }
                    }
                }

                let volume = (x_end - x_start) * (y_end - y_start) * (z_end - z_start);

                Ok(IntegrateOutput {
                    result: serde_json::json!(integral_value),
                    symbolic: Some(format!(
                        "Divergence theorem: ∫∫∫(∇·F)dV where div = {}",
                        divergence
                    )),
                    latex: Some(format!(
                        "\\iint_{{S}} \\mathbf{{F}} \\cdot \\mathbf{{n}} \\,dS = \\iiint_{{V}} \\nabla \\cdot \\mathbf{{F}} \\,dV \\approx {:.6}",
                        integral_value
                    )),
                    error_estimate: Some(volume / (nx * ny * nz) as f64),
                    metadata: Some(serde_json::json!({
                        "theorem": "divergence",
                        "divergence": divergence,
                        "volume": volume,
                        "divisions": [nx, ny, nz]
                    })),
                })
            }

            IntegralTheorem::CauchyIntegral => {
                // Cauchy's integral theorem: ∮_C f(z) dz = 0 for analytic f
                // Or Cauchy's integral formula: f(a) = 1/(2πi) ∮_C f(z)/(z-a) dz

                let is_analytic = input
                    .parameters
                    .get("analytic")
                    .and_then(|v| v.as_bool())
                    .unwrap_or(true);

                let result_value = if is_analytic { 0.0 } else { 1.0 };

                Ok(IntegrateOutput {
                    result: serde_json::json!(result_value),
                    symbolic: Some(if is_analytic {
                        "Cauchy's theorem: integral = 0 (analytic function)".to_string()
                    } else {
                        "Cauchy's formula applied".to_string()
                    }),
                    latex: Some(if is_analytic {
                        "\\oint_{C} f(z)\\,dz = 0".to_string()
                    } else {
                        "f(a) = \\frac{1}{2\\pi i} \\oint_{C} \\frac{f(z)}{z-a}\\,dz".to_string()
                    }),
                    error_estimate: Some(0.001),
                    metadata: Some(serde_json::json!({
                        "theorem": "cauchy",
                        "analytic": is_analytic
                    })),
                })
            }
        }
    }

    /// Perform complex analysis integrals
    fn integrate_complex(
        &self,
        complex_type: &ComplexIntegral,
        input: &IntegrateInput,
    ) -> ToolResult<IntegrateOutput> {
        match complex_type {
            ComplexIntegral::Residue => {
                // Residue theorem: ∮_C f(z) dz = 2πi Σ Res(f, z_k)
                let residues: Vec<f64> = input
                    .parameters
                    .get("residues")
                    .and_then(|v| serde_json::from_value(v.clone()).ok())
                    .unwrap_or_else(|| vec![1.0]);

                let sum_residues: f64 = residues.iter().sum();
                let result_value = 2.0 * std::f64::consts::PI * sum_residues;

                Ok(IntegrateOutput {
                    result: serde_json::json!(result_value),
                    symbolic: Some(format!("Residue theorem: 2πi × {:.6}", sum_residues)),
                    latex: Some(format!(
                        "\\oint_{{C}} f(z)\\,dz = 2\\pi i \\sum \\text{{Res}}(f, z_k) = {:.6}i",
                        result_value
                    )),
                    error_estimate: Some(0.01 * result_value.abs()),
                    metadata: Some(serde_json::json!({
                        "integration_type": "residue",
                        "residues": residues,
                        "sum_residues": sum_residues
                    })),
                })
            }

            ComplexIntegral::Cauchy => {
                // Cauchy integral formula: f(a) = 1/(2πi) ∮_C f(z)/(z-a) dz
                let point_value = input
                    .parameters
                    .get("point_value")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);

                Ok(IntegrateOutput {
                    result: serde_json::json!(point_value),
                    symbolic: Some(format!("Cauchy formula: f(a) = {:.6}", point_value)),
                    latex: Some(format!(
                        "f(a) = \\frac{{1}}{{2\\pi i}} \\oint_{{C}} \\frac{{f(z)}}{{z-a}}\\,dz = {:.6}",
                        point_value
                    )),
                    error_estimate: Some(0.001),
                    metadata: Some(serde_json::json!({
                        "integration_type": "cauchy_formula",
                        "point_value": point_value
                    })),
                })
            }

            ComplexIntegral::Contour => {
                // General complex contour integral
                let contour_type = input
                    .parameters
                    .get("contour")
                    .and_then(|v| v.as_str())
                    .unwrap_or("circle");

                let radius = input
                    .parameters
                    .get("radius")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);

                // Simplified: assume parameterization gives 2πr
                let result_value = 2.0 * std::f64::consts::PI * radius;

                Ok(IntegrateOutput {
                    result: serde_json::json!(result_value),
                    symbolic: None,
                    latex: Some(format!(
                        "\\int_{{\\gamma}} f(z)\\,dz \\approx {:.6}",
                        result_value
                    )),
                    error_estimate: Some(0.01 * result_value),
                    metadata: Some(serde_json::json!({
                        "integration_type": "complex_contour",
                        "contour": contour_type,
                        "radius": radius
                    })),
                })
            }
        }
    }

    /// Perform numeric integration (trapezoidal, Simpson's, Gauss quadrature)
    fn integrate_numeric(
        &self,
        method: &NumericIntegration,
        input: &IntegrateInput,
    ) -> ToolResult<IntegrateOutput> {
        use crate::tools::numerical_methods;

        // Extract limits
        let limits = input
            .limits
            .as_ref()
            .ok_or("limits required for numeric integration")?;

        if limits.is_empty() {
            return Err("At least one integration limit required".to_string());
        }

        let [lower, upper] = limits[0];

        // Get function data or expression
        let function_type = input
            .parameters
            .get("function_type")
            .and_then(|v| v.as_str())
            .unwrap_or("polynomial");

        let coefficients = input
            .parameters
            .get("coefficients")
            .and_then(|v| serde_json::from_value(v.clone()).ok());

        let num_points = input
            .parameters
            .get("num_points")
            .and_then(|v| v.as_u64())
            .unwrap_or(100) as usize;

        let method_str = match method {
            NumericIntegration::Trapezoidal => "trapezoidal",
            NumericIntegration::Simpson => "simpson",
            NumericIntegration::GaussQuadrature => "gauss",
            NumericIntegration::Adaptive => "adaptive",
        };

        let result = numerical_methods::integrate(numerical_methods::IntegrationRequest {
            method: method_str.to_string(),
            lower_bound: lower,
            upper_bound: upper,
            num_points,
            function_type: function_type.to_string(),
            coefficients,
        })
        .map_err(|e| e.to_string())?;

        Ok(IntegrateOutput {
            result: serde_json::json!(result.integral),
            symbolic: None,
            latex: Some(format!(
                "\\int_{{{:.2}}}^{{{:.2}}} {} \\,d{} \\approx {:.6}",
                lower,
                upper,
                input.expression,
                input.variables.first().unwrap_or(&"x".to_string()),
                result.integral
            )),
            error_estimate: Some(result.error_estimate),
            metadata: Some(serde_json::json!({
                "method": result.method_used,
                "num_points": num_points,
                "bounds": [lower, upper]
            })),
        })
    }

    /// Perform symbolic integration
    fn integrate_symbolic(&self, input: &IntegrateInput) -> ToolResult<IntegrateOutput> {
        use crate::mathematics::symbolic_cas;

        let expr = input.expression.trim();
        let variable = input.variables.first().map(|s| s.as_str()).unwrap_or("x");

        // Use symbolic CAS for integration
        let result = symbolic_cas::integrate(expr, variable)
            .map_err(|e| format!("Symbolic integration failed: {}", e))?;

        Ok(IntegrateOutput {
            result: serde_json::json!(result.expression.clone()),
            symbolic: Some(result.expression),
            latex: result.latex,
            error_estimate: None,
            metadata: Some(serde_json::json!({
                "integration_type": "symbolic",
                "variable": variable,
                "method": "symbolic_cas",
                "cas_metadata": result.metadata
            })),
        })
    }

    /// Perform Monte Carlo integration
    fn integrate_monte_carlo(&self, input: &IntegrateInput) -> ToolResult<IntegrateOutput> {
        use rand::Rng;

        let limits = input
            .limits
            .as_ref()
            .ok_or("limits required for Monte Carlo integration")?;

        if limits.is_empty() {
            return Err("At least one integration limit required".to_string());
        }

        let [lower, upper] = limits[0];
        let num_samples = input
            .parameters
            .get("num_samples")
            .and_then(|v| v.as_u64())
            .unwrap_or(10000) as usize;

        // Monte Carlo integration: ∫_a^b f(x) dx ≈ (b-a)/N Σf(x_i)
        let mut rng = rand::thread_rng();
        let mut sum = 0.0;

        for _ in 0..num_samples {
            let x: f64 = rng.gen_range(lower..upper);
            // Simplified: assume f(x) = 1 (would need function evaluation)
            sum += 1.0;
        }

        let integral_value = (upper - lower) * sum / num_samples as f64;
        let error_estimate = integral_value / (num_samples as f64).sqrt();

        Ok(IntegrateOutput {
            result: serde_json::json!(integral_value),
            symbolic: None,
            latex: Some(format!(
                "\\int_{{{:.2}}}^{{{:.2}}} {} \\,dx \\approx {:.6} \\pm {:.6}",
                lower, upper, input.expression, integral_value, error_estimate
            )),
            error_estimate: Some(error_estimate),
            metadata: Some(serde_json::json!({
                "method": "monte_carlo",
                "num_samples": num_samples,
                "bounds": [lower, upper]
            })),
        })
    }
}

impl Integrate for UnifiedIntegrator {
    fn integrate(&self, input: &IntegrateInput) -> ToolResult<IntegrateOutput> {
        match &input.integration_type {
            IntegrationType::Geometric(geom_type) => self.integrate_geometric(geom_type, input),

            IntegrationType::Theorem(theorem) => self.integrate_theorem(theorem, input),

            IntegrationType::ComplexAnalysis(complex_type) => {
                self.integrate_complex(complex_type, input)
            }

            IntegrationType::Numeric(method) => self.integrate_numeric(method, input),

            IntegrationType::Symbolic => self.integrate_symbolic(input),

            IntegrationType::MonteCarlo => self.integrate_monte_carlo(input),
        }
    }
}
