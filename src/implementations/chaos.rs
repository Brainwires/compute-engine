//! Unified Chaos Theory implementation
//!
//! Routes chaos theory requests to fractal, attractor, and Lyapunov modules
//! in mathematics/chaos/

use crate::engine::*;

pub struct UnifiedChaos;

impl UnifiedChaos {
    pub fn new() -> Self {
        Self
    }

    /// Handle fractal operations
    fn handle_fractal(&self, fractal_type: &FractalType, input: &ChaosInput) -> ToolResult<ChaosOutput> {
        use crate::mathematics::chaos::{fractals, Complex};

        let max_iter = input.iterations.unwrap_or(100);
        let escape_radius = input.parameters
            .get("escape_radius")
            .and_then(|v| v.as_f64())
            .unwrap_or(2.0);

        match fractal_type {
            FractalType::Mandelbrot => {
                // Get the complex point c
                let c_re = input.parameters
                    .get("c_re")
                    .or_else(|| input.parameters.get("re"))
                    .and_then(|v| v.as_f64())
                    .unwrap_or(-0.5);
                let c_im = input.parameters
                    .get("c_im")
                    .or_else(|| input.parameters.get("im"))
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);

                let c = Complex::new(c_re, c_im);
                let result = fractals::mandelbrot(c, max_iter, escape_radius);

                let mut values = std::collections::HashMap::new();
                values.insert("iterations".to_string(), result.iterations as f64);
                values.insert("escaped".to_string(), if result.escaped { 1.0 } else { 0.0 });
                values.insert("final_magnitude".to_string(), result.final_magnitude);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "iterations": result.iterations,
                        "escaped": result.escaped,
                        "final_magnitude": result.final_magnitude,
                        "in_set": !result.escaped
                    }),
                    trajectory: None,
                    values: Some(values),
                    metadata: Some(serde_json::json!({
                        "fractal": "mandelbrot",
                        "c": {"re": c_re, "im": c_im},
                        "max_iter": max_iter,
                        "escape_radius": escape_radius
                    })),
                })
            }

            FractalType::Julia => {
                // Julia set: z_0 varies, c is fixed
                let z_re = input.parameters
                    .get("z_re")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                let z_im = input.parameters
                    .get("z_im")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                let c_re = input.parameters
                    .get("c_re")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(-0.4);
                let c_im = input.parameters
                    .get("c_im")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.6);

                let z0 = Complex::new(z_re, z_im);
                let c = Complex::new(c_re, c_im);
                let result = fractals::julia(z0, c, max_iter, escape_radius);

                let mut values = std::collections::HashMap::new();
                values.insert("iterations".to_string(), result.iterations as f64);
                values.insert("escaped".to_string(), if result.escaped { 1.0 } else { 0.0 });
                values.insert("final_magnitude".to_string(), result.final_magnitude);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "iterations": result.iterations,
                        "escaped": result.escaped,
                        "final_magnitude": result.final_magnitude,
                        "in_set": !result.escaped
                    }),
                    trajectory: None,
                    values: Some(values),
                    metadata: Some(serde_json::json!({
                        "fractal": "julia",
                        "z0": {"re": z_re, "im": z_im},
                        "c": {"re": c_re, "im": c_im}
                    })),
                })
            }

            FractalType::BurningShip => {
                let c_re = input.parameters
                    .get("c_re")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(-0.5);
                let c_im = input.parameters
                    .get("c_im")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(-0.5);

                let c = Complex::new(c_re, c_im);
                let result = fractals::burning_ship(c, max_iter, escape_radius);

                let mut values = std::collections::HashMap::new();
                values.insert("iterations".to_string(), result.iterations as f64);
                values.insert("escaped".to_string(), if result.escaped { 1.0 } else { 0.0 });

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "iterations": result.iterations,
                        "escaped": result.escaped,
                        "final_magnitude": result.final_magnitude
                    }),
                    trajectory: None,
                    values: Some(values),
                    metadata: Some(serde_json::json!({
                        "fractal": "burning_ship",
                        "c": {"re": c_re, "im": c_im}
                    })),
                })
            }

            FractalType::KochSnowflake => {
                let order = input.iterations.unwrap_or(4).min(8); // Limit order to prevent huge output
                let points = fractals::koch_snowflake(order);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "points": points,
                        "n_points": points.len()
                    }),
                    trajectory: Some(points.iter().map(|&(x, y)| vec![x, y]).collect()),
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("n_points".to_string(), points.len() as f64);
                        v.insert("order".to_string(), order as f64);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "fractal": "koch_snowflake",
                        "order": order,
                        "theoretical_dimension": 1.2619 // log(4)/log(3)
                    })),
                })
            }

            FractalType::SierpinskiTriangle => {
                Err("Sierpinski Triangle not yet implemented".to_string())
            }

            FractalType::DragonCurve => {
                Err("Dragon Curve not yet implemented".to_string())
            }

            FractalType::Cantor => {
                Err("Cantor set not yet implemented".to_string())
            }
        }
    }

    /// Handle attractor operations
    fn handle_attractor(&self, attractor_type: &AttractorType, input: &ChaosInput) -> ToolResult<ChaosOutput> {
        use crate::mathematics::chaos::{attractors, Point3D};

        let steps = input.iterations.unwrap_or(10000);
        let dt = input.parameters
            .get("dt")
            .and_then(|v| v.as_f64())
            .unwrap_or(0.01);

        let initial = Point3D::new(
            input.parameters.get("x0").and_then(|v| v.as_f64()).unwrap_or(1.0),
            input.parameters.get("y0").and_then(|v| v.as_f64()).unwrap_or(1.0),
            input.parameters.get("z0").and_then(|v| v.as_f64()).unwrap_or(1.0),
        );

        match attractor_type {
            AttractorType::Lorenz => {
                let config = attractors::LorenzConfig {
                    sigma: input.parameters.get("sigma").and_then(|v| v.as_f64()).unwrap_or(10.0),
                    rho: input.parameters.get("rho").and_then(|v| v.as_f64()).unwrap_or(28.0),
                    beta: input.parameters.get("beta").and_then(|v| v.as_f64()).unwrap_or(8.0 / 3.0),
                };

                let result = attractors::lorenz_attractor(initial, &config, dt, steps);

                let trajectory: Vec<Vec<f64>> = result.points.iter()
                    .map(|p| vec![p.x, p.y, p.z])
                    .collect();

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "n_points": result.points.len(),
                        "times": result.times.first().zip(result.times.last())
                            .map(|(&t0, &tn)| serde_json::json!({"start": t0, "end": tn}))
                    }),
                    trajectory: Some(trajectory),
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("sigma".to_string(), config.sigma);
                        v.insert("rho".to_string(), config.rho);
                        v.insert("beta".to_string(), config.beta);
                        v.insert("n_steps".to_string(), steps as f64);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "attractor": "lorenz",
                        "config": {"sigma": config.sigma, "rho": config.rho, "beta": config.beta},
                        "dt": dt,
                        "steps": steps,
                        "initial": {"x": initial.x, "y": initial.y, "z": initial.z}
                    })),
                })
            }

            AttractorType::Rossler => {
                let config = attractors::RosslerConfig {
                    a: input.parameters.get("a").and_then(|v| v.as_f64()).unwrap_or(0.2),
                    b: input.parameters.get("b").and_then(|v| v.as_f64()).unwrap_or(0.2),
                    c: input.parameters.get("c").and_then(|v| v.as_f64()).unwrap_or(5.7),
                };

                let result = attractors::rossler_attractor(initial, &config, dt, steps);

                let trajectory: Vec<Vec<f64>> = result.points.iter()
                    .map(|p| vec![p.x, p.y, p.z])
                    .collect();

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "n_points": result.points.len()
                    }),
                    trajectory: Some(trajectory),
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("a".to_string(), config.a);
                        v.insert("b".to_string(), config.b);
                        v.insert("c".to_string(), config.c);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "attractor": "rossler",
                        "config": {"a": config.a, "b": config.b, "c": config.c}
                    })),
                })
            }

            AttractorType::Henon => {
                Err("Henon attractor not yet implemented".to_string())
            }

            AttractorType::Chua => {
                Err("Chua attractor not yet implemented".to_string())
            }

            AttractorType::Thomas => {
                Err("Thomas attractor not yet implemented".to_string())
            }

            AttractorType::Aizawa => {
                Err("Aizawa attractor not yet implemented".to_string())
            }

            AttractorType::Chen => {
                Err("Chen attractor not yet implemented".to_string())
            }
        }
    }

    /// Handle Lyapunov exponent operations
    fn handle_lyapunov(&self, method: &LyapunovMethod, input: &ChaosInput) -> ToolResult<ChaosOutput> {
        use crate::mathematics::chaos::lyapunov;

        let iterations = input.iterations.unwrap_or(1000);
        let transient = input.parameters
            .get("transient")
            .and_then(|v| v.as_u64())
            .map(|v| v as usize)
            .unwrap_or(500);

        match method {
            LyapunovMethod::Map1D => {
                // Calculate Lyapunov exponent for logistic map
                let r = input.parameters
                    .get("r")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(3.9);

                let exponent = lyapunov::lyapunov_logistic_map(r, iterations, transient);
                let is_chaotic = exponent > 0.0;

                let mut values = std::collections::HashMap::new();
                values.insert("lyapunov_exponent".to_string(), exponent);
                values.insert("is_chaotic".to_string(), if is_chaotic { 1.0 } else { 0.0 });
                values.insert("r".to_string(), r);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "lyapunov_exponent": exponent,
                        "is_chaotic": is_chaotic,
                        "interpretation": if exponent > 0.0 {
                            "Chaotic (positive exponent indicates sensitive dependence on initial conditions)"
                        } else if exponent < -0.01 {
                            "Periodic (negative exponent indicates convergent dynamics)"
                        } else {
                            "Edge of chaos (near-zero exponent)"
                        }
                    }),
                    trajectory: None,
                    values: Some(values),
                    metadata: Some(serde_json::json!({
                        "method": "lyapunov_1d",
                        "map": "logistic",
                        "r": r,
                        "iterations": iterations,
                        "transient": transient
                    })),
                })
            }

            LyapunovMethod::Spectrum1D | LyapunovMethod::MaxLyapunov => {
                // For 1D, same as Map1D
                let r = input.parameters
                    .get("r")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(3.9);

                let exponent = lyapunov::lyapunov_logistic_map(r, iterations, transient);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "largest_exponent": exponent,
                        "is_chaotic": exponent > 0.0
                    }),
                    trajectory: None,
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("largest_exponent".to_string(), exponent);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "method": "max_lyapunov"
                    })),
                })
            }

            LyapunovMethod::Spectrum3D => {
                // 3D Lyapunov spectrum requires full implementation
                Err("3D Lyapunov spectrum calculation requires full Jacobian - use individual attractor analysis".to_string())
            }
        }
    }

    /// Handle bifurcation diagram operations
    fn handle_bifurcation(&self, bif_type: &BifurcationType, input: &ChaosInput) -> ToolResult<ChaosOutput> {
        use crate::mathematics::chaos::attractors;

        match bif_type {
            BifurcationType::LogisticMap => {
                let r_min = input.parameters.get("r_min").and_then(|v| v.as_f64()).unwrap_or(2.5);
                let r_max = input.parameters.get("r_max").and_then(|v| v.as_f64()).unwrap_or(4.0);
                let num_r = input.parameters.get("num_r").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or(100);
                let iterations = input.iterations.unwrap_or(100);
                let transient = input.parameters.get("transient").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or(500);

                let results = attractors::logistic_map_bifurcation(r_min, r_max, num_r, iterations, transient);

                // Flatten for visualization
                let mut r_values = Vec::new();
                let mut x_values = Vec::new();
                for (r, xs) in &results {
                    for x in xs {
                        r_values.push(*r);
                        x_values.push(*x);
                    }
                }

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "r_values": r_values,
                        "x_values": x_values,
                        "n_points": r_values.len()
                    }),
                    trajectory: Some(results.iter().map(|(r, xs)| {
                        let mut v = vec![*r];
                        v.extend(xs.iter().take(5)); // First 5 values per r
                        v
                    }).collect()),
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("r_min".to_string(), r_min);
                        v.insert("r_max".to_string(), r_max);
                        v.insert("n_r".to_string(), num_r as f64);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "type": "bifurcation_diagram",
                        "map": "logistic",
                        "r_range": [r_min, r_max],
                        "num_r": num_r
                    })),
                })
            }

            BifurcationType::PeriodDoubling => {
                // Period doubling cascade in logistic map
                // Find r values where period doubles
                Err("Period doubling detection not yet implemented - use bifurcation diagram".to_string())
            }

            BifurcationType::Pitchfork => {
                Err("Pitchfork bifurcation not yet implemented".to_string())
            }

            BifurcationType::SaddleNode => {
                Err("Saddle-node bifurcation not yet implemented".to_string())
            }

            BifurcationType::Hopf => {
                Err("Hopf bifurcation not yet implemented".to_string())
            }
        }
    }

    /// Handle dimension calculation operations
    fn handle_dimension(&self, method: &DimensionMethod, input: &ChaosInput) -> ToolResult<ChaosOutput> {
        use crate::mathematics::chaos::{fractals, lyapunov};

        match method {
            DimensionMethod::BoxCounting => {
                // Extract points from parameters
                let points: Vec<(f64, f64)> = input.parameters
                    .get("points")
                    .and_then(|v| v.as_array())
                    .map(|arr| {
                        arr.iter().filter_map(|p| {
                            p.as_array().and_then(|coords| {
                                if coords.len() >= 2 {
                                    Some((coords[0].as_f64()?, coords[1].as_f64()?))
                                } else {
                                    None
                                }
                            })
                        }).collect()
                    })
                    .unwrap_or_default();

                if points.is_empty() {
                    return Err("points parameter required for box counting dimension".to_string());
                }

                let min_box = input.parameters.get("min_box").and_then(|v| v.as_f64()).unwrap_or(0.01);
                let max_box = input.parameters.get("max_box").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let num_scales = input.parameters.get("num_scales").and_then(|v| v.as_u64()).map(|v| v as usize).unwrap_or(20);

                let dimension = fractals::box_counting_dimension(&points, min_box, max_box, num_scales);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "box_counting_dimension": dimension,
                        "n_points": points.len()
                    }),
                    trajectory: None,
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("dimension".to_string(), dimension);
                        v.insert("n_points".to_string(), points.len() as f64);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "method": "box_counting",
                        "min_box": min_box,
                        "max_box": max_box,
                        "num_scales": num_scales
                    })),
                })
            }

            DimensionMethod::KaplanYorke => {
                // Kaplan-Yorke dimension from Lyapunov exponents
                let exponents: Vec<f64> = input.parameters
                    .get("exponents")
                    .and_then(|v| v.as_array())
                    .map(|arr| arr.iter().filter_map(|v| v.as_f64()).collect())
                    .unwrap_or_default();

                if exponents.is_empty() {
                    return Err("exponents parameter required for Kaplan-Yorke dimension".to_string());
                }

                let dimension = lyapunov::kaplan_yorke_dimension(&exponents);

                Ok(ChaosOutput {
                    result: serde_json::json!({
                        "kaplan_yorke_dimension": dimension,
                        "exponents": exponents
                    }),
                    trajectory: None,
                    values: Some({
                        let mut v = std::collections::HashMap::new();
                        v.insert("dimension".to_string(), dimension);
                        v
                    }),
                    metadata: Some(serde_json::json!({
                        "method": "kaplan_yorke"
                    })),
                })
            }

            DimensionMethod::CorrelationDimension => {
                Err("Correlation dimension not yet implemented".to_string())
            }

            DimensionMethod::HausdorffDimension => {
                Err("Hausdorff dimension not yet implemented".to_string())
            }
        }
    }
}

impl Default for UnifiedChaos {
    fn default() -> Self {
        Self::new()
    }
}

impl Chaos for UnifiedChaos {
    fn chaos(&self, input: &ChaosInput) -> ToolResult<ChaosOutput> {
        match &input.operation {
            ChaosOp::Fractal(fractal_type) => self.handle_fractal(fractal_type, input),
            ChaosOp::Attractor(attractor_type) => self.handle_attractor(attractor_type, input),
            ChaosOp::Lyapunov(method) => self.handle_lyapunov(method, input),
            ChaosOp::Bifurcation(bif_type) => self.handle_bifurcation(bif_type, input),
            ChaosOp::Dimension(method) => self.handle_dimension(method, input),
        }
    }
}
