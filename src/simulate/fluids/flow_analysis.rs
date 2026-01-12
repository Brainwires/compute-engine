use serde::Serialize;
use serde_json::{Value, json};
use std::error::Error;

pub struct FlowAnalyzer;

impl FlowAnalyzer {
    pub fn new() -> Self {
        Self
    }

    pub fn analyze(
        &self,
        velocity_field: Value,
        grid_spacing: Value,
        analysis_types: Vec<String>,
    ) -> Result<Value, Box<dyn Error>> {
        // Parse velocity field
        let u_velocity = velocity_field["u_velocity"]
            .as_array()
            .ok_or("Missing u_velocity field")?;
        let v_velocity = velocity_field["v_velocity"]
            .as_array()
            .ok_or("Missing v_velocity field")?;

        // Parse grid spacing
        let dx = grid_spacing["dx"]
            .as_f64()
            .ok_or("Missing dx in grid_spacing")?;
        let dy = grid_spacing["dy"]
            .as_f64()
            .ok_or("Missing dy in grid_spacing")?;

        // Convert JSON arrays to 2D vectors
        let u: Result<Vec<Vec<f64>>, _> = u_velocity
            .iter()
            .map(|row| {
                row.as_array()
                    .ok_or("Invalid u_velocity row format")?
                    .iter()
                    .map(|val| val.as_f64().ok_or("Invalid u_velocity value"))
                    .collect()
            })
            .collect();
        let u = u?;

        let v: Result<Vec<Vec<f64>>, _> = v_velocity
            .iter()
            .map(|row| {
                row.as_array()
                    .ok_or("Invalid v_velocity row format")?
                    .iter()
                    .map(|val| val.as_f64().ok_or("Invalid v_velocity value"))
                    .collect()
            })
            .collect();
        let v = v?;

        let rows = u.len();
        let cols = u[0].len();

        // Validate dimensions
        if v.len() != rows || v[0].len() != cols {
            return Err("Velocity field dimensions don't match".into());
        }

        let mut results = serde_json::Map::new();

        // Perform requested analyses
        for analysis_type in &analysis_types {
            match analysis_type.as_str() {
                "velocity_magnitude" => {
                    let magnitude = self.calculate_velocity_magnitude(&u, &v);
                    results.insert("velocity_magnitude".to_string(), json!(magnitude));
                }
                "vorticity" => {
                    let vorticity = self.calculate_vorticity(&u, &v, dx, dy)?;
                    results.insert("vorticity".to_string(), json!(vorticity));
                }
                "divergence" => {
                    let divergence = self.calculate_divergence(&u, &v, dx, dy)?;
                    results.insert("divergence".to_string(), json!(divergence));
                }
                "stream_function" => {
                    let stream_function = self.calculate_stream_function(&u, &v, dx, dy)?;
                    results.insert("stream_function".to_string(), json!(stream_function));
                }
                "pressure_field" => {
                    // For pressure field calculation, we would need additional information
                    // This is a placeholder implementation
                    results.insert("pressure_field".to_string(),
                                 json!({"error": "Pressure field calculation requires additional boundary conditions and fluid properties"}));
                }
                "kinetic_energy" => {
                    let kinetic_energy = self.calculate_kinetic_energy(&u, &v);
                    results.insert(
                        "kinetic_energy".to_string(),
                        json!({
                            "total": kinetic_energy.total,
                            "field": kinetic_energy.field
                        }),
                    );
                }
                "strain_rate" => {
                    let strain_rate = self.calculate_strain_rate(&u, &v, dx, dy)?;
                    results.insert("strain_rate".to_string(), json!(strain_rate));
                }
                _ => {
                    results.insert(
                        analysis_type.clone(),
                        json!(format!("Analysis type '{}' not implemented", analysis_type)),
                    );
                }
            }
        }

        Ok(json!({
            "analysis_results": results,
            "grid_info": {
                "dimensions": [rows, cols],
                "grid_spacing": {"dx": dx, "dy": dy}
            },
            "analysis_types_requested": analysis_types
        }))
    }

    fn calculate_velocity_magnitude(&self, u: &Vec<Vec<f64>>, v: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let rows = u.len();
        let cols = u[0].len();
        let mut magnitude = vec![vec![0.0; cols]; rows];

        for i in 0..rows {
            for j in 0..cols {
                magnitude[i][j] = (u[i][j].powi(2) + v[i][j].powi(2)).sqrt();
            }
        }

        magnitude
    }

    fn calculate_vorticity(
        &self,
        u: &Vec<Vec<f64>>,
        v: &Vec<Vec<f64>>,
        dx: f64,
        dy: f64,
    ) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
        let rows = u.len();
        let cols = u[0].len();

        if rows < 3 || cols < 3 {
            return Err("Grid too small for vorticity calculation (need at least 3x3)".into());
        }

        let mut vorticity = vec![vec![0.0; cols]; rows];

        // Calculate vorticity using central differences (interior points only)
        for i in 1..rows - 1 {
            for j in 1..cols - 1 {
                let dvdx = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dx);
                let dudy = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dy);
                vorticity[i][j] = dvdx - dudy;
            }
        }

        Ok(vorticity)
    }

    fn calculate_divergence(
        &self,
        u: &Vec<Vec<f64>>,
        v: &Vec<Vec<f64>>,
        dx: f64,
        dy: f64,
    ) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
        let rows = u.len();
        let cols = u[0].len();

        if rows < 3 || cols < 3 {
            return Err("Grid too small for divergence calculation (need at least 3x3)".into());
        }

        let mut divergence = vec![vec![0.0; cols]; rows];

        // Calculate divergence using central differences (interior points only)
        for i in 1..rows - 1 {
            for j in 1..cols - 1 {
                let dudx = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dx);
                let dvdy = (v[i + 1][j] - v[i - 1][j]) / (2.0 * dy);
                divergence[i][j] = dudx + dvdy;
            }
        }

        Ok(divergence)
    }

    fn calculate_stream_function(
        &self,
        u: &Vec<Vec<f64>>,
        v: &Vec<Vec<f64>>,
        dx: f64,
        dy: f64,
    ) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
        let rows = u.len();
        let cols = u[0].len();
        let mut psi = vec![vec![0.0; cols]; rows];

        // Calculate stream function using line integration
        // ∂ψ/∂y = u, ∂ψ/∂x = -v

        // Start from bottom-left corner (psi = 0)
        // Integrate along y-direction first
        for i in 1..rows {
            for j in 0..cols {
                psi[i][j] = psi[i - 1][j] + u[i][j] * dy;
            }
        }

        // Check consistency with v-component (optional validation)
        let mut max_inconsistency = 0.0_f64;
        for i in 0..rows {
            for j in 1..cols - 1 {
                let dpsidx = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * dx);
                let inconsistency = (dpsidx + v[i][j]).abs();
                max_inconsistency = max_inconsistency.max(inconsistency);
            }
        }

        if max_inconsistency > 0.1 {
            eprintln!(
                "Warning: Stream function calculation has high inconsistency: {}",
                max_inconsistency
            );
        }

        Ok(psi)
    }

    fn calculate_kinetic_energy(
        &self,
        u: &Vec<Vec<f64>>,
        v: &Vec<Vec<f64>>,
    ) -> KineticEnergyResult {
        let rows = u.len();
        let cols = u[0].len();
        let mut field = vec![vec![0.0; cols]; rows];
        let mut total = 0.0;

        for i in 0..rows {
            for j in 0..cols {
                let ke = 0.5 * (u[i][j].powi(2) + v[i][j].powi(2));
                field[i][j] = ke;
                total += ke;
            }
        }

        // Normalize total by number of grid points
        total /= (rows * cols) as f64;

        KineticEnergyResult { total, field }
    }

    fn calculate_strain_rate(
        &self,
        u: &Vec<Vec<f64>>,
        v: &Vec<Vec<f64>>,
        dx: f64,
        dy: f64,
    ) -> Result<StrainRateResult, Box<dyn Error>> {
        let rows = u.len();
        let cols = u[0].len();

        if rows < 3 || cols < 3 {
            return Err("Grid too small for strain rate calculation (need at least 3x3)".into());
        }

        let mut s11 = vec![vec![0.0; cols]; rows]; // ∂u/∂x
        let mut s22 = vec![vec![0.0; cols]; rows]; // ∂v/∂y
        let mut s12 = vec![vec![0.0; cols]; rows]; // 0.5 * (∂u/∂y + ∂v/∂x)
        let mut magnitude = vec![vec![0.0; cols]; rows];

        // Calculate strain rate components (interior points only)
        for i in 1..rows - 1 {
            for j in 1..cols - 1 {
                let dudx = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dx);
                let dvdy = (v[i + 1][j] - v[i - 1][j]) / (2.0 * dy);
                let dudy = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dy);
                let dvdx = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dx);

                s11[i][j] = dudx;
                s22[i][j] = dvdy;
                s12[i][j] = 0.5 * (dudy + dvdx);

                // Magnitude of strain rate tensor: |S| = sqrt(2 * S_ij * S_ij)
                magnitude[i][j] = (2.0
                    * (s11[i][j].powi(2) + s22[i][j].powi(2) + 2.0 * s12[i][j].powi(2)))
                .sqrt();
            }
        }

        Ok(StrainRateResult {
            s11,
            s22,
            s12,
            magnitude,
        })
    }
}

struct KineticEnergyResult {
    total: f64,
    field: Vec<Vec<f64>>,
}

#[derive(Serialize)]
struct StrainRateResult {
    s11: Vec<Vec<f64>>,
    s22: Vec<Vec<f64>>,
    s12: Vec<Vec<f64>>,
    magnitude: Vec<Vec<f64>>,
}
