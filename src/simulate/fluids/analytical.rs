use serde_json::{Value, json};
use std::error::Error;
use std::f64::consts::PI;

pub struct AnalyticalSolutions;

impl AnalyticalSolutions {
    pub fn new() -> Self {
        Self
    }

    pub fn solve(&self, solution_type: &str, parameters: Value) -> Result<Value, Box<dyn Error>> {
        match solution_type {
            "couette_flow" => self.couette_flow(parameters),
            "poiseuille_flow" => self.poiseuille_flow(parameters),
            "stagnation_point_flow" => self.stagnation_point_flow(parameters),
            "rotating_flow" => self.rotating_flow(parameters),
            "potential_vortex" => self.potential_vortex(parameters),
            "rankine_vortex" => self.rankine_vortex(parameters),
            _ => Err(format!("Unknown analytical solution type: {}", solution_type).into()),
        }
    }

    fn couette_flow(&self, params: Value) -> Result<Value, Box<dyn Error>> {
        let plate_velocity = params["plate_velocity"]
            .as_f64()
            .ok_or("Missing plate_velocity parameter")?;
        let plate_separation = params["plate_separation"]
            .as_f64()
            .ok_or("Missing plate_separation parameter")?;
        let num_points = params["num_points"].as_u64().unwrap_or(100) as usize;

        // Generate y coordinates
        let mut y_coords = Vec::with_capacity(num_points);
        let mut velocities = Vec::with_capacity(num_points);
        let mut shear_stress = Vec::with_capacity(num_points);

        let viscosity = params["viscosity"].as_f64().unwrap_or(1e-6);
        let shear_rate = plate_velocity / plate_separation;

        for i in 0..num_points {
            let y = (i as f64 / (num_points - 1) as f64) * plate_separation;
            let u = (plate_velocity / plate_separation) * y;
            let tau = viscosity * shear_rate; // constant shear stress

            y_coords.push(y);
            velocities.push(u);
            shear_stress.push(tau);
        }

        Ok(json!({
            "solution_type": "couette_flow",
            "description": "Flow between two parallel plates with the top plate moving",
            "analytical_solution": {
                "velocity_profile": format!("u(y) = ({}/{})*y", plate_velocity, plate_separation),
                "velocity_gradient": format!("du/dy = {}", shear_rate),
                "shear_stress": format!("τ = μ * {} = {}", shear_rate, viscosity * shear_rate)
            },
            "numerical_data": {
                "y_coordinates": y_coords,
                "u_velocity": velocities,
                "shear_stress": shear_stress
            },
            "parameters": {
                "plate_velocity": plate_velocity,
                "plate_separation": plate_separation,
                "viscosity": viscosity,
                "shear_rate": shear_rate
            }
        }))
    }

    fn poiseuille_flow(&self, params: Value) -> Result<Value, Box<dyn Error>> {
        let pressure_gradient = params["pressure_gradient"]
            .as_f64()
            .ok_or("Missing pressure_gradient parameter")?;
        let channel_height = params["channel_height"]
            .as_f64()
            .ok_or("Missing channel_height parameter")?;
        let viscosity = params["viscosity"]
            .as_f64()
            .ok_or("Missing viscosity parameter")?;
        let num_points = params["num_points"].as_u64().unwrap_or(100) as usize;

        // Calculate key parameters
        let max_velocity = -pressure_gradient * channel_height.powi(2) / (8.0 * viscosity);
        let flow_rate = (2.0 * max_velocity * channel_height) / 3.0;
        let wall_shear_stress = -pressure_gradient * channel_height / 2.0;

        // Generate profile data
        let mut y_coords = Vec::with_capacity(num_points);
        let mut velocities = Vec::with_capacity(num_points);
        let mut shear_stress = Vec::with_capacity(num_points);

        for i in 0..num_points {
            let y = (i as f64 / (num_points - 1) as f64) * channel_height;
            let u = (-pressure_gradient / (2.0 * viscosity)) * y * (channel_height - y);
            let tau = -pressure_gradient * (channel_height / 2.0 - y);

            y_coords.push(y);
            velocities.push(u);
            shear_stress.push(tau);
        }

        Ok(json!({
            "solution_type": "poiseuille_flow",
            "description": "Pressure-driven flow between parallel plates",
            "analytical_solution": {
                "velocity_profile": format!("u(y) = (-dp/dx)/(2μ) * y * (H - y)"),
                "maximum_velocity": max_velocity,
                "flow_rate_per_width": flow_rate,
                "wall_shear_stress": wall_shear_stress
            },
            "numerical_data": {
                "y_coordinates": y_coords,
                "u_velocity": velocities,
                "shear_stress": shear_stress
            },
            "parameters": {
                "pressure_gradient": pressure_gradient,
                "channel_height": channel_height,
                "viscosity": viscosity
            },
            "derived_quantities": {
                "reynolds_number": flow_rate * channel_height / viscosity,
                "friction_factor": wall_shear_stress / (0.5 * flow_rate.powi(2))
            }
        }))
    }

    fn stagnation_point_flow(&self, params: Value) -> Result<Value, Box<dyn Error>> {
        let strain_rate = params["strain_rate"]
            .as_f64()
            .ok_or("Missing strain_rate parameter")?;
        let domain_size = params["domain_size"].as_f64().unwrap_or(2.0);
        let num_points = params["num_points"].as_u64().unwrap_or(50) as usize;

        // Generate 2D velocity field
        let mut x_coords = Vec::new();
        let mut y_coords = Vec::new();
        let mut u_field = Vec::new();
        let mut v_field = Vec::new();
        let mut pressure_field = Vec::new();

        let density = params["density"].as_f64().unwrap_or(1.0);

        for i in 0..num_points {
            let mut u_row = Vec::new();
            let mut v_row = Vec::new();
            let mut p_row = Vec::new();

            let y = (i as f64 / (num_points - 1) as f64 - 0.5) * domain_size;
            y_coords.push(y);

            if i == 0 {
                for j in 0..num_points {
                    let x = (j as f64 / (num_points - 1) as f64 - 0.5) * domain_size;
                    x_coords.push(x);
                }
            }

            for j in 0..num_points {
                let x = x_coords[j];

                // Stagnation point flow: u = ax, v = -ay
                let u = strain_rate * x;
                let v = -strain_rate * y;

                // Pressure field: p = p0 - (ρa²/2)(x² + y²)
                let p = -0.5 * density * strain_rate.powi(2) * (x.powi(2) + y.powi(2));

                u_row.push(u);
                v_row.push(v);
                p_row.push(p);
            }

            u_field.push(u_row);
            v_field.push(v_row);
            pressure_field.push(p_row);
        }

        Ok(json!({
            "solution_type": "stagnation_point_flow",
            "description": "Two-dimensional stagnation point flow",
            "analytical_solution": {
                "velocity_u": format!("u(x,y) = {} * x", strain_rate),
                "velocity_v": format!("v(x,y) = {} * y", -strain_rate),
                "pressure": format!("p(x,y) = -ρa²/2 * (x² + y²)")
            },
            "numerical_data": {
                "x_coordinates": x_coords,
                "y_coordinates": y_coords,
                "u_velocity": u_field,
                "v_velocity": v_field,
                "pressure": pressure_field
            },
            "parameters": {
                "strain_rate": strain_rate,
                "domain_size": domain_size,
                "density": density
            },
            "flow_characteristics": {
                "stagnation_points": [{"x": 0.0, "y": 0.0}],
                "flow_type": "irrotational",
                "vorticity": 0.0
            }
        }))
    }

    fn rotating_flow(&self, params: Value) -> Result<Value, Box<dyn Error>> {
        let angular_velocity = params["angular_velocity"]
            .as_f64()
            .ok_or("Missing angular_velocity parameter")?;
        let radius_max = params["radius_max"].as_f64().unwrap_or(1.0);
        let num_points = params["num_points"].as_u64().unwrap_or(50) as usize;

        // Generate cylindrical coordinates
        let mut r_coords = Vec::new();
        let mut theta_coords = Vec::new();
        let mut velocity_theta = Vec::new();

        // Cartesian field for visualization
        let mut x_coords = Vec::new();
        let mut y_coords = Vec::new();
        let mut u_field = Vec::new();
        let mut v_field = Vec::new();

        // Generate radial data
        for i in 0..num_points {
            let r = (i as f64 / (num_points - 1) as f64) * radius_max;
            let v_theta = angular_velocity * r; // solid body rotation

            r_coords.push(r);
            velocity_theta.push(v_theta);
        }

        // Generate 2D Cartesian field
        for i in 0..num_points {
            let mut u_row = Vec::new();
            let mut v_row = Vec::new();

            let y = (i as f64 / (num_points - 1) as f64 - 0.5) * 2.0 * radius_max;
            y_coords.push(y);

            if i == 0 {
                for j in 0..num_points {
                    let x = (j as f64 / (num_points - 1) as f64 - 0.5) * 2.0 * radius_max;
                    x_coords.push(x);
                }
            }

            for j in 0..num_points {
                let x = x_coords[j];

                // Convert to cylindrical and back for rotating flow
                let r = (x.powi(2) + y.powi(2)).sqrt();

                if r > 0.0 {
                    // Solid body rotation: v_θ = ωr
                    let v_theta = angular_velocity * r;

                    // Convert to Cartesian: u = -v_θ sin(θ), v = v_θ cos(θ)
                    let sin_theta = y / r;
                    let cos_theta = x / r;

                    let u = -v_theta * sin_theta;
                    let v = v_theta * cos_theta;

                    u_row.push(u);
                    v_row.push(v);
                } else {
                    u_row.push(0.0);
                    v_row.push(0.0);
                }
            }

            u_field.push(u_row);
            v_field.push(v_row);
        }

        // Generate angular data
        for i in 0..36 {
            // 10-degree increments
            let theta = (i as f64 * 10.0) * PI / 180.0;
            theta_coords.push(theta);
        }

        Ok(json!({
            "solution_type": "rotating_flow",
            "description": "Solid body rotation flow",
            "analytical_solution": {
                "velocity_theta": format!("v_θ(r) = {} * r", angular_velocity),
                "velocity_radial": "v_r = 0",
                "vorticity": format!("ω = 2 * {} = {}", angular_velocity, 2.0 * angular_velocity)
            },
            "cylindrical_data": {
                "r_coordinates": r_coords,
                "theta_coordinates": theta_coords,
                "velocity_theta": velocity_theta
            },
            "cartesian_data": {
                "x_coordinates": x_coords,
                "y_coordinates": y_coords,
                "u_velocity": u_field,
                "v_velocity": v_field
            },
            "parameters": {
                "angular_velocity": angular_velocity,
                "radius_max": radius_max
            },
            "flow_characteristics": {
                "flow_type": "rotational",
                "vorticity": 2.0 * angular_velocity,
                "circulation": 2.0 * PI * angular_velocity * radius_max.powi(2)
            }
        }))
    }

    fn potential_vortex(&self, params: Value) -> Result<Value, Box<dyn Error>> {
        let circulation = params["circulation"]
            .as_f64()
            .ok_or("Missing circulation parameter")?;
        let core_radius = params["core_radius"].as_f64().unwrap_or(0.1);
        let radius_max = params["radius_max"].as_f64().unwrap_or(2.0);
        let num_points = params["num_points"].as_u64().unwrap_or(50) as usize;

        // Generate radial velocity profile
        let mut r_coords = Vec::new();
        let mut velocity_theta = Vec::new();

        for i in 1..num_points {
            // Start from i=1 to avoid r=0
            let r = core_radius + (i as f64 / (num_points - 1) as f64) * (radius_max - core_radius);
            let v_theta = circulation / (2.0 * PI * r);

            r_coords.push(r);
            velocity_theta.push(v_theta);
        }

        Ok(json!({
            "solution_type": "potential_vortex",
            "description": "Irrotational vortex flow",
            "analytical_solution": {
                "velocity_theta": format!("v_θ(r) = Γ/(2πr)"),
                "velocity_radial": "v_r = 0",
                "circulation": circulation,
                "vorticity": "ω = 0 (except at origin)"
            },
            "numerical_data": {
                "r_coordinates": r_coords,
                "velocity_theta": velocity_theta
            },
            "parameters": {
                "circulation": circulation,
                "core_radius": core_radius,
                "radius_max": radius_max
            },
            "flow_characteristics": {
                "flow_type": "irrotational",
                "singularity": "Line vortex at r = 0"
            }
        }))
    }

    fn rankine_vortex(&self, params: Value) -> Result<Value, Box<dyn Error>> {
        let circulation = params["circulation"]
            .as_f64()
            .ok_or("Missing circulation parameter")?;
        let core_radius = params["core_radius"]
            .as_f64()
            .ok_or("Missing core_radius parameter")?;
        let radius_max = params["radius_max"].as_f64().unwrap_or(3.0);
        let num_points = params["num_points"].as_u64().unwrap_or(100) as usize;

        let mut r_coords = Vec::new();
        let mut velocity_theta = Vec::new();
        let mut vorticity = Vec::new();

        // Core vorticity
        let omega_core = circulation / (PI * core_radius.powi(2));

        for i in 0..num_points {
            let r = (i as f64 / (num_points - 1) as f64) * radius_max;

            let (v_theta, omega) = if r <= core_radius {
                // Inside core: solid body rotation
                let v_theta = omega_core * r / 2.0;
                (v_theta, omega_core)
            } else {
                // Outside core: potential vortex
                let v_theta = circulation / (2.0 * PI * r);
                (v_theta, 0.0)
            };

            r_coords.push(r);
            velocity_theta.push(v_theta);
            vorticity.push(omega);
        }

        Ok(json!({
            "solution_type": "rankine_vortex",
            "description": "Combined vortex with solid core and potential flow outside",
            "analytical_solution": {
                "inner_region": format!("v_θ(r) = (Γ/2πR²) * r, r ≤ {}", core_radius),
                "outer_region": format!("v_θ(r) = Γ/(2πr), r > {}", core_radius),
                "core_vorticity": omega_core
            },
            "numerical_data": {
                "r_coordinates": r_coords,
                "velocity_theta": velocity_theta,
                "vorticity": vorticity
            },
            "parameters": {
                "circulation": circulation,
                "core_radius": core_radius,
                "radius_max": radius_max,
                "core_vorticity": omega_core
            },
            "flow_characteristics": {
                "transition_radius": core_radius,
                "maximum_velocity": circulation / (2.0 * PI * core_radius)
            }
        }))
    }
}
