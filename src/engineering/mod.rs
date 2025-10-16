/**
 * Engineering Module
 *
 * Implements practical engineering formulas:
 * - Acoustics (SPL, Doppler, impedance, reverberation)
 * - Materials Science (stress-strain, fracture mechanics)
 * - Fluid Mechanics Advanced (Bernoulli, Poiseuille, drag)
 * - Control Theory (PID tuning basics)
 */
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineeringInput {
    pub discipline: EngineeringDiscipline,
    pub parameters: EngineeringParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineeringDiscipline {
    Acoustics,
    Materials,
    FluidMechanics,
    ControlTheory,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineeringParams {
    // Acoustics
    pub pressure_rms: Option<f64>,           // Pa
    pub reference_pressure: Option<f64>,     // Pa (default 20 μPa)
    pub velocity_source: Option<f64>,        // m/s
    pub velocity_observer: Option<f64>,      // m/s
    pub sound_speed: Option<f64>,            // m/s
    pub frequency: Option<f64>,              // Hz
    pub room_volume: Option<f64>,            // m³
    pub absorption_coefficient: Option<f64>, // 0-1

    // Materials
    pub stress: Option<f64>,          // Pa
    pub strain: Option<f64>,          // dimensionless
    pub youngs_modulus: Option<f64>,  // Pa
    pub yield_strength: Option<f64>,  // Pa
    pub crack_length: Option<f64>,    // m
    pub geometry_factor: Option<f64>, // Y (dimensionless)

    // Fluid mechanics
    pub velocity: Option<f64>,             // m/s
    pub pressure_1: Option<f64>,           // Pa
    pub pressure_2: Option<f64>,           // Pa
    pub height_1: Option<f64>,             // m
    pub height_2: Option<f64>,             // m
    pub density: Option<f64>,              // kg/m³
    pub viscosity: Option<f64>,            // Pa·s
    pub radius: Option<f64>,               // m
    pub length: Option<f64>,               // m
    pub flow_rate: Option<f64>,            // m³/s
    pub drag_coefficient: Option<f64>,     // Cd
    pub cross_sectional_area: Option<f64>, // m²

    // Control theory
    pub kp: Option<f64>, // Proportional gain
    pub ki: Option<f64>, // Integral gain
    pub kd: Option<f64>, // Derivative gain
    pub setpoint: Option<f64>,
    pub process_variable: Option<f64>,
    pub time_constant: Option<f64>, // τ
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineeringResult {
    pub value: f64,
    pub unit: String,
    pub formula_used: String,
    pub classification: Option<String>,
    pub interpretation: String,
    pub additional: Option<std::collections::HashMap<String, f64>>,
}

const REFERENCE_PRESSURE_ACOUSTICS: f64 = 20e-6; // 20 μPa
const G_ACCEL: f64 = 9.80665; // m/s²

pub fn calculate_engineering(input: EngineeringInput) -> Result<EngineeringResult, String> {
    match input.discipline {
        EngineeringDiscipline::Acoustics => calculate_acoustics(&input.parameters),
        EngineeringDiscipline::Materials => calculate_materials(&input.parameters),
        EngineeringDiscipline::FluidMechanics => calculate_fluid_mechanics(&input.parameters),
        EngineeringDiscipline::ControlTheory => calculate_control(&input.parameters),
    }
}

/// Sound Pressure Level and acoustic calculations
fn calculate_acoustics(params: &EngineeringParams) -> Result<EngineeringResult, String> {
    let mut additional = std::collections::HashMap::new();

    // SPL calculation: SPL = 20·log₁₀(P/P₀)
    if let Some(p_rms) = params.pressure_rms {
        let p_ref = params
            .reference_pressure
            .unwrap_or(REFERENCE_PRESSURE_ACOUSTICS);
        let spl = 20.0 * (p_rms / p_ref).log10();

        additional.insert("pressure_pa".to_string(), p_rms);

        let classification = match spl {
            s if s < 30.0 => "Whisper quiet",
            s if s < 50.0 => "Quiet",
            s if s < 70.0 => "Moderate",
            s if s < 85.0 => "Loud",
            s if s < 100.0 => "Very loud - hearing damage risk",
            s if s < 120.0 => "Dangerous - immediate damage risk",
            _ => "Pain threshold - severe damage",
        };

        return Ok(EngineeringResult {
            value: spl,
            unit: "dB SPL".to_string(),
            formula_used: "SPL = 20·log₁₀(P/P₀)".to_string(),
            classification: Some(classification.to_string()),
            interpretation: format!("{:.1} dB SPL - {}", spl, classification),
            additional: Some(additional),
        });
    }

    // Doppler effect: f' = f·(v + v_observer)/(v - v_source)
    if let (Some(f), Some(v_sound)) = (params.frequency, params.sound_speed) {
        let v_obs = params.velocity_observer.unwrap_or(0.0);
        let v_src = params.velocity_source.unwrap_or(0.0);

        let f_observed = f * (v_sound + v_obs) / (v_sound - v_src);
        let shift = f_observed - f;
        let shift_percent = (shift / f) * 100.0;

        additional.insert("original_frequency".to_string(), f);
        additional.insert("frequency_shift".to_string(), shift);
        additional.insert("shift_percent".to_string(), shift_percent);

        return Ok(EngineeringResult {
            value: f_observed,
            unit: "Hz".to_string(),
            formula_used: "Doppler: f' = f·(v + v_obs)/(v - v_src)".to_string(),
            classification: None,
            interpretation: format!("{:.1} Hz ({:+.1}% shift)", f_observed, shift_percent),
            additional: Some(additional),
        });
    }

    // Sabine reverberation time: RT₆₀ = 0.161·V/A
    if let (Some(volume), Some(alpha)) = (params.room_volume, params.absorption_coefficient) {
        // Approximate total absorption assuming uniform coefficient
        let total_absorption = alpha * volume.powf(2.0 / 3.0) * 6.0; // Rough estimate of surface area

        let rt60 = 0.161 * volume / total_absorption;

        additional.insert("volume_m3".to_string(), volume);
        additional.insert("absorption_coefficient".to_string(), alpha);

        let room_type = if rt60 < 0.5 {
            "Very dry (recording studio)"
        } else if rt60 < 1.0 {
            "Dry (lecture hall)"
        } else if rt60 < 2.0 {
            "Medium (auditorium)"
        } else {
            "Reverberant (cathedral)"
        };

        return Ok(EngineeringResult {
            value: rt60,
            unit: "seconds".to_string(),
            formula_used: "Sabine: RT₆₀ = 0.161·V/A".to_string(),
            classification: Some(room_type.to_string()),
            interpretation: format!("{:.2}s reverberation - {}", rt60, room_type),
            additional: Some(additional),
        });
    }

    Err("Insufficient acoustic parameters".to_string())
}

/// Materials science: stress-strain, fracture mechanics
fn calculate_materials(params: &EngineeringParams) -> Result<EngineeringResult, String> {
    let mut additional = std::collections::HashMap::new();

    // Hooke's law: σ = E·ε
    if let Some(e_mod) = params.youngs_modulus {
        if let Some(strain) = params.strain {
            let stress = e_mod * strain;

            let yield_str = params.yield_strength.unwrap_or(e_mod * 0.002); // Approximate
            let safety_factor = yield_str / stress;

            additional.insert("stress_pa".to_string(), stress);
            additional.insert("safety_factor".to_string(), safety_factor);

            let regime = if stress < yield_str * 0.5 {
                "Elastic (safe)"
            } else if stress < yield_str {
                "Elastic (approaching yield)"
            } else {
                "Plastic deformation!"
            };

            return Ok(EngineeringResult {
                value: stress / 1e6,
                unit: "MPa".to_string(),
                formula_used: "Hooke's law: σ = E·ε".to_string(),
                classification: Some(regime.to_string()),
                interpretation: format!("{:.1} MPa - {}", stress / 1e6, regime),
                additional: Some(additional),
            });
        }

        if let Some(stress) = params.stress {
            let strain = stress / e_mod;

            return Ok(EngineeringResult {
                value: strain * 100.0,
                unit: "% strain".to_string(),
                formula_used: "Hooke's law: ε = σ/E".to_string(),
                classification: None,
                interpretation: format!("{:.3}% strain", strain * 100.0),
                additional: None,
            });
        }
    }

    // Fracture mechanics: K = Y·σ·√(π·a)
    if let (Some(stress), Some(crack_len), Some(y)) =
        (params.stress, params.crack_length, params.geometry_factor)
    {
        let k_factor = y * stress * (PI * crack_len).sqrt();

        // Typical fracture toughness values (MPa·√m)
        let k_ic_steel = 50e6; // ~50 MPa·√m in Pa·√m
        let safety_factor = k_ic_steel / k_factor;

        additional.insert("stress_intensity_factor".to_string(), k_factor / 1e6);
        additional.insert("safety_factor_vs_steel".to_string(), safety_factor);

        let criticality = if k_factor < k_ic_steel * 0.3 {
            "Safe"
        } else if k_factor < k_ic_steel * 0.7 {
            "Monitor crack growth"
        } else {
            "Critical - failure risk!"
        };

        return Ok(EngineeringResult {
            value: k_factor / 1e6,
            unit: "MPa·√m".to_string(),
            formula_used: "Fracture: K = Y·σ·√(π·a)".to_string(),
            classification: Some(criticality.to_string()),
            interpretation: format!("{:.1} MPa·√m - {}", k_factor / 1e6, criticality),
            additional: Some(additional),
        });
    }

    Err("Insufficient materials parameters".to_string())
}

/// Advanced fluid mechanics
fn calculate_fluid_mechanics(params: &EngineeringParams) -> Result<EngineeringResult, String> {
    let mut additional = std::collections::HashMap::new();

    // Bernoulli's equation: P₁ + ½ρv₁² + ρgh₁ = P₂ + ½ρv₂² + ρgh₂
    if let (Some(rho), Some(p1), Some(v1), Some(h1)) = (
        params.density,
        params.pressure_1,
        params.velocity,
        params.height_1,
    ) {
        let h2 = params.height_2.unwrap_or(h1);
        let total_pressure_1 = p1 + 0.5 * rho * v1 * v1 + rho * G_ACCEL * h1;

        // If P2 given, calculate v2
        if let Some(p2) = params.pressure_2 {
            let dynamic_p2 = total_pressure_1 - p2 - rho * G_ACCEL * h2;
            let v2 = (2.0 * dynamic_p2 / rho).max(0.0).sqrt();

            additional.insert("velocity_2".to_string(), v2);
            additional.insert("velocity_change".to_string(), v2 - v1);

            return Ok(EngineeringResult {
                value: v2,
                unit: "m/s".to_string(),
                formula_used: "Bernoulli: P + ½ρv² + ρgh = const".to_string(),
                classification: None,
                interpretation: format!("Exit velocity: {:.2} m/s", v2),
                additional: Some(additional),
            });
        }
    }

    // Poiseuille's law: Q = (π·ΔP·r⁴)/(8·η·L)
    if let (Some(delta_p), Some(r), Some(length), Some(eta)) = (
        params.pressure_1,
        params.radius,
        params.length,
        params.viscosity,
    ) {
        let flow_rate = (PI * delta_p * r.powi(4)) / (8.0 * eta * length);

        // Calculate Reynolds number
        let v_avg = flow_rate / (PI * r * r);
        let re = (params.density.unwrap_or(1000.0) * v_avg * 2.0 * r) / eta;

        additional.insert("flow_rate_m3_per_s".to_string(), flow_rate);
        additional.insert("average_velocity".to_string(), v_avg);
        additional.insert("reynolds_number".to_string(), re);

        let flow_regime = if re < 2300.0 {
            "Laminar (Poiseuille valid)"
        } else {
            "Turbulent (Poiseuille invalid!)"
        };

        return Ok(EngineeringResult {
            value: flow_rate * 1000.0,
            unit: "L/s".to_string(),
            formula_used: "Poiseuille: Q = (π·ΔP·r⁴)/(8·η·L)".to_string(),
            classification: Some(flow_regime.to_string()),
            interpretation: format!("{:.3} L/s - {}", flow_rate * 1000.0, flow_regime),
            additional: Some(additional),
        });
    }

    // Drag force: F_D = ½·ρ·v²·C_D·A
    if let (Some(rho), Some(v), Some(cd), Some(area)) = (
        params.density,
        params.velocity,
        params.drag_coefficient,
        params.cross_sectional_area,
    ) {
        let drag_force = 0.5 * rho * v * v * cd * area;

        // Calculate terminal velocity if needed
        additional.insert("drag_force_n".to_string(), drag_force);
        additional.insert("drag_power_w".to_string(), drag_force * v);

        return Ok(EngineeringResult {
            value: drag_force,
            unit: "N".to_string(),
            formula_used: "Drag: F_D = ½·ρ·v²·C_D·A".to_string(),
            classification: None,
            interpretation: format!("{:.1} N drag force at {:.1} m/s", drag_force, v),
            additional: Some(additional),
        });
    }

    Err("Insufficient fluid mechanics parameters".to_string())
}

/// Basic PID control calculations
fn calculate_control(params: &EngineeringParams) -> Result<EngineeringResult, String> {
    let setpoint = params.setpoint.ok_or("Setpoint required")?;
    let pv = params.process_variable.ok_or("Process variable required")?;

    let error = setpoint - pv;

    let mut additional = std::collections::HashMap::new();

    // Basic PID output calculation (simplified, single time step)
    let kp = params.kp.unwrap_or(1.0);
    let ki = params.ki.unwrap_or(0.0);
    let kd = params.kd.unwrap_or(0.0);

    let p_term = kp * error;

    additional.insert("error".to_string(), error);
    additional.insert("proportional_term".to_string(), p_term);
    additional.insert("kp".to_string(), kp);

    // Ziegler-Nichols tuning suggestion (for first-order plus dead time)
    if let Some(tau) = params.time_constant {
        let kp_suggested = 0.9 * tau;
        let ti_suggested = 3.33 * tau;
        let td_suggested = 0.0; // Conservative

        additional.insert("zn_kp_suggested".to_string(), kp_suggested);
        additional.insert("zn_ti_suggested".to_string(), ti_suggested);
        additional.insert("zn_ki_suggested".to_string(), kp_suggested / ti_suggested);

        return Ok(EngineeringResult {
            value: kp_suggested,
            unit: "gain".to_string(),
            formula_used: "Ziegler-Nichols PID tuning".to_string(),
            classification: Some("PI control recommended".to_string()),
            interpretation: format!("Suggested Kp={:.2}, Ti={:.1}s", kp_suggested, ti_suggested),
            additional: Some(additional),
        });
    }

    Ok(EngineeringResult {
        value: p_term,
        unit: "control output".to_string(),
        formula_used: "PID: u = Kp·e + Ki·∫e + Kd·de/dt".to_string(),
        classification: None,
        interpretation: format!("Error={:.2}, P-term={:.2}", error, p_term),
        additional: Some(additional),
    })
}

impl Default for EngineeringParams {
    fn default() -> Self {
        Self {
            pressure_rms: None,
            reference_pressure: None,
            velocity_source: None,
            velocity_observer: None,
            sound_speed: None,
            frequency: None,
            room_volume: None,
            absorption_coefficient: None,
            stress: None,
            strain: None,
            youngs_modulus: None,
            yield_strength: None,
            crack_length: None,
            geometry_factor: None,
            velocity: None,
            pressure_1: None,
            pressure_2: None,
            height_1: None,
            height_2: None,
            density: None,
            viscosity: None,
            radius: None,
            length: None,
            flow_rate: None,
            drag_coefficient: None,
            cross_sectional_area: None,
            kp: None,
            ki: None,
            kd: None,
            setpoint: None,
            process_variable: None,
            time_constant: None,
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    // Acoustics Tests (6 tests)
    #[test]
    fn test_acoustics_spl_quiet() {
        let params = EngineeringParams {
            pressure_rms: Some(0.002), // 2 mPa
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // SPL = 20·log₁₀(0.002/20e-6) = 20·log₁₀(100) = 40 dB
        assert!((result.value - 40.0).abs() < 1.0);
        assert!(result.classification.unwrap().contains("Quiet"));
    }

    #[test]
    fn test_acoustics_spl_loud() {
        let params = EngineeringParams {
            pressure_rms: Some(0.2), // 200 mPa
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // SPL = 20·log₁₀(0.2/20e-6) = 20·log₁₀(10000) = 80 dB
        assert!((result.value - 80.0).abs() < 1.0);
        assert!(result.classification.unwrap().contains("Loud"));
    }

    #[test]
    fn test_acoustics_doppler_approaching() {
        let params = EngineeringParams {
            frequency: Some(1000.0),     // 1 kHz
            sound_speed: Some(343.0),    // Speed of sound in air
            velocity_source: Some(20.0), // Source moving toward observer (positive)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // f' = f·(v + v_obs)/(v - v_src) = 1000·343/(343-20) = 1000·343/323 ≈ 1062 Hz
        assert!(result.value > 1000.0); // Frequency should increase
    }

    #[test]
    fn test_acoustics_doppler_receding() {
        let params = EngineeringParams {
            frequency: Some(1000.0),
            sound_speed: Some(343.0),
            velocity_source: Some(-20.0), // Source moving away (negative)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // f' = 1000·343/(343-(-20)) = 1000·343/363 ≈ 945 Hz
        assert!(result.value < 1000.0); // Frequency should decrease
    }

    #[test]
    fn test_acoustics_reverberation_dry() {
        let params = EngineeringParams {
            room_volume: Some(100.0),          // 100 m³ room
            absorption_coefficient: Some(0.8), // High absorption (carpet, curtains)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // RT60 should be short for dry room
        assert!(result.value < 1.0);
        assert!(result.classification.is_some());
    }

    #[test]
    fn test_acoustics_reverberation_reverberant() {
        let params = EngineeringParams {
            room_volume: Some(10000.0),        // Large cathedral
            absorption_coefficient: Some(0.1), // Low absorption (stone walls)
            ..Default::default()
        };
        let result = calculate_acoustics(&params).unwrap();
        // RT60 should be long for reverberant space
        assert!(result.value > 2.0);
    }

    // Materials Tests (5 tests)
    #[test]
    fn test_materials_hookes_law_elastic() {
        let params = EngineeringParams {
            youngs_modulus: Some(200e9), // Steel: 200 GPa
            strain: Some(0.0005),        // 0.05% strain (lower to stay in safe zone)
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // σ = E·ε = 200e9 · 0.0005 = 100e6 Pa = 100 MPa
        assert!((result.value - 100.0).abs() < 1.0);
        assert!(result.classification.unwrap().contains("safe"));
    }

    #[test]
    fn test_materials_hookes_law_strain_from_stress() {
        let params = EngineeringParams {
            youngs_modulus: Some(70e9), // Aluminum: 70 GPa
            stress: Some(140e6),        // 140 MPa
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // ε = σ/E = 140e6/70e9 = 0.002 = 0.2%
        assert!((result.value - 0.2).abs() < 0.01);
    }

    #[test]
    fn test_materials_fracture_mechanics_safe() {
        let params = EngineeringParams {
            stress: Some(100e6),       // 100 MPa
            crack_length: Some(0.001), // 1 mm crack
            geometry_factor: Some(1.0),
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // K = Y·σ·√(π·a) = 1.0 · 100e6 · √(π·0.001) ≈ 5.6 MPa·√m
        assert!(result.value < 20.0); // Safe range
        assert!(result.classification.unwrap().contains("Safe"));
    }

    #[test]
    fn test_materials_fracture_mechanics_critical() {
        let params = EngineeringParams {
            stress: Some(500e6),         // 500 MPa (very high)
            crack_length: Some(0.01),    // 10 mm crack
            geometry_factor: Some(1.12), // Edge crack
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // K = 1.12 · 500e6 · √(π·0.01) ≈ 99 MPa·√m (very high!)
        assert!(result.value > 40.0);
    }

    #[test]
    fn test_materials_yield_strength_check() {
        let params = EngineeringParams {
            youngs_modulus: Some(200e9),
            strain: Some(0.0025),        // 0.25% strain
            yield_strength: Some(250e6), // 250 MPa yield
            ..Default::default()
        };
        let result = calculate_materials(&params).unwrap();
        // σ = 200e9 · 0.0025 = 500e6 Pa = 500 MPa (exceeds yield!)
        assert!(result.value > 400.0);
        assert!(result.classification.unwrap().contains("Plastic"));
    }

    // Fluid Mechanics Tests (5 tests)
    #[test]
    fn test_fluid_bernoulli_exit_velocity() {
        let params = EngineeringParams {
            density: Some(1000.0),      // Water
            pressure_1: Some(200000.0), // 200 kPa
            pressure_2: Some(101325.0), // Atmospheric
            velocity: Some(1.0),        // 1 m/s inlet
            height_1: Some(0.0),
            height_2: Some(0.0), // Same height
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // P₁ + ½ρv₁² = P₂ + ½ρv₂²
        // 200000 + 500 = 101325 + 500·v₂²
        // v₂² = (200000 + 500 - 101325) / 500 ≈ 198.35
        // v₂ ≈ 14.08 m/s
        assert!(result.value > 10.0 && result.value < 20.0);
    }

    #[test]
    fn test_fluid_poiseuille_laminar() {
        let params = EngineeringParams {
            pressure_1: Some(1000.0), // 1 kPa pressure drop
            radius: Some(0.001),      // 1 mm radius pipe (smaller for laminar flow)
            length: Some(1.0),        // 1 meter long
            viscosity: Some(0.001),   // Water viscosity
            density: Some(1000.0),
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // Q = (π·ΔP·r⁴)/(8·η·L), Re = 250 (laminar)
        assert!(result.value > 0.0); // Positive flow rate
        assert!(result.classification.unwrap().contains("Laminar"));
    }

    #[test]
    fn test_fluid_drag_force() {
        let params = EngineeringParams {
            density: Some(1.225),            // Air at sea level
            velocity: Some(30.0),            // 30 m/s ≈ 108 km/h
            drag_coefficient: Some(0.3),     // Streamlined car
            cross_sectional_area: Some(2.0), // 2 m² frontal area
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // F_D = ½·ρ·v²·C_D·A = 0.5·1.225·900·0.3·2.0 = 330.75 N
        assert!((result.value - 330.75).abs() < 10.0);
    }

    #[test]
    fn test_fluid_drag_high_speed() {
        let params = EngineeringParams {
            density: Some(1.225),
            velocity: Some(100.0),       // 100 m/s ≈ 360 km/h
            drag_coefficient: Some(0.5), // Less streamlined
            cross_sectional_area: Some(1.5),
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        // F_D = 0.5·1.225·10000·0.5·1.5 = 4593.75 N
        assert!(result.value > 4000.0);
    }

    #[test]
    fn test_fluid_poiseuille_reynolds() {
        let params = EngineeringParams {
            pressure_1: Some(1000.0),
            radius: Some(0.01), // 10 mm radius
            length: Some(2.0),
            viscosity: Some(0.001),
            density: Some(1000.0),
            ..Default::default()
        };
        let result = calculate_fluid_mechanics(&params).unwrap();
        let additional = result.additional.unwrap();
        // Should calculate Reynolds number
        assert!(additional.contains_key("reynolds_number"));
    }

    // Control Theory Tests (4 tests)
    #[test]
    fn test_control_pid_proportional() {
        let params = EngineeringParams {
            setpoint: Some(100.0),
            process_variable: Some(80.0),
            kp: Some(2.0),
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // error = 100 - 80 = 20
        // P-term = 2.0 · 20 = 40
        assert!((result.value - 40.0).abs() < 0.1);
    }

    #[test]
    fn test_control_pid_negative_error() {
        let params = EngineeringParams {
            setpoint: Some(50.0),
            process_variable: Some(70.0), // Overshoot
            kp: Some(1.5),
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // error = 50 - 70 = -20
        // P-term = 1.5 · (-20) = -30
        assert!((result.value - (-30.0)).abs() < 0.1);
    }

    #[test]
    fn test_control_ziegler_nichols_tuning() {
        let params = EngineeringParams {
            setpoint: Some(100.0),
            process_variable: Some(100.0), // At setpoint
            time_constant: Some(10.0),     // τ = 10 seconds
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // Kp = 0.9·τ = 9.0
        assert!((result.value - 9.0).abs() < 0.1);
        assert!(result.interpretation.contains("Kp"));
    }

    #[test]
    fn test_control_pid_zero_error() {
        let params = EngineeringParams {
            setpoint: Some(75.0),
            process_variable: Some(75.0), // Perfect match
            kp: Some(3.0),
            ..Default::default()
        };
        let result = calculate_control(&params).unwrap();
        // error = 0, P-term = 0
        assert!((result.value).abs() < 0.001);
    }
}
