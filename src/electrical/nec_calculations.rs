//! NEC (National Electrical Code) Calculations
//!
//! Practical electrical calculations for electricians and contractors
//! based on NEC requirements and industry standards.

use std::f64::consts::PI;

/// Calculate conductor ampacity adjustment for ambient temperature
///
/// Based on NEC Table 310.15(B)(2)(a) - Ambient Temperature Correction Factors
///
/// # Arguments
/// * `base_ampacity` - Base ampacity from NEC tables
/// * `ambient_temp_c` - Ambient temperature in Celsius
/// * `conductor_temp_rating` - Conductor temperature rating (60, 75, or 90°C)
///
/// # Returns
/// Adjusted ampacity in amperes
pub fn temperature_corrected_ampacity(
    base_ampacity: f64,
    ambient_temp_c: f64,
    conductor_temp_rating: u32,
) -> f64 {
    let correction_factor = match conductor_temp_rating {
        60 => {
            if ambient_temp_c <= 30.0 {
                1.0
            } else if ambient_temp_c <= 35.0 {
                0.94
            } else if ambient_temp_c <= 40.0 {
                0.82
            } else if ambient_temp_c <= 45.0 {
                0.71
            } else if ambient_temp_c <= 50.0 {
                0.58
            } else {
                0.41
            }
        }
        75 => {
            if ambient_temp_c <= 30.0 {
                1.0
            } else if ambient_temp_c <= 35.0 {
                0.96
            } else if ambient_temp_c <= 40.0 {
                0.88
            } else if ambient_temp_c <= 45.0 {
                0.82
            } else if ambient_temp_c <= 50.0 {
                0.75
            } else {
                0.67
            }
        }
        90 => {
            if ambient_temp_c <= 30.0 {
                1.0
            } else if ambient_temp_c <= 35.0 {
                0.96
            } else if ambient_temp_c <= 40.0 {
                0.91
            } else if ambient_temp_c <= 45.0 {
                0.87
            } else if ambient_temp_c <= 50.0 {
                0.82
            } else {
                0.76
            }
        }
        _ => 1.0,
    };

    base_ampacity * correction_factor
}

/// Calculate ampacity adjustment for number of conductors in conduit
///
/// Based on NEC Table 310.15(B)(3)(a) - Adjustment Factors
///
/// # Arguments
/// * `base_ampacity` - Base ampacity
/// * `num_current_carrying` - Number of current-carrying conductors
///
/// # Returns
/// Adjusted ampacity
pub fn conduit_fill_ampacity(base_ampacity: f64, num_current_carrying: u32) -> f64 {
    let adjustment_factor = match num_current_carrying {
        1..=3 => 1.0,
        4..=6 => 0.8,
        7..=9 => 0.7,
        10..=20 => 0.5,
        21..=30 => 0.45,
        31..=40 => 0.4,
        _ => 0.35,
    };

    base_ampacity * adjustment_factor
}

/// Calculate voltage drop in a circuit
///
/// VD = 2 * I * R * L (for single-phase)
/// VD = √3 * I * R * L (for three-phase)
///
/// # Arguments
/// * `current` - Load current in amperes
/// * `length_feet` - One-way length in feet
/// * `conductor_resistance` - Conductor resistance in ohms per 1000 feet
/// * `is_three_phase` - True for three-phase, false for single-phase
///
/// # Returns
/// Voltage drop in volts
pub fn voltage_drop(
    current: f64,
    length_feet: f64,
    conductor_resistance_per_1000ft: f64,
    is_three_phase: bool,
) -> f64 {
    let resistance = conductor_resistance_per_1000ft * length_feet / 1000.0;
    if is_three_phase {
        3_f64.sqrt() * current * resistance
    } else {
        2.0 * current * resistance
    }
}

/// Calculate percent voltage drop
///
/// %VD = (Voltage Drop / Source Voltage) * 100
///
/// NEC recommends ≤3% for branch circuits, ≤5% total (feeder + branch)
pub fn percent_voltage_drop(voltage_drop: f64, source_voltage: f64) -> f64 {
    (voltage_drop / source_voltage) * 100.0
}

/// Calculate minimum conductor size for given voltage drop
///
/// Returns circular mils (CM) required
pub fn conductor_size_for_voltage_drop(
    current: f64,
    length_feet: f64,
    source_voltage: f64,
    max_percent_drop: f64,
    is_three_phase: bool,
) -> f64 {
    let max_vd = source_voltage * max_percent_drop / 100.0;
    let k = if is_three_phase { 3_f64.sqrt() } else { 2.0 };

    // CM = (K * I * L) / VD
    // K = 12.9 for copper, 21.2 for aluminum
    let k_copper = 12.9;
    (k_copper * k * current * length_feet) / max_vd
}

/// Convert circular mils to AWG size (approximation)
///
/// Returns approximate AWG gauge number
pub fn circular_mils_to_awg(circular_mils: f64) -> u32 {
    // AWG formula: d(n) = 0.005 * 92^((36-n)/39) inches
    // CM = π/4 * (d * 1000)^2
    let diameter_inches = (circular_mils * 4.0 / (PI * 1_000_000.0)).sqrt();
    let awg = 36.0 - 39.0 * (diameter_inches / 0.005).log2() / 92_f64.log2();
    awg.round() as u32
}

/// Calculate conduit fill percentage
///
/// NEC Article 344.22 - Maximum fill: 53% for 2 conductors, 31% for 3+
///
/// # Arguments
/// * `num_conductors` - Number of conductors
/// * `conductor_area_sqin` - Cross-sectional area of each conductor
/// * `conduit_area_sqin` - Internal area of conduit
///
/// # Returns
/// Fill percentage
pub fn conduit_fill_percentage(
    num_conductors: u32,
    conductor_area_sqin: f64,
    conduit_area_sqin: f64,
) -> f64 {
    let total_conductor_area = num_conductors as f64 * conductor_area_sqin;
    (total_conductor_area / conduit_area_sqin) * 100.0
}

/// Check if conduit fill meets NEC requirements
pub fn is_conduit_fill_compliant(fill_percentage: f64, num_conductors: u32) -> bool {
    match num_conductors {
        1 => fill_percentage <= 53.0,
        2 => fill_percentage <= 31.0,
        _ => fill_percentage <= 40.0,
    }
}

/// Calculate branch circuit load (continuous + non-continuous)
///
/// NEC 210.19(A)(1) - Branch circuits shall have an ampacity of not less than
/// the maximum load to be served + 125% of continuous load
///
/// # Arguments
/// * `continuous_load` - Continuous load in amperes (≥3 hours)
/// * `non_continuous_load` - Non-continuous load in amperes
///
/// # Returns
/// Minimum required ampacity
pub fn branch_circuit_ampacity(continuous_load: f64, non_continuous_load: f64) -> f64 {
    (continuous_load * 1.25) + non_continuous_load
}

/// Calculate feeder load with demand factors
///
/// # Arguments
/// * `connected_load` - Total connected load in amperes
/// * `demand_factor` - Demand factor (0.0 to 1.0)
///
/// # Returns
/// Feeder load in amperes
pub fn feeder_demand_load(connected_load: f64, demand_factor: f64) -> f64 {
    connected_load * demand_factor
}

/// Calculate general lighting load (NEC Table 220.12)
///
/// # Arguments
/// * `area_sq_ft` - Floor area in square feet
/// * `building_type` - "dwelling", "office", "warehouse", etc.
///
/// # Returns
/// Lighting load in volt-amperes
pub fn general_lighting_load(area_sq_ft: f64, building_type: &str) -> f64 {
    let va_per_sqft = match building_type {
        "dwelling" => 3.0,
        "office" => 3.5,
        "warehouse" => 0.25,
        "school" => 3.0,
        "store" => 3.0,
        _ => 3.0, // Default
    };
    area_sq_ft * va_per_sqft
}

/// Calculate residential dwelling unit load (NEC 220.82)
///
/// Standard method for single-family dwellings
///
/// # Arguments
/// * `area_sq_ft` - Living area
/// * `num_small_appliance_circuits` - Number of 20A small appliance circuits
/// * `num_laundry_circuits` - Number of laundry circuits
/// * `largest_motor_hp` - Largest motor load in horsepower
/// * `additional_loads_va` - Other loads (HVAC, water heater, etc.)
///
/// # Returns
/// Total calculated load in volt-amperes
pub fn dwelling_unit_load_calculation(
    area_sq_ft: f64,
    num_small_appliance_circuits: u32,
    num_laundry_circuits: u32,
    largest_motor_hp: f64,
    additional_loads_va: f64,
) -> f64 {
    // General lighting: 3 VA/sq ft
    let lighting = area_sq_ft * 3.0;

    // Small appliance: 1500 VA per circuit
    let small_appliance = (num_small_appliance_circuits as f64) * 1500.0;

    // Laundry: 1500 VA per circuit
    let laundry = (num_laundry_circuits as f64) * 1500.0;

    // First 10,000 VA at 100%
    let first_10k = 10_000.0_f64.min(lighting + small_appliance + laundry);

    // Remainder at 40%
    let remainder = (lighting + small_appliance + laundry - 10_000.0).max(0.0) * 0.4;

    // Motor load (add 25% of largest motor)
    let motor_va = largest_motor_hp * 746.0; // 1 HP = 746 W
    let motor_load = motor_va * 1.25;

    first_10k + remainder + motor_load + additional_loads_va
}

/// Calculate required service size in amperes
///
/// # Arguments
/// * `total_load_va` - Total calculated load in volt-amperes
/// * `voltage` - Service voltage (typically 120/240V or 120/208V)
///
/// # Returns
/// Minimum service size in amperes
pub fn minimum_service_size(total_load_va: f64, voltage: f64) -> f64 {
    total_load_va / voltage
}

/// Round service size up to standard rating
///
/// Standard service sizes: 100A, 125A, 150A, 200A, 225A, 400A, 600A, etc.
pub fn standard_service_size(calculated_amperes: f64) -> u32 {
    let standard_sizes = [100, 125, 150, 200, 225, 400, 600, 800, 1000, 1200];
    for &size in &standard_sizes {
        if calculated_amperes <= size as f64 {
            return size;
        }
    }
    1200 // Maximum
}

/// Calculate motor full-load current (approximate, for copper wire sizing)
///
/// Based on NEC Tables 430.247-250
///
/// # Arguments
/// * `horsepower` - Motor horsepower
/// * `voltage` - Motor voltage
/// * `is_three_phase` - True for three-phase motor
///
/// # Returns
/// Full-load current in amperes
pub fn motor_full_load_current(horsepower: f64, voltage: u32, is_three_phase: bool) -> f64 {
    // Simplified approximations
    match (voltage, is_three_phase) {
        (115, false) => horsepower * 10.0,
        (230, false) => horsepower * 5.0,
        (208, true) => horsepower * 2.5,
        (230, true) => horsepower * 2.3,
        (460, true) => horsepower * 1.15,
        _ => horsepower * 5.0, // Default
    }
}

/// Calculate overcurrent protection device size for motors
///
/// NEC 430.52 - Motor branch-circuit short-circuit and ground-fault protection
/// Typically 175% to 250% of FLA depending on motor type
///
/// # Arguments
/// * `full_load_amps` - Motor full-load amperes
/// * `motor_type` - "ac", "dc", or "wound_rotor"
///
/// # Returns
/// Maximum OCPD rating in amperes
pub fn motor_ocpd_size(full_load_amps: f64, motor_type: &str) -> f64 {
    let multiplier = match motor_type {
        "ac" => 2.5,           // 250% for non-time delay fuse
        "dc" => 1.5,           // 150%
        "wound_rotor" => 1.5,  // 150%
        _ => 2.5,
    };
    full_load_amps * multiplier
}

/// Calculate transformer secondary current
///
/// # Arguments
/// * `kva` - Transformer rating in kVA
/// * `secondary_voltage` - Secondary voltage
/// * `is_three_phase` - True for three-phase
///
/// # Returns
/// Secondary current in amperes
pub fn transformer_secondary_current(kva: f64, secondary_voltage: f64, is_three_phase: bool) -> f64 {
    let va = kva * 1000.0;
    if is_three_phase {
        va / (3_f64.sqrt() * secondary_voltage)
    } else {
        va / secondary_voltage
    }
}

/// Calculate transformer primary overcurrent protection
///
/// NEC 450.3(B) - Primary OCPD for transformers
/// 125% for transformers 1000V and below
pub fn transformer_primary_ocpd(primary_current: f64) -> f64 {
    primary_current * 1.25
}

