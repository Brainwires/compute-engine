//! Network Analysis
//!
//! Circuit network analysis including Thévenin/Norton equivalents,
//! mesh analysis, nodal analysis, and superposition theorem.

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// Thévenin equivalent circuit
#[derive(Debug, Clone)]
pub struct TheveninEquivalent {
    /// Thévenin voltage (open-circuit voltage)
    pub v_th: f64,
    /// Thévenin resistance (equivalent resistance)
    pub r_th: f64,
}

/// Norton equivalent circuit
#[derive(Debug, Clone)]
pub struct NortonEquivalent {
    /// Norton current (short-circuit current)
    pub i_n: f64,
    /// Norton resistance (equivalent resistance)
    pub r_n: f64,
}

/// Convert Thévenin to Norton equivalent
///
/// I_N = V_TH / R_TH
/// R_N = R_TH
pub fn thevenin_to_norton(thevenin: &TheveninEquivalent) -> NortonEquivalent {
    NortonEquivalent {
        i_n: thevenin.v_th / thevenin.r_th,
        r_n: thevenin.r_th,
    }
}

/// Convert Norton to Thévenin equivalent
///
/// V_TH = I_N * R_N
/// R_TH = R_N
pub fn norton_to_thevenin(norton: &NortonEquivalent) -> TheveninEquivalent {
    TheveninEquivalent {
        v_th: norton.i_n * norton.r_n,
        r_th: norton.r_n,
    }
}

/// Calculate maximum power transfer to load
///
/// P_max = V_TH² / (4 * R_TH)  when R_L = R_TH
pub fn maximum_power_transfer(v_th: f64, r_th: f64) -> f64 {
    v_th * v_th / (4.0 * r_th)
}

/// Calculate optimal load resistance for maximum power transfer
///
/// R_L = R_TH for maximum power transfer
pub fn optimal_load_resistance(r_th: f64) -> f64 {
    r_th
}

/// Calculate power delivered to load given Thévenin equivalent
///
/// P_L = (V_TH² * R_L) / (R_TH + R_L)²
pub fn power_to_load(v_th: f64, r_th: f64, r_load: f64) -> f64 {
    let i = v_th / (r_th + r_load);
    i * i * r_load
}

/// Solve DC mesh analysis using mesh equations
///
/// Solves: [R] * [I] = [V]
/// where [R] is the resistance matrix, [I] is mesh currents, [V] is voltage sources
///
/// # Arguments
/// * `resistance_matrix` - N×N matrix of resistances
/// * `voltage_vector` - N-element vector of voltage sources
///
/// # Returns
/// Vector of mesh currents
///
/// # Example
/// For two meshes with R1=2Ω, R2=3Ω, R3=1Ω (shared), V1=10V, V2=5V:
/// ```
/// use nalgebra::{dmatrix, dvector};
/// use computational_engine::electrical::network_analysis::mesh_analysis;
///
/// let r = dmatrix![
///     5.0, -1.0;  // Mesh 1: R1 + R3, -R3
///     -1.0, 4.0   // Mesh 2: -R3, R2 + R3
/// ];
/// let v = dvector![10.0, 5.0];
/// let currents = mesh_analysis(&r, &v);
/// ```
pub fn mesh_analysis(resistance_matrix: &DMatrix<f64>, voltage_vector: &DVector<f64>) -> DVector<f64> {
    resistance_matrix
        .clone()
        .lu()
        .solve(voltage_vector)
        .expect("Failed to solve mesh equations")
}

/// Solve DC nodal analysis using nodal equations
///
/// Solves: [G] * [V] = [I]
/// where [G] is the conductance matrix, [V] is node voltages, [I] is current sources
///
/// # Arguments
/// * `conductance_matrix` - N×N matrix of conductances
/// * `current_vector` - N-element vector of current sources
///
/// # Returns
/// Vector of node voltages
pub fn nodal_analysis(
    conductance_matrix: &DMatrix<f64>,
    current_vector: &DVector<f64>,
) -> DVector<f64> {
    conductance_matrix
        .clone()
        .lu()
        .solve(current_vector)
        .expect("Failed to solve nodal equations")
}

/// Calculate Y-parameters (admittance parameters) for two-port network
///
/// I1 = Y11*V1 + Y12*V2
/// I2 = Y21*V1 + Y22*V2
#[derive(Debug, Clone)]
pub struct YParameters {
    pub y11: Complex64,
    pub y12: Complex64,
    pub y21: Complex64,
    pub y22: Complex64,
}

/// Calculate Z-parameters (impedance parameters) for two-port network
///
/// V1 = Z11*I1 + Z12*I2
/// V2 = Z21*I1 + Z22*I2
#[derive(Debug, Clone)]
pub struct ZParameters {
    pub z11: Complex64,
    pub z12: Complex64,
    pub z21: Complex64,
    pub z22: Complex64,
}

/// Convert Y-parameters to Z-parameters
pub fn y_to_z_parameters(y: &YParameters) -> ZParameters {
    let det_y = y.y11 * y.y22 - y.y12 * y.y21;
    ZParameters {
        z11: y.y22 / det_y,
        z12: -y.y12 / det_y,
        z21: -y.y21 / det_y,
        z22: y.y11 / det_y,
    }
}

/// Convert Z-parameters to Y-parameters
pub fn z_to_y_parameters(z: &ZParameters) -> YParameters {
    let det_z = z.z11 * z.z22 - z.z12 * z.z21;
    YParameters {
        y11: z.z22 / det_z,
        y12: -z.z12 / det_z,
        y21: -z.z21 / det_z,
        y22: z.z11 / det_z,
    }
}

/// Calculate ABCD parameters (transmission parameters) from Z-parameters
///
/// V1 = A*V2 - B*I2
/// I1 = C*V2 - D*I2
#[derive(Debug, Clone)]
pub struct ABCDParameters {
    pub a: Complex64,
    pub b: Complex64,
    pub c: Complex64,
    pub d: Complex64,
}

pub fn z_to_abcd_parameters(z: &ZParameters) -> ABCDParameters {
    ABCDParameters {
        a: z.z11 / z.z21,
        b: (z.z11 * z.z22 - z.z12 * z.z21) / z.z21,
        c: Complex64::new(1.0, 0.0) / z.z21,
        d: z.z22 / z.z21,
    }
}

/// Calculate voltage gain of two-port network
///
/// A_v = V2 / V1 = Z21 / (Z11 + Z_source) for voltage amplifier
pub fn voltage_gain(z_params: &ZParameters, z_source: Complex64, z_load: Complex64) -> Complex64 {
    let det_z = z_params.z11 * z_params.z22 - z_params.z12 * z_params.z21;
    let numerator = z_params.z21 * z_load;
    let denominator = (z_params.z11 + z_source) * (z_params.z22 + z_load) - det_z;
    numerator / denominator
}

/// Calculate current gain of two-port network
///
/// A_i = I2 / I1
pub fn current_gain(z_params: &ZParameters, _z_source: Complex64, z_load: Complex64) -> Complex64 {
    let numerator = -z_params.z21;
    let denominator = z_params.z22 + z_load;
    numerator / denominator
}

/// Calculate power gain of two-port network (in dB)
pub fn power_gain_db(voltage_gain: Complex64, z_in: Complex64, z_out: Complex64) -> f64 {
    let power_ratio = voltage_gain.norm_sqr() * (z_in.re / z_out.re);
    10.0 * power_ratio.log10()
}

/// Calculate input impedance of loaded two-port network
///
/// Z_in = Z11 - (Z12 * Z21) / (Z22 + Z_L)
pub fn input_impedance_two_port(z_params: &ZParameters, z_load: Complex64) -> Complex64 {
    z_params.z11 - (z_params.z12 * z_params.z21) / (z_params.z22 + z_load)
}

/// Calculate output impedance of two-port network
///
/// Z_out = Z22 - (Z12 * Z21) / (Z11 + Z_S)
pub fn output_impedance_two_port(z_params: &ZParameters, z_source: Complex64) -> Complex64 {
    z_params.z22 - (z_params.z12 * z_params.z21) / (z_params.z11 + z_source)
}

/// Calculate delta-to-wye (Δ → Y) transformation
///
/// # Arguments
/// * `r_ab`, `r_bc`, `r_ca` - Delta resistances
///
/// # Returns
/// (R_a, R_b, R_c) - Wye resistances
pub fn delta_to_wye(r_ab: f64, r_bc: f64, r_ca: f64) -> (f64, f64, f64) {
    let sum = r_ab + r_bc + r_ca;
    let r_a = (r_ab * r_ca) / sum;
    let r_b = (r_ab * r_bc) / sum;
    let r_c = (r_bc * r_ca) / sum;
    (r_a, r_b, r_c)
}

/// Calculate wye-to-delta (Y → Δ) transformation
///
/// # Arguments
/// * `r_a`, `r_b`, `r_c` - Wye resistances
///
/// # Returns
/// (R_ab, R_bc, R_ca) - Delta resistances
pub fn wye_to_delta(r_a: f64, r_b: f64, r_c: f64) -> (f64, f64, f64) {
    let sum_products = r_a * r_b + r_b * r_c + r_c * r_a;
    let r_ab = sum_products / r_c;
    let r_bc = sum_products / r_a;
    let r_ca = sum_products / r_b;
    (r_ab, r_bc, r_ca)
}

