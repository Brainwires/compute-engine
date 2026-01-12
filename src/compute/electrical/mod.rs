//! Electrical Circuit Analysis Module
//!
//! Comprehensive circuit analysis for DC and AC circuits, including:
//! - DC circuits (Ohm's law, series/parallel, voltage/current dividers)
//! - AC circuits (impedance, phasors, RLC analysis)
//! - Circuit theorems (Thévenin, Norton, superposition)
//! - Network analysis (mesh, nodal, loop analysis)
//! - Transfer functions and frequency response
//! - Transient analysis (RL, RC, RLC circuits)
//! - Power calculations (real, reactive, apparent, power factor)

pub mod dc_analysis;
pub mod ac_analysis;
pub mod impedance;
pub mod network_analysis;
pub mod transient_analysis;
pub mod power_analysis;
pub mod filter_design;
pub mod nec_calculations;

pub use dc_analysis::*;
pub use ac_analysis::*;
pub use impedance::*;
pub use network_analysis::*;
pub use transient_analysis::*;
pub use power_analysis::*;
pub use filter_design::*;
pub use nec_calculations::*;

use num_complex::Complex64;

/// Complex number for AC analysis (voltage, current, impedance)
pub type ComplexValue = Complex64;

/// Represents a circuit component
#[derive(Debug, Clone)]
pub enum Component {
    Resistor { resistance: f64 },
    Capacitor { capacitance: f64 },
    Inductor { inductance: f64 },
    VoltageSource { voltage: f64, phase: f64 },
    CurrentSource { current: f64, phase: f64 },
}

/// Circuit node for network analysis
#[derive(Debug, Clone)]
pub struct Node {
    pub id: usize,
    pub voltage: Option<f64>,
}

/// Circuit branch connecting two nodes
#[derive(Debug, Clone)]
pub struct Branch {
    pub from_node: usize,
    pub to_node: usize,
    pub component: Component,
}

#[cfg(test)]
#[path = "../../../tests/unit/compute/electrical/electrical_tests.rs"]
mod tests;
