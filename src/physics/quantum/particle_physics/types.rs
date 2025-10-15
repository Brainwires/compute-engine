//! Particle physics type definitions

use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Particle {
    pub id: String,
    pub particle_type: ParticleType,
    pub mass: f64, // GeV/c²
    pub charge: f64, // in units of elementary charge
    pub spin: f64, // ℏ units
    pub position: Vector3D,
    pub momentum: Vector3D,
    pub energy: f64, // GeV
    pub lifetime: Option<f64>, // seconds
    pub quantum_numbers: QuantumNumbers,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub enum ParticleType {
    Electron,
    Muon,
    Tau,
    ElectronNeutrino,
    MuonNeutrino,
    TauNeutrino,
    Up,
    Down,
    Charm,
    Strange,
    Top,
    Bottom,
    Photon,
    W,
    Z,
    Gluon,
    Higgs,
    Proton,
    Neutron,
    Pion,
    Kaon,
    Custom(String),
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct QuantumNumbers {
    pub baryon_number: f64,
    pub lepton_number: f64,
    pub strangeness: f64,
    pub charm: f64,
    pub beauty: f64,
    pub truth: f64,
    pub isospin: f64,
    pub color_charge: Option<String>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ParticleInteractionResult {
    pub initial_particles: Vec<Particle>,
    pub final_particles: Vec<Particle>,
    pub interaction_type: String,
    pub conservation_analysis: ConservationAnalysis,
    pub cross_section: f64, // in barns
    pub matrix_element: Option<Complex64>,
    pub feynman_diagrams: Vec<FeynmanDiagram>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ConservationAnalysis {
    pub energy_conserved: bool,
    pub momentum_conserved: bool,
    pub charge_conserved: bool,
    pub baryon_number_conserved: bool,
    pub lepton_number_conserved: bool,
    pub violations: Vec<String>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct DecayChannel {
    pub parent: ParticleType,
    pub products: Vec<ParticleType>,
    pub branching_ratio: f64,
    pub lifetime: f64,
    pub decay_mode: String,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct FeynmanDiagram {
    pub vertices: Vec<InteractionVertex>,
    pub propagators: Vec<Propagator>,
    pub external_lines: Vec<ExternalLine>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct InteractionVertex {
    pub vertex_type: String, // QED, QCD, weak, etc.
    pub incoming: Vec<String>,
    pub outgoing: Vec<String>,
    pub coupling_constant: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Propagator {
    pub particle_type: ParticleType,
    pub momentum: FourMomentum,
    pub propagator_type: String, // "fermion", "boson", "photon"
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ExternalLine {
    pub particle: Particle,
    pub direction: String, // "incoming" or "outgoing"
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct FourMomentum {
    pub e: f64,  // Energy
    pub px: f64, // x-momentum
    pub py: f64, // y-momentum
    pub pz: f64, // z-momentum
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SimulationStep {
    pub time: f64,
    pub particles: Vec<Particle>,
    pub interactions: Vec<String>,
    pub total_energy: f64,
    pub total_momentum: Vector3D,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ParticleState {
    pub particle: Particle,
    pub wave_function: Vec<Complex64>,
    pub probability_density: Vec<f64>,
}
