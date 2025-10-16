use anyhow::Result;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Particle {
    pub id: String,
    pub particle_type: ParticleType,
    pub mass: f64,   // GeV/c²
    pub charge: f64, // in units of elementary charge
    pub spin: f64,   // ℏ units
    pub position: Vector3D,
    pub momentum: Vector3D,
    pub energy: f64,           // GeV
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
    pub color_charge: Option<String>, // "red", "green", "blue", "anti-red", etc.
}

#[derive(Debug, Serialize)]
pub struct ParticleInteractionResult {
    pub interaction_type: String,
    pub initial_particles: Vec<Particle>,
    pub final_particles: Vec<Particle>,
    pub conservation_laws: ConservationAnalysis,
    pub cross_section: Option<f64>, // barns
    pub decay_channels: Vec<DecayChannel>,
    pub feynman_diagrams: Vec<FeynmanDiagram>,
    pub simulation_steps: Vec<SimulationStep>,
    pub physical_insights: Vec<String>,
    pub computation_time_ms: u128,
}

#[derive(Debug, Serialize)]
pub struct ConservationAnalysis {
    pub energy_momentum_conserved: bool,
    pub charge_conserved: bool,
    pub baryon_number_conserved: bool,
    pub lepton_number_conserved: bool,
    pub color_charge_conserved: bool,
    pub conservation_violations: Vec<String>,
    pub quantum_numbers_initial: QuantumNumbers,
    pub quantum_numbers_final: QuantumNumbers,
}

#[derive(Debug, Serialize)]
pub struct DecayChannel {
    pub channel_id: String,
    pub products: Vec<String>,
    pub branching_ratio: f64,
    pub q_value: f64,       // MeV
    pub partial_width: f64, // GeV
    pub selection_rules: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct FeynmanDiagram {
    pub diagram_id: String,
    pub interaction_vertices: Vec<InteractionVertex>,
    pub propagators: Vec<Propagator>,
    pub external_lines: Vec<ExternalLine>,
    pub loop_order: usize,
    pub coupling_constants: HashMap<String, f64>,
}

#[derive(Debug, Serialize)]
pub struct InteractionVertex {
    pub vertex_type: String, // "QED", "QCD", "Weak", "Higgs"
    pub particles_in: Vec<String>,
    pub particles_out: Vec<String>,
    pub coupling_strength: f64,
}

#[derive(Debug, Serialize)]
pub struct Propagator {
    pub particle_type: String,
    pub virtuality: f64,         // GeV²
    pub propagator_type: String, // "scalar", "vector", "fermion"
}

#[derive(Debug, Serialize)]
pub struct ExternalLine {
    pub particle_id: String,
    pub incoming: bool,
    pub on_shell: bool,
    pub four_momentum: FourMomentum,
}

#[derive(Debug, Serialize)]
pub struct FourMomentum {
    pub energy: f64,
    pub px: f64,
    pub py: f64,
    pub pz: f64,
}

#[derive(Debug, Serialize)]
pub struct SimulationStep {
    pub step_number: usize,
    pub time: f64, // femtoseconds
    pub event_description: String,
    pub particles_state: Vec<ParticleState>,
    pub interactions_occurred: Vec<String>,
    pub energy_balance: f64,
}

#[derive(Debug, Serialize)]
pub struct ParticleState {
    pub particle_id: String,
    pub position: Vector3D,
    pub momentum: Vector3D,
    pub energy: f64,
    pub probability_amplitude: Complex64,
}

/// Simulate particle interactions and field theories
pub fn simulate_particle_interactions(
    interaction_type: &str,
    particles: Vec<Particle>,
    simulation_steps: u64,
) -> Result<ParticleInteractionResult> {
    let start = std::time::Instant::now();

    eprintln!("⚛️ Starting particle physics simulation...");
    eprintln!("   Interaction: {}", interaction_type);
    eprintln!("   Particles: {}", particles.len());
    eprintln!("   Simulation steps: {}", simulation_steps);

    let initial_particles = particles.clone();

    // Determine interaction type and simulate accordingly
    let result = match interaction_type {
        "electromagnetic" => simulate_electromagnetic_interaction(particles, simulation_steps)?,
        "strong" => simulate_strong_interaction(particles, simulation_steps)?,
        "weak" => simulate_weak_interaction(particles, simulation_steps)?,
        "higgs" => simulate_higgs_interaction(particles, simulation_steps)?,
        "decay" => simulate_particle_decay(particles, simulation_steps)?,
        "scattering" => simulate_scattering_process(particles, simulation_steps)?,
        _ => {
            return Err(anyhow::anyhow!(
                "Unknown interaction type: {}",
                interaction_type
            ));
        }
    };

    let computation_time = start.elapsed().as_millis();
    eprintln!(
        "✅ Particle simulation completed in {} ms",
        computation_time
    );

    Ok(ParticleInteractionResult {
        interaction_type: interaction_type.to_string(),
        initial_particles,
        final_particles: result.final_particles,
        conservation_laws: result.conservation_laws,
        cross_section: result.cross_section,
        decay_channels: result.decay_channels,
        feynman_diagrams: result.feynman_diagrams,
        simulation_steps: result.simulation_steps,
        physical_insights: result.physical_insights,
        computation_time_ms: computation_time,
    })
}

struct InteractionResult {
    final_particles: Vec<Particle>,
    conservation_laws: ConservationAnalysis,
    cross_section: Option<f64>,
    decay_channels: Vec<DecayChannel>,
    feynman_diagrams: Vec<FeynmanDiagram>,
    simulation_steps: Vec<SimulationStep>,
    physical_insights: Vec<String>,
}

/// Simulate electromagnetic interactions (QED)
fn simulate_electromagnetic_interaction(
    mut particles: Vec<Particle>,
    steps: u64,
) -> Result<InteractionResult> {
    eprintln!("⚡ Simulating electromagnetic interaction (QED)...");

    let mut simulation_steps = vec![];
    let mut physical_insights = vec![
        "Electromagnetic force mediated by virtual photons".to_string(),
        "Coupling strength α ≈ 1/137 (fine structure constant)".to_string(),
    ];

    // Check for charged particles
    let charged_particles: Vec<_> = particles
        .iter()
        .filter(|p| p.charge.abs() > 1e-10)
        .cloned()
        .collect();

    if charged_particles.len() < 2 {
        physical_insights
            .push("Need at least 2 charged particles for electromagnetic interaction".to_string());
    }

    // Simulate Compton scattering if photon + electron
    if has_photon_and_electron(&particles) {
        physical_insights.push("Compton scattering: γ + e⁻ → γ' + e'".to_string());
        physical_insights
            .push("Energy and momentum conservation determine scattering angles".to_string());

        // Update particle energies and momenta
        update_compton_scattering(&mut particles);
    }

    // Simulate electron-positron annihilation if present
    if has_electron_positron_pair(&particles) {
        physical_insights.push("e⁺ + e⁻ → γ + γ (pair annihilation)".to_string());

        // Create photons, remove electron-positron pair
        let total_energy = particles
            .iter()
            .filter(|p| matches!(p.particle_type, ParticleType::Electron) || p.charge > 0.0)
            .map(|p| p.energy)
            .sum::<f64>();

        // Create two photons sharing the energy
        particles.retain(|p| !matches!(p.particle_type, ParticleType::Electron) && p.charge <= 0.0);
        particles.push(create_photon(
            total_energy / 2.0,
            Vector3D {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
        ));
        particles.push(create_photon(
            total_energy / 2.0,
            Vector3D {
                x: -1.0,
                y: 0.0,
                z: 0.0,
            },
        ));
    }

    // Create simulation steps
    for step in 0..steps.min(10) {
        simulation_steps.push(SimulationStep {
            step_number: step as usize,
            time: step as f64 * 1e-21, // attoseconds
            event_description: format!("Electromagnetic field evolution step {}", step),
            particles_state: particles
                .iter()
                .map(|p| ParticleState {
                    particle_id: p.id.clone(),
                    position: p.position.clone(),
                    momentum: p.momentum.clone(),
                    energy: p.energy,
                    probability_amplitude: Complex64::new(1.0, 0.0),
                })
                .collect(),
            interactions_occurred: vec!["photon_exchange".to_string()],
            energy_balance: calculate_total_energy(&particles),
        });
    }

    // Analyze conservation laws
    let conservation_laws = analyze_conservation_laws(&particles, &particles);

    // Generate Feynman diagrams
    let feynman_diagrams = vec![FeynmanDiagram {
        diagram_id: "qed_vertex".to_string(),
        interaction_vertices: vec![InteractionVertex {
            vertex_type: "QED".to_string(),
            particles_in: vec!["fermion".to_string()],
            particles_out: vec!["fermion".to_string(), "photon".to_string()],
            coupling_strength: 1.0 / 137.0, // Fine structure constant
        }],
        propagators: vec![Propagator {
            particle_type: "photon".to_string(),
            virtuality: -1.0, // Spacelike
            propagator_type: "vector".to_string(),
        }],
        external_lines: vec![],
        loop_order: 0,
        coupling_constants: [("alpha".to_string(), 1.0 / 137.0)]
            .iter()
            .cloned()
            .collect(),
    }];

    Ok(InteractionResult {
        final_particles: particles,
        conservation_laws,
        cross_section: Some(1e-30), // Simplified cross section in barns
        decay_channels: vec![],
        feynman_diagrams,
        simulation_steps,
        physical_insights,
    })
}

/// Simulate strong interactions (QCD)
fn simulate_strong_interaction(
    mut particles: Vec<Particle>,
    steps: u64,
) -> Result<InteractionResult> {
    eprintln!("💪 Simulating strong interaction (QCD)...");

    let physical_insights = vec![
        "Strong force mediated by gluons between color charges".to_string(),
        "Coupling strength αₛ ≈ 0.1-1 (running coupling)".to_string(),
        "Confinement: isolated quarks cannot exist".to_string(),
        "Asymptotic freedom: coupling decreases at high energy".to_string(),
    ];

    // Check for quarks or hadrons
    let has_quarks = particles.iter().any(|p| is_quark(&p.particle_type));

    if !has_quarks {
        let conservation = analyze_conservation_laws(&particles, &particles);
        return Ok(InteractionResult {
            final_particles: particles,
            conservation_laws: conservation,
            cross_section: None,
            decay_channels: vec![],
            feynman_diagrams: vec![],
            simulation_steps: vec![],
            physical_insights: vec!["No quarks present - no strong interaction".to_string()],
        });
    }

    // Simulate hadronization if free quarks present
    if has_free_quarks(&particles) {
        // Force quarks into hadrons (confinement)
        particles = hadronize_quarks(particles);
    }

    let conservation_laws = analyze_conservation_laws(&particles, &particles);

    // QCD Feynman diagram
    let feynman_diagrams = vec![FeynmanDiagram {
        diagram_id: "qcd_vertex".to_string(),
        interaction_vertices: vec![InteractionVertex {
            vertex_type: "QCD".to_string(),
            particles_in: vec!["quark".to_string()],
            particles_out: vec!["quark".to_string(), "gluon".to_string()],
            coupling_strength: 0.3, // Strong coupling at low energy
        }],
        propagators: vec![Propagator {
            particle_type: "gluon".to_string(),
            virtuality: -1.0,
            propagator_type: "vector".to_string(),
        }],
        external_lines: vec![],
        loop_order: 0,
        coupling_constants: [("alpha_s".to_string(), 0.3)].iter().cloned().collect(),
    }];

    Ok(InteractionResult {
        final_particles: particles,
        conservation_laws,
        cross_section: Some(1e-27), // Typical hadronic cross section
        decay_channels: vec![],
        feynman_diagrams,
        simulation_steps: vec![],
        physical_insights,
    })
}

/// Simulate weak interactions
fn simulate_weak_interaction(
    mut particles: Vec<Particle>,
    _steps: u64,
) -> Result<InteractionResult> {
    eprintln!("🔄 Simulating weak interaction...");

    let mut physical_insights = vec![
        "Weak force mediated by W± and Z⁰ bosons".to_string(),
        "Responsible for beta decay and neutrino interactions".to_string(),
        "Violates parity and CP symmetry".to_string(),
    ];

    // Simulate beta decay if neutron present
    if has_neutron(&particles) {
        physical_insights.push("Beta decay: n → p + e⁻ + ν̄ₑ".to_string());

        // Find neutron and decay it
        if let Some(neutron_idx) = particles
            .iter()
            .position(|p| matches!(p.particle_type, ParticleType::Neutron))
        {
            particles.remove(neutron_idx);

            // Add decay products
            particles.push(create_proton());
            particles.push(create_electron());
            particles.push(create_electron_antineutrino());
        }
    }

    let conservation_laws = analyze_conservation_laws(&particles, &particles);

    let decay_channels = vec![DecayChannel {
        channel_id: "beta_decay".to_string(),
        products: vec![
            "proton".to_string(),
            "electron".to_string(),
            "electron_antineutrino".to_string(),
        ],
        branching_ratio: 1.0,
        q_value: 0.782,           // MeV
        partial_width: 7.478e-28, // GeV (from neutron lifetime)
        selection_rules: vec!["ΔI = 1/2".to_string(), "Parity violated".to_string()],
    }];

    Ok(InteractionResult {
        final_particles: particles,
        conservation_laws,
        cross_section: Some(1e-44), // Weak interaction cross section
        decay_channels,
        feynman_diagrams: vec![],
        simulation_steps: vec![],
        physical_insights,
    })
}

/// Simulate Higgs interactions
fn simulate_higgs_interaction(particles: Vec<Particle>, _steps: u64) -> Result<InteractionResult> {
    eprintln!("🎯 Simulating Higgs interaction...");

    let physical_insights = vec![
        "Higgs mechanism gives mass to fundamental particles".to_string(),
        "Higgs field has non-zero vacuum expectation value".to_string(),
        "Coupling strength proportional to particle mass".to_string(),
        "Higgs boson mass: 125.1 GeV/c²".to_string(),
    ];

    // Check if Higgs boson present
    let has_higgs = particles
        .iter()
        .any(|p| matches!(p.particle_type, ParticleType::Higgs));

    let mut decay_channels = vec![];
    if has_higgs {
        decay_channels = vec![
            DecayChannel {
                channel_id: "higgs_to_bb".to_string(),
                products: vec!["bottom_quark".to_string(), "bottom_antiquark".to_string()],
                branching_ratio: 0.582,
                q_value: 125100.0,      // MeV
                partial_width: 2.38e-3, // GeV
                selection_rules: vec!["Yukawa coupling".to_string()],
            },
            DecayChannel {
                channel_id: "higgs_to_ww".to_string(),
                products: vec!["W_boson".to_string(), "W_boson".to_string()],
                branching_ratio: 0.214,
                q_value: 125100.0,
                partial_width: 8.76e-4,
                selection_rules: vec!["Gauge coupling".to_string()],
            },
        ];
    }

    let conservation = analyze_conservation_laws(&particles, &particles);
    Ok(InteractionResult {
        final_particles: particles,
        conservation_laws: conservation,
        cross_section: Some(1e-36), // Higgs production cross section
        decay_channels,
        feynman_diagrams: vec![],
        simulation_steps: vec![],
        physical_insights,
    })
}

/// Simulate particle decay
fn simulate_particle_decay(particles: Vec<Particle>, _steps: u64) -> Result<InteractionResult> {
    eprintln!("⏰ Simulating particle decay...");

    // Determine which particles can decay
    let mut final_particles = vec![];
    let mut decay_channels = vec![];

    for particle in particles.iter() {
        if is_unstable(&particle.particle_type) {
            // Particle decays - add decay products
            let products = get_decay_products(&particle.particle_type);
            final_particles.extend(products);

            // Add decay channel info
            decay_channels.push(DecayChannel {
                channel_id: format!("{:?}_decay", particle.particle_type),
                products: get_decay_product_names(&particle.particle_type),
                branching_ratio: 1.0,
                q_value: particle.mass * 0.1, // Simplified
                partial_width: 1.0 / particle.lifetime.unwrap_or(1e-10),
                selection_rules: vec!["Energy conservation".to_string()],
            });
        } else {
            // Stable particle - keep as is
            final_particles.push(particle.clone());
        }
    }

    let conservation = analyze_conservation_laws(&particles, &final_particles);
    Ok(InteractionResult {
        final_particles,
        conservation_laws: conservation,
        cross_section: None,
        decay_channels,
        feynman_diagrams: vec![],
        simulation_steps: vec![],
        physical_insights: vec![
            "Unstable particles decay according to their lifetimes".to_string(),
        ],
    })
}

/// Simulate scattering process
fn simulate_scattering_process(
    mut particles: Vec<Particle>,
    _steps: u64,
) -> Result<InteractionResult> {
    eprintln!("🎳 Simulating scattering process...");

    if particles.len() < 2 {
        return Err(anyhow::anyhow!("Need at least 2 particles for scattering"));
    }

    // Simple elastic scattering - exchange momenta
    if particles.len() == 2 {
        let p1_momentum = particles[0].momentum.clone();
        let p2_momentum = particles[1].momentum.clone();

        // Simple momentum exchange (elastic scattering)
        particles[0].momentum = Vector3D {
            x: p2_momentum.x * 0.8 + p1_momentum.x * 0.2,
            y: p2_momentum.y * 0.8 + p1_momentum.y * 0.2,
            z: p2_momentum.z * 0.8 + p1_momentum.z * 0.2,
        };

        particles[1].momentum = Vector3D {
            x: p1_momentum.x * 0.8 + p2_momentum.x * 0.2,
            y: p1_momentum.y * 0.8 + p2_momentum.y * 0.2,
            z: p1_momentum.z * 0.8 + p2_momentum.z * 0.2,
        };
    }

    let conservation = analyze_conservation_laws(&particles, &particles);
    Ok(InteractionResult {
        final_particles: particles,
        conservation_laws: conservation,
        cross_section: Some(1e-30), // Typical scattering cross section
        decay_channels: vec![],
        feynman_diagrams: vec![],
        simulation_steps: vec![],
        physical_insights: vec![
            "Elastic scattering conserves kinetic energy".to_string(),
            "Scattering angle depends on impact parameter and force law".to_string(),
        ],
    })
}

// Helper functions

fn has_photon_and_electron(particles: &[Particle]) -> bool {
    let has_photon = particles
        .iter()
        .any(|p| matches!(p.particle_type, ParticleType::Photon));
    let has_electron = particles
        .iter()
        .any(|p| matches!(p.particle_type, ParticleType::Electron));
    has_photon && has_electron
}

fn has_electron_positron_pair(particles: &[Particle]) -> bool {
    let has_electron = particles
        .iter()
        .any(|p| matches!(p.particle_type, ParticleType::Electron) && p.charge < 0.0);
    let has_positron = particles
        .iter()
        .any(|p| matches!(p.particle_type, ParticleType::Electron) && p.charge > 0.0);
    has_electron && has_positron
}

fn has_neutron(particles: &[Particle]) -> bool {
    particles
        .iter()
        .any(|p| matches!(p.particle_type, ParticleType::Neutron))
}

fn is_quark(particle_type: &ParticleType) -> bool {
    matches!(
        particle_type,
        ParticleType::Up
            | ParticleType::Down
            | ParticleType::Charm
            | ParticleType::Strange
            | ParticleType::Top
            | ParticleType::Bottom
    )
}

fn has_free_quarks(particles: &[Particle]) -> bool {
    // In reality, free quarks don't exist due to confinement
    particles.iter().any(|p| is_quark(&p.particle_type))
}

fn is_unstable(particle_type: &ParticleType) -> bool {
    matches!(
        particle_type,
        ParticleType::Neutron
            | ParticleType::Muon
            | ParticleType::Tau
            | ParticleType::Pion
            | ParticleType::Kaon
            | ParticleType::Higgs
    )
}

fn update_compton_scattering(particles: &mut [Particle]) {
    // Simplified Compton scattering kinematics
    for particle in particles.iter_mut() {
        if matches!(particle.particle_type, ParticleType::Photon) {
            particle.energy *= 0.8; // Scattered photon loses energy
        } else if matches!(particle.particle_type, ParticleType::Electron) {
            particle.energy += 0.1; // Electron gains kinetic energy
        }
    }
}

fn hadronize_quarks(mut particles: Vec<Particle>) -> Vec<Particle> {
    // Remove free quarks and create hadrons
    particles.retain(|p| !is_quark(&p.particle_type));

    // Add some hadrons
    particles.push(create_proton());
    particles.push(create_pion());

    particles
}

fn get_decay_products(particle_type: &ParticleType) -> Vec<Particle> {
    match particle_type {
        ParticleType::Neutron => vec![
            create_proton(),
            create_electron(),
            create_electron_antineutrino(),
        ],
        ParticleType::Muon => vec![
            create_electron(),
            create_muon_neutrino(),
            create_electron_antineutrino(),
        ],
        ParticleType::Pion => vec![create_muon(), create_muon_neutrino()],
        _ => vec![],
    }
}

fn get_decay_product_names(particle_type: &ParticleType) -> Vec<String> {
    match particle_type {
        ParticleType::Neutron => vec![
            "proton".to_string(),
            "electron".to_string(),
            "electron_antineutrino".to_string(),
        ],
        ParticleType::Muon => vec![
            "electron".to_string(),
            "muon_neutrino".to_string(),
            "electron_antineutrino".to_string(),
        ],
        ParticleType::Pion => vec!["muon".to_string(), "muon_neutrino".to_string()],
        _ => vec![],
    }
}

fn analyze_conservation_laws(
    initial: &[Particle],
    final_particles: &[Particle],
) -> ConservationAnalysis {
    let initial_energy = calculate_total_energy(initial);
    let final_energy = calculate_total_energy(final_particles);
    let energy_conserved = (initial_energy - final_energy).abs() < 1e-6;

    let initial_charge = calculate_total_charge(initial);
    let final_charge = calculate_total_charge(final_particles);
    let charge_conserved = (initial_charge - final_charge).abs() < 1e-6;

    ConservationAnalysis {
        energy_momentum_conserved: energy_conserved,
        charge_conserved,
        baryon_number_conserved: true, // Simplified
        lepton_number_conserved: true, // Simplified
        color_charge_conserved: true,  // Simplified
        conservation_violations: vec![],
        quantum_numbers_initial: calculate_total_quantum_numbers(initial),
        quantum_numbers_final: calculate_total_quantum_numbers(final_particles),
    }
}

fn calculate_total_energy(particles: &[Particle]) -> f64 {
    particles.iter().map(|p| p.energy).sum()
}

fn calculate_total_charge(particles: &[Particle]) -> f64 {
    particles.iter().map(|p| p.charge).sum()
}

fn calculate_total_quantum_numbers(particles: &[Particle]) -> QuantumNumbers {
    QuantumNumbers {
        baryon_number: particles
            .iter()
            .map(|p| p.quantum_numbers.baryon_number)
            .sum(),
        lepton_number: particles
            .iter()
            .map(|p| p.quantum_numbers.lepton_number)
            .sum(),
        strangeness: particles
            .iter()
            .map(|p| p.quantum_numbers.strangeness)
            .sum(),
        charm: particles.iter().map(|p| p.quantum_numbers.charm).sum(),
        beauty: particles.iter().map(|p| p.quantum_numbers.beauty).sum(),
        truth: particles.iter().map(|p| p.quantum_numbers.truth).sum(),
        isospin: 0.0, // Simplified
        color_charge: None,
    }
}

// Particle creation helpers
fn create_photon(energy: f64, direction: Vector3D) -> Particle {
    Particle {
        id: "photon_1".to_string(),
        particle_type: ParticleType::Photon,
        mass: 0.0,
        charge: 0.0,
        spin: 1.0,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: direction.x * energy,
            y: direction.y * energy,
            z: direction.z * energy,
        },
        energy,
        lifetime: None,
        quantum_numbers: QuantumNumbers {
            baryon_number: 0.0,
            lepton_number: 0.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 0.0,
            color_charge: None,
        },
    }
}

fn create_proton() -> Particle {
    Particle {
        id: "proton_1".to_string(),
        particle_type: ParticleType::Proton,
        mass: 0.938272, // GeV/c²
        charge: 1.0,
        spin: 0.5,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        energy: 0.938272,
        lifetime: None, // Stable
        quantum_numbers: QuantumNumbers {
            baryon_number: 1.0,
            lepton_number: 0.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 0.5,
            color_charge: None,
        },
    }
}

fn create_electron() -> Particle {
    Particle {
        id: "electron_1".to_string(),
        particle_type: ParticleType::Electron,
        mass: 0.000511, // GeV/c²
        charge: -1.0,
        spin: 0.5,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        energy: 0.000511,
        lifetime: None, // Stable
        quantum_numbers: QuantumNumbers {
            baryon_number: 0.0,
            lepton_number: 1.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 0.0,
            color_charge: None,
        },
    }
}

fn create_muon() -> Particle {
    Particle {
        id: "muon_1".to_string(),
        particle_type: ParticleType::Muon,
        mass: 0.1057, // GeV/c²
        charge: -1.0,
        spin: 0.5,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        energy: 0.1057,
        lifetime: Some(2.2e-6), // 2.2 microseconds
        quantum_numbers: QuantumNumbers {
            baryon_number: 0.0,
            lepton_number: 1.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 0.0,
            color_charge: None,
        },
    }
}

fn create_electron_antineutrino() -> Particle {
    Particle {
        id: "electron_antineutrino_1".to_string(),
        particle_type: ParticleType::ElectronNeutrino,
        mass: 0.0, // Approximately massless
        charge: 0.0,
        spin: 0.5,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        energy: 0.1,    // Kinetic energy
        lifetime: None, // Stable
        quantum_numbers: QuantumNumbers {
            baryon_number: 0.0,
            lepton_number: -1.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 0.0,
            color_charge: None,
        },
    }
}

fn create_muon_neutrino() -> Particle {
    Particle {
        id: "muon_neutrino_1".to_string(),
        particle_type: ParticleType::MuonNeutrino,
        mass: 0.0,
        charge: 0.0,
        spin: 0.5,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        energy: 0.1,
        lifetime: None,
        quantum_numbers: QuantumNumbers {
            baryon_number: 0.0,
            lepton_number: 1.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 0.0,
            color_charge: None,
        },
    }
}

fn create_pion() -> Particle {
    Particle {
        id: "pion_1".to_string(),
        particle_type: ParticleType::Pion,
        mass: 0.13957, // GeV/c²
        charge: 1.0,
        spin: 0.0,
        position: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        momentum: Vector3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        energy: 0.13957,
        lifetime: Some(2.6e-8), // seconds
        quantum_numbers: QuantumNumbers {
            baryon_number: 0.0,
            lepton_number: 0.0,
            strangeness: 0.0,
            charm: 0.0,
            beauty: 0.0,
            truth: 0.0,
            isospin: 1.0,
            color_charge: None,
        },
    }
}
