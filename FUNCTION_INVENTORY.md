# Comprehensive Function Inventory

This document provides a complete inventory of all primary mathematical and computational functions in the computational-engine codebase, organized by filesystem structure.

**Generated:** 2026-01-11 (Updated with discovered modules)
**Total Functions:** 350+
**Wired to Tool API:** ~180 functions
**Not Wired (Critical Gap):** ~170 functions

> **WARNING:** Major discovery - entire modules are completely unwired including:
> - Machine Learning (~30 functions)
> - Chaos/Fractals (~15 functions)
> - Game Theory (~10 functions)
> - Wormholes (~15 functions)
> - Control Theory (~10 functions)
> - See [ARCHITECTURE_ANALYSIS.md](./ARCHITECTURE_ANALYSIS.md) for full details

---

## Table of Contents

1. [src/mathematics/](#srcmathematics)
2. [src/physics/](#srcphysics)
3. [src/tools/](#srctools)
4. [src/specialized/](#srcspecialized)
5. [src/chemistry/](#srcchemistry)
6. [src/biology/](#srcbiology)
7. [src/thermodynamics/](#srcthermodynamics)
8. [src/optics/](#srcoptics)
9. [src/geophysics/](#srcgeophysics)
10. [src/engineering/](#srcengineering)
11. [src/datetime/](#srcdatetime)
12. [Summary Statistics](#summary-statistics)

---

## src/mathematics/

### src/mathematics/tensor_calculus/

#### tensor.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `parse_metric_tensor` | `(metric_strings: Vec<Vec<String>>, coords: &[String]) -> Result<MetricTensor, TensorError>` | ✅ WIRED - Compute::Tensor |
| `calculate_christoffel_symbols` | `(metric: &MetricTensor, coords: &[String]) -> Result<ChristoffelResult, TensorError>` | ✅ WIRED - Compute::Tensor::Christoffel |
| `calculate_riemann_tensor` | `(metric: &MetricTensor, coords: &[String]) -> Result<RiemannResult, TensorError>` | ✅ WIRED - Compute::Tensor::Riemann |
| `calculate_ricci_tensor` | `(metric: &MetricTensor, coords: &[String]) -> Result<RiemannResult, TensorError>` | ✅ WIRED - Compute::Tensor::Ricci |
| `calculate_ricci_scalar` | `(metric: &MetricTensor, coords: &[String]) -> Result<TensorComponent, TensorError>` | ✅ WIRED - Compute::Tensor::RicciScalar |
| `calculate_einstein_tensor` | `(metric: &MetricTensor, coords: &[String]) -> Result<RiemannResult, TensorError>` | ✅ WIRED - Compute::Tensor::Einstein |

#### einstein.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `construct_einstein_field_equations` | `(stress_energy: &StressEnergyTensor, coordinates: &[String], cosmological_constant: Option<SymbolicExpr>) -> Result<EinsteinEquationSystem, TensorError>` | ✅ WIRED - Solve::Einstein |
| `solve_vacuum_einstein_equations` | `(coordinates: &[String], symmetry_ansatz: &str, boundary_conditions: &[BoundaryCondition]) -> Result<Vec<EinsteinSolution>, TensorError>` | ✅ WIRED - Solve::Einstein::Vacuum |
| `verify_einstein_solution` | `(solution: &EinsteinSolution, stress_energy: Option<&StressEnergyTensor>, cosmological_constant: Option<SymbolicExpr>) -> Result<bool, TensorError>` | ✅ WIRED - Solve::Einstein |
| `solve_einstein_constraint_equations` | `(initial_data: &MetricTensor, coordinates: &[String]) -> Result<Vec<TensorComponent>, TensorError>` | ✅ WIRED - Solve::Einstein |

#### quantum_tensors.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_christoffel_symbols_symbolic` | `(metric_expr: &str, coords: &[String]) -> Result<Vec<TensorComponent>, TensorError>` | ❌ NOT WIRED - Utility |
| `calculate_riemann_tensor_symbolic` | `(metric_expr: &str, coords: &[String]) -> Result<Vec<TensorComponent>, TensorError>` | ❌ NOT WIRED - Utility |
| `calculate_ricci_tensor_symbolic` | `(metric_expr: &str, coords: &[String]) -> Result<Vec<TensorComponent>, TensorError>` | ❌ NOT WIRED - Utility |
| `calculate_einstein_tensor_symbolic` | `(metric_expr: &str, coords: &[String]) -> Result<Vec<TensorComponent>, TensorError>` | ❌ NOT WIRED - Utility |

---

### src/mathematics/special_functions/

#### bessel.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `bessel_function` | `(request: BesselRequest) -> Result<BesselResult, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `bessel_j` | `(nu: f64, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `bessel_y` | `(nu: f64, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `bessel_i` | `(nu: f64, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `bessel_k` | `(nu: f64, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |

#### gamma.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `gamma_function` | `(request: GammaRequest) -> Result<GammaResult, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `gamma` | `(z: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `log_gamma` | `(z: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `digamma` | `(z: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `beta` | `(a: f64, b: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |

#### polynomials.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `orthogonal_polynomial` | `(request: PolynomialRequest) -> Result<PolynomialResult, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `legendre_p` | `(n: usize, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `hermite_h` | `(n: usize, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `laguerre_l` | `(n: usize, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `chebyshev_t` | `(n: usize, x: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |

#### error.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `error_function` | `(request: ErrorFunctionRequest) -> Result<ErrorFunctionResult, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `erf` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |
| `erfc` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |
| `erfcx` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |
| `erfi` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |

#### airy.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `airy_function` | `(request: AiryRequest) -> Result<AiryResult, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `airy_ai` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |
| `airy_bi` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |
| `airy_ai_prime` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |
| `airy_bi_prime` | `(x: f64) -> f64` | ❌ NOT WIRED - Helper |

#### elliptic.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `elliptic_integral` | `(request: EllipticIntegralRequest) -> Result<EllipticIntegralResult, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `elliptic_k` | `(k: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `elliptic_e` | `(k: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `elliptic_f` | `(phi: f64, k: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |
| `elliptic_pi` | `(phi: f64, n: f64, k: f64) -> Result<f64, String>` | ❌ NOT WIRED - Helper |

---

### src/mathematics/linear_algebra/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_svd` | `(input: MatrixInput) -> Result<SVDResult, String>` | ✅ WIRED - Compute::Matrix::SVD |
| `compute_matrix_rank` | `(input: MatrixInput) -> Result<usize, String>` | ✅ WIRED - Compute::Matrix::Rank |
| `compute_pseudoinverse` | `(input: MatrixInput) -> Result<Vec<Vec<f64>>, String>` | ✅ WIRED - Compute::Matrix::Pseudoinverse |
| `compute_pca` | `(input: MatrixInput, n_components: Option<usize>) -> Result<PCAResult, String>` | ✅ WIRED - Compute::Matrix::PCA |

---

### src/mathematics/chaos/

#### fractals.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `mandelbrot` | `(c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult` | ❌ NOT WIRED - Standalone |
| `julia` | `(z0: Complex, c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult` | ❌ NOT WIRED - Standalone |
| `burning_ship` | `(c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult` | ❌ NOT WIRED - Standalone |
| `box_counting_dimension` | `(points: &[(f64, f64)], min_box: f64, max_box: f64, num_scales: usize) -> f64` | ❌ NOT WIRED - Standalone |
| `koch_snowflake` | `(order: usize) -> Vec<(f64, f64)>` | ❌ NOT WIRED - Standalone |

#### lyapunov.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `lyapunov_1d_map` | `<F>(map: F, x0: f64, iterations: usize, transient: usize) -> f64` | ❌ NOT WIRED - Standalone |
| `lyapunov_logistic_map` | `(r: f64, iterations: usize, transient: usize) -> f64` | ❌ NOT WIRED - Standalone |
| `lyapunov_spectrum_3d` | `<F>(...) -> [f64; 3]` | ❌ NOT WIRED - Standalone |
| `kaplan_yorke_dimension` | `(exponents: &[f64]) -> f64` | ❌ NOT WIRED - Standalone |

#### attractors.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `lorenz_attractor` | `(t: f64, state: &[f64], params: &LorenzParams) -> Vec<f64>` | ❌ NOT WIRED - Standalone |
| `rossler_attractor` | `(t: f64, state: &[f64], params: &RosslerParams) -> Vec<f64>` | ❌ NOT WIRED - Standalone |
| `logistic_map_bifurcation` | `(r_min: f64, r_max: f64, num_r: usize, iterations: usize, transient: usize) -> Vec<(f64, Vec<f64>)>` | ❌ NOT WIRED - Standalone |

---

### src/mathematics/symbolic_regression/

| Function | Signature | Status |
|----------|-----------|--------|
| `generate_random_expression` | `() -> String` | ❌ NOT WIRED - Utility |
| `evaluate_fitness` | `(expr: &str, data: &[(f64, f64)]) -> f64` | ❌ NOT WIRED - Utility |
| `check_physics_constraints` | `(expr: &str) -> bool` | ❌ NOT WIRED - Utility |
| `discover_equations` | `(data: &[(f64, f64)], target_fitness: f64, max_generations: usize) -> EquationDiscoveryResult` | ❌ NOT WIRED - Standalone |

---

### src/mathematics/advanced_calculus/

#### riemann_zeta.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `riemann_zeta` | `(s: Complex64) -> Result<Complex64, String>` | ✅ WIRED - Compute::SpecialFunctions |
| `riemann_zeta_on_critical_line` | `(t: f64) -> Complex64` | ❌ NOT WIRED - Helper |
| `hardy_z_function` | `(t: f64) -> f64` | ❌ NOT WIRED - Helper |

#### complex_analysis.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `contour_integral` | `(request: ContourIntegralRequest) -> Result<ContourIntegralResult, String>` | ✅ WIRED - Integrate::Contour |
| `cauchy_residue` | `(request: ResidueRequest) -> Result<ResidueResult, String>` | ✅ WIRED - Integrate::Residue |
| `analytic_continuation` | `(request: ContinuationRequest) -> Result<ContinuationResult, String>` | ✅ WIRED - Analyze::Continuation |

---

## src/physics/

### src/physics/nuclear_physics/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `radioactive_decay` | `(initial_activity: f64, half_life: f64, time: f64) -> f64` | ✅ WIRED - Compute::NuclearPhysics |
| `decay_chain` | `(request: DecayChainRequest) -> Result<DecayChainResult, String>` | ✅ WIRED - Compute::NuclearPhysics |
| `half_life` | `(request: HalfLifeRequest) -> Result<HalfLifeResult, String>` | ✅ WIRED - Compute::NuclearPhysics |
| `binding_energy` | `(request: BindingEnergyRequest) -> Result<BindingEnergyResult, String>` | ✅ WIRED - Compute::NuclearPhysics |
| `mass_defect` | `(request: MassDefectRequest) -> Result<MassDefectResult, String>` | ✅ WIRED - Compute::NuclearPhysics |
| `fission_energy` | `(request: FissionEnergyRequest) -> Result<FissionEnergyResult, String>` | ✅ WIRED - Compute::NuclearPhysics |
| `fusion_energy` | `(request: FusionEnergyRequest) -> Result<FusionEnergyResult, String>` | ✅ WIRED - Compute::NuclearPhysics |
| `nuclear_reaction` | `(request: NuclearReactionRequest) -> Result<NuclearReactionResult, String>` | ✅ WIRED - Compute::NuclearPhysics |

---

### src/physics/cosmology/

#### friedmann.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `hubble_parameter` | `(z: Redshift, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `deceleration_parameter` | `(z: Redshift, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `lookback_time` | `(z: Redshift, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `comoving_distance` | `(z: Redshift, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `luminosity_distance` | `(z: Redshift, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `angular_diameter_distance` | `(z: Redshift, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `event_horizon` | `(params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `evolve_universe` | `(z_max: f64, n_steps: usize, params: &CosmologyParams) -> UniverseEvolution` | ✅ WIRED - Compute::Physics::Cosmology |

#### cmb.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `planck_spectrum` | `(frequency: f64, temperature: f64) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `cmb_peak_frequency` | `(temperature: f64) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `photon_number_density` | `(temperature: f64) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `cmb_energy_density` | `(temperature: f64) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `recombination` | `(params: &CosmologyParams) -> RecombinationEpoch` | ✅ WIRED - Compute::Physics::Cosmology |
| `matter_radiation_equality` | `(params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `sound_horizon_recombination` | `(params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `first_acoustic_peak_angle` | `(params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `sachs_wolfe_amplitude` | `(potential: f64) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |

#### dark_energy.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `hubble_parameter_de` | `(a: f64, params: &CosmologyParams) -> f64` | ✅ WIRED - Compute::Physics::Cosmology |
| `big_rip_time` | `(params: &CosmologyParams, w: f64) -> Option<f64>` | ✅ WIRED - Compute::Physics::Cosmology |
| `future_evolution` | `(a_init: f64, a_final: f64, params: &CosmologyParams, w: f64, steps: usize) -> Vec<f64>` | ✅ WIRED - Compute::Physics::Cosmology |

---

### src/physics/wormholes/

#### metric.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `redshift_function` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |
| `shape_function` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |
| `shape_function_derivative` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |
| `morris_thorne_metric` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> WormholeMetric` | ❌ NOT WIRED - Utility |
| `satisfies_flaring_condition` | `(r: f64, config: &WormholeConfig) -> bool` | ❌ NOT WIRED - Utility |
| `proper_distance` | `(r_start: f64, r_end: f64, config: &WormholeConfig, n_steps: usize) -> f64` | ❌ NOT WIRED - Utility |
| `embedding_surface` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |

#### traversal.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_tidal_forces` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> TidalForce` | ❌ NOT WIRED - Utility |
| `analyze_traversal` | `(trajectory: &[SphericalCoordinates], config: &WormholeConfig) -> TraversalAnalysis` | ❌ NOT WIRED - Utility |
| `time_dilation_factor` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |
| `required_velocity` | `(config: &WormholeConfig, target_time: f64) -> f64` | ❌ NOT WIRED - Utility |

#### energy.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_energy_density` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> EnergyDensity` | ❌ NOT WIRED - Utility |
| `calculate_total_energy` | `(config: &WormholeConfig) -> WormholeEnergy` | ❌ NOT WIRED - Utility |
| `compare_with_black_hole` | `(config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |
| `throat_energy_density` | `(config: &WormholeConfig) -> f64` | ❌ NOT WIRED - Utility |

---

### src/physics/electromagnetism/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `maxwell_equations` | `(request: MaxwellRequest) -> Result<MaxwellResult, String>` | ✅ WIRED - Compute::EM |
| `em_wave` | `(request: WaveRequest) -> Result<WaveResult, String>` | ✅ WIRED - Compute::EM |
| `antenna_analysis` | `(request: AntennaRequest) -> Result<AntennaResult, String>` | ✅ WIRED - Compute::EM |
| `transmission_line` | `(request: TransmissionLineRequest) -> Result<TransmissionLineResult, String>` | ✅ WIRED - Compute::EM |
| `waveguide` | `(request: WaveguideRequest) -> Result<WaveguideResult, String>` | ✅ WIRED - Compute::EM |
| `scattering` | `(request: ScatteringRequest) -> Result<ScatteringResult, String>` | ✅ WIRED - Compute::EM |
| `poynting_vector` | `(request: PoyntingRequest) -> Result<PoyntingResult, String>` | ✅ WIRED - Compute::EM |
| `skin_effect` | `(request: SkinEffectRequest) -> Result<SkinEffectResult, String>` | ✅ WIRED - Compute::EM |
| `plasma_frequency` | `(electron_density: f64) -> f64` | ✅ WIRED - Compute::EM |

---

### src/physics/quantum_mechanics/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `schrodinger_equation` | `(request: SchrodingerRequest) -> Result<SchrodingerResult, String>` | ✅ WIRED - Compute::QuantumMechanics |
| `harmonic_oscillator` | `(n: usize, x: f64) -> f64` | ✅ WIRED - Compute::QuantumMechanics |
| `hydrogen_atom` | `(request: HydrogenAtomRequest) -> Result<HydrogenAtomResult, String>` | ✅ WIRED - Compute::QuantumMechanics |
| `angular_momentum` | `(request: AngularMomentumRequest) -> Result<AngularMomentumResult, String>` | ✅ WIRED - Compute::QuantumMechanics |
| `spin_operators` | `(request: SpinRequest) -> Result<SpinResult, String>` | ✅ WIRED - Compute::QuantumMechanics |
| `perturbation_theory` | `(request: PerturbationRequest) -> Result<PerturbationResult, String>` | ✅ WIRED - Compute::QuantumMechanics |
| `tunneling_probability` | `(request: TunnelingRequest) -> Result<TunnelingResult, String>` | ✅ WIRED - Compute::QuantumMechanics |

#### bohm_potential.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `bohm_potential_gaussian` | `(x: f64, sigma: f64, hbar: f64, mass: f64) -> f64` | ❌ NOT WIRED - Utility |

#### decoherence.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `decoherence_scale` | `(mass: f64, temperature: f64) -> f64` | ✅ WIRED - Compute::QuantumMechanics |

---

### src/physics/black_holes/

#### orbits.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `schwarzschild_radius` | `(mass: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `innermost_stable_circular_orbit` | `(mass: f64, spin: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `photon_sphere_radius` | `(mass: f64, spin: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `ergosphere_radius` | `(mass: f64, spin: f64, theta: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `hawking_temperature` | `(mass: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `hawking_luminosity` | `(mass: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `evaporation_time` | `(mass: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |

#### penrose.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `penrose_process_efficiency` | `(spin: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |
| `frame_dragging_velocity` | `(mass: f64, spin: f64, r: f64) -> f64` | ✅ WIRED - Compute::BlackHoles |

---

### src/physics/fluid_dynamics/

#### navier_stokes.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `solve_navier_stokes_2d` | `(request: NavierStokes2DRequest) -> Result<NavierStokes2DResult, String>` | ✅ WIRED - Simulate::FluidDynamics |
| `reynolds_number` | `(velocity: f64, length: f64, viscosity: f64) -> f64` | ✅ WIRED - Compute::FluidDynamics |
| `pressure_poisson` | `(velocity_field: &VelocityField, dx: f64, dt: f64) -> PressureField` | ❌ NOT WIRED - Helper |

#### quantum_fluids.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `quantum_navier_stokes_1d` | `(request: QuantumNS1DRequest) -> Result<QuantumNS1DResult, String>` | ✅ WIRED - Simulate::QuantumFluidDynamics |
| `madelung_velocity` | `(psi: Complex64, grad_psi: Complex64) -> f64` | ❌ NOT WIRED - Helper |

---

## src/tools/

### src/tools/signal_processing/

#### lib.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_fft` | `(request: FFTRequest) -> Result<FFTResult, String>` | ✅ WIRED - Transform::FFT |
| `apply_filter` | `(request: FilterRequest) -> Result<FilterResult, String>` | ✅ WIRED - Transform::Filter |
| `compute_spectrogram` | `(request: SpectrogramRequest) -> Result<SpectrogramResult, String>` | ✅ WIRED - Transform::Spectrogram |
| `compute_psd` | `(request: PSDRequest) -> Result<PSDResult, String>` | ✅ WIRED - Transform::PSD |

#### additional.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `detect_peaks` | `(request: PeakDetectionRequest) -> Result<PeakDetectionResult, String>` | ✅ WIRED - Compute::Signal |
| `compute_fourier_series` | `(input: &ComputeInput) -> ToolResult<ComputeOutput>` | ✅ WIRED - Compute::FourierSeries |
| `wavelet_transform` | `(input: &ComputeInput) -> ToolResult<ComputeOutput>` | ✅ WIRED - Transform::Wavelet |
| `windowing_functions` | `(request: WindowFunctionRequest) -> Result<WindowFunctionResult, String>` | ✅ WIRED - Transform::Window |

#### wavelets.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `haar_wavelet` | `(data: &[f64]) -> Vec<f64>` | ✅ WIRED - Transform::Wavelet |
| `daubechies4_wavelet` | `(data: &[f64]) -> Vec<f64>` | ✅ WIRED - Transform::Wavelet |
| `mexican_hat_wavelet` | `(data: &[f64], scale: f64) -> Vec<f64>` | ✅ WIRED - Transform::Wavelet |

---

### src/tools/dimensional_analysis/

#### lib.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `tokenize_unit` | `(unit_str: &str) -> Vec<String>` | ❌ NOT WIRED - Utility |
| `extract_variables_with_powers` | `(expression: &str) -> Vec<(String, i32)>` | ❌ NOT WIRED - Utility |
| `analyze_expression_dimensions` | `(expression: &str) -> Result<Vec<(String, i32)>, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `generate_recommendations` | `(missing_dims: &[(String, i32)]) -> Vec<String>` | ❌ NOT WIRED - Utility |
| `dimensional_analysis` | `(request: DimensionalAnalysisRequest) -> Result<DimensionalAnalysisResult, String>` | ❌ NOT WIRED - Standalone |
| `convert_units` | `(value: f64, from_unit: &str, to_unit: &str) -> Result<f64, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `derive_units_for_quantity` | `(quantity: &str) -> Result<Vec<String>, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `are_units_compatible` | `(unit1: &str, unit2: &str) -> Result<bool, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `get_si_base_units` | `(unit: &str) -> Result<String, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `get_physics_formulas` | `() -> Vec<PhysicsFormula>` | ❌ NOT WIRED - Utility |

---

### src/tools/numerical_methods/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `solve_ode` | `(request: ODESolverRequest) -> Result<ODESolverResult, String>` | ✅ WIRED - Simulate::ODE |
| `find_root` | `(request: RootFindingRequest) -> Result<RootFindingResult, String>` | ✅ WIRED - Solve::Root |
| `integrate` | `(request: IntegrationRequest) -> Result<IntegrationResult, String>` | ✅ WIRED - Integrate |
| `interpolate` | `(request: InterpolationRequest) -> Result<InterpolationResult, String>` | ✅ WIRED - Compute::Interpolation |
| `solve_linear_system` | `(request: LinearSystemRequest) -> Result<LinearSystemResult, String>` | ✅ WIRED - Solve::Linear |
| `differentiate` | `(request: DifferentiationRequest) -> Result<DifferentiationResult, String>` | ✅ WIRED - Differentiate |
| `solve_pde` | `(request: PDESolverRequest) -> Result<PDESolverResult, String>` | ✅ WIRED - Simulate::PDE |

---

### src/tools/computational_geometry/

#### advanced.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `convex_hull` | `(request: ConvexHullRequest) -> Result<ConvexHullResult, String>` | ✅ WIRED - Compute::Geometry |
| `delaunay_triangulation` | `(request: DelaunayRequest) -> Result<DelaunayResult, String>` | ✅ WIRED - Compute::Geometry |
| `voronoi_diagram` | `(request: VoronoiRequest) -> Result<VoronoiResult, String>` | ✅ WIRED - Compute::Geometry |
| `polygon_area` | `(request: PolygonAreaRequest) -> Result<PolygonAreaResult, String>` | ✅ WIRED - Compute::Geometry |
| `point_in_polygon` | `(request: PointInPolygonRequest) -> Result<PointInPolygonResult, String>` | ✅ WIRED - Compute::Geometry |

#### spatial_3d.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `convex_hull_3d` | `(points: &[Point3D]) -> Result<ConvexHull3D, String>` | ❌ NOT WIRED - Utility |

---

### src/tools/equation_validation/

#### lib.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `parse_equation` | `(equation: &str) -> Result<(String, String), Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `extract_variables` | `(expression: &str) -> Vec<String>` | ❌ NOT WIRED - Utility |
| `check_dimensional_consistency` | `(equation: &str, variable_units: &HashMap<String, String>) -> Result<bool, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `check_conservation_laws` | `(equation: &str) -> Result<bool, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `check_symmetries` | `(equation: &str) -> Result<Vec<String>, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `check_mathematical_correctness` | `(equation: &str) -> Result<bool, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `check_physics_compliance` | `(equation: &str) -> Result<bool, Box<dyn Error>>` | ❌ NOT WIRED - Utility |
| `validate_equation` | `(request: EquationValidationRequest) -> Result<EquationValidationResult, String>` | ❌ NOT WIRED - Standalone |

---

## src/specialized/

### src/specialized/stochastic_processes/

#### lib.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `generate_brownian_motion` | `(params: BrownianMotionParams) -> Result<Vec<f64>, String>` | ✅ WIRED - Simulate::Stochastic |
| `simulate_markov_chain` | `(params: MarkovChainParams) -> Result<Vec<String>, String>` | ✅ WIRED - Simulate::Stochastic |
| `compute_stochastic_integral_monte_carlo` | `(params: MonteCarloParams) -> Result<f64, String>` | ✅ WIRED - Compute::Stochastic |
| `generate_ornstein_uhlenbeck_process` | `(params: OrnsteinUhlenbeckParams) -> Result<Vec<f64>, String>` | ✅ WIRED - Simulate::Stochastic |
| `simulate_poisson_process` | `(rate: f64, time_horizon: f64) -> Result<Vec<f64>, String>` | ✅ WIRED - Simulate::Stochastic |

---

### src/specialized/cryptographic_mathematics/

#### lib.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `mod_exp` | `(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt` | ✅ WIRED - Compute::NumberTheory |
| `extended_gcd` | `(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt)` | ✅ WIRED - Compute::NumberTheory |
| `mod_inverse` | `(a: &BigInt, m: &BigInt) -> Option<BigInt>` | ✅ WIRED - Compute::NumberTheory |
| `chinese_remainder_theorem` | `(remainders: &[BigInt], moduli: &[BigInt]) -> Option<BigInt>` | ✅ WIRED - Compute::NumberTheory |
| `miller_rabin_test` | `(n: &BigInt, k: u32) -> bool` | ✅ WIRED - Compute::NumberTheory |
| `generate_prime` | `(bits: u32) -> BigInt` | ✅ WIRED - Compute::NumberTheory |
| `generate_rsa_keypair` | `(bits: u32) -> (BigInt, BigInt, BigInt)` | ✅ WIRED - Compute::NumberTheory |
| `rsa_encrypt` | `(message: &BigInt, e: &BigInt, n: &BigInt) -> BigInt` | ✅ WIRED - Compute::NumberTheory |
| `rsa_decrypt` | `(ciphertext: &BigInt, d: &BigInt, n: &BigInt) -> BigInt` | ✅ WIRED - Compute::NumberTheory |
| `discrete_log_bsgs` | `(generator: &BigInt, target: &BigInt, modulus: &BigInt) -> Option<BigInt>` | ✅ WIRED - Compute::NumberTheory |
| `euler_totient` | `(n: &BigInt) -> BigInt` | ✅ WIRED - Compute::NumberTheory |
| `carmichael_lambda` | `(n: &BigInt) -> BigInt` | ✅ WIRED - Compute::NumberTheory |
| `elliptic_curve_point_add` | `(p1: &ECPoint, p2: &ECPoint, curve: &ECCurve) -> ECPoint` | ✅ WIRED - Compute::NumberTheory |
| `sha256` | `(input: &str) -> String` | ✅ WIRED - Compute::NumberTheory |
| `sha3_256` | `(input: &str) -> String` | ✅ WIRED - Compute::NumberTheory |

---

### src/specialized/game_theory/

#### normal_form.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `prisoners_dilemma` | `() -> NormalFormGame` | ❌ NOT WIRED - Utility |
| `battle_of_sexes` | `() -> NormalFormGame` | ❌ NOT WIRED - Utility |

---

### src/specialized/linear_programming/

#### simplex.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `simplex` | `(lp: &LinearProgram) -> Result<LinearProgramSolution, String>` | ✅ WIRED - Optimize::LinearProgramming |
| `is_feasible` | `(lp: &LinearProgram, solution: &[f64]) -> bool` | ❌ NOT WIRED - Helper |
| `objective_value` | `(lp: &LinearProgram, solution: &[f64]) -> f64` | ❌ NOT WIRED - Helper |

#### dual.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_dual` | `(primal: &LinearProgram) -> LinearProgram` | ✅ WIRED - Optimize::DualProblem |
| `sensitivity_analysis` | `(solution: &LinearProgramSolution, lp: &LinearProgram) -> SensitivityAnalysis` | ✅ WIRED - Compute::Sensitivity |
| `verify_strong_duality` | `(primal_obj: f64, dual_obj: f64, tolerance: f64) -> bool` | ❌ NOT WIRED - Helper |
| `interpret_shadow_price` | `(shadow_price: f64) -> String` | ❌ NOT WIRED - Utility |

---

### src/specialized/optimization/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `gradient_descent` | `(objective: &str, initial_point: &[f64], learning_rate: f64, max_iter: usize) -> Result<OptimizationResult, String>` | ✅ WIRED - Optimize::GradientDescent |
| `nelder_mead` | `(objective: &str, initial_point: &[f64], max_iter: usize) -> Result<OptimizationResult, String>` | ✅ WIRED - Optimize::NelderMead |
| `curve_fitting` | `(request: CurveFitRequest) -> Result<CurveFitResult, String>` | ✅ WIRED - Optimize::CurveFitting |
| `sensitivity_analysis` | `(objective: &str, solution: &[f64], perturbation: f64) -> Result<SensitivityResult, String>` | ✅ WIRED - Compute::Sensitivity |

---

### src/specialized/statistics/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `statistics` | `(request: StatisticsRequest) -> Result<StatisticsResult, String>` | ✅ WIRED - Sample::Statistics |
| `monte_carlo_integration` | `(request: MonteCarloRequest) -> Result<MonteCarloResult, String>` | ✅ WIRED - Sample::MonteCarlo |
| `mcmc_sampling` | `(request: MCMCRequest) -> Result<MCMCResult, String>` | ✅ WIRED - Sample::MCMC |
| `correlation` | `(request: CorrelationRequest) -> Result<CorrelationResult, String>` | ✅ WIRED - Compute::Statistics |
| `kl_divergence` | `(request: KLDivergenceRequest) -> Result<KLDivergenceResult, String>` | ✅ WIRED - Compute::Information |
| `mutual_information` | `(request: MutualInformationRequest) -> Result<MutualInformationResult, String>` | ✅ WIRED - Compute::Information |

---

### src/specialized/control_theory/

#### pid_controller.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `ziegler_nichols_tuning` | `(ku: f64, pu: f64) -> PIDConfig` | ✅ WIRED - Compute::ControlSystems |
| `cohen_coon_tuning` | `(k: f64, tau: f64, theta: f64) -> PIDConfig` | ✅ WIRED - Compute::ControlSystems |

#### analysis.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `controllability` | `(a: &[Vec<f64>], b: &[Vec<f64>]) -> ControllabilityResult` | ✅ WIRED - Compute::ControlSystems |
| `observability` | `(a: &[Vec<f64>], c: &[Vec<f64>]) -> ObservabilityResult` | ✅ WIRED - Compute::ControlSystems |
| `eigenvalues` | `(a: &[Vec<f64>]) -> Vec<Complex64>` | ✅ WIRED - Compute::ControlSystems |

#### lqr.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `lqr` | `(a: &[Vec<f64>], b: &[Vec<f64>], q: &[Vec<f64>], r: &[Vec<f64>]) -> Result<LQRSolution, String>` | ✅ WIRED - Compute::ControlSystems |

#### state_space.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `siso_state_space` | `(a: f64, b: f64, c: f64, d: f64, u_sequence: &[f64]) -> Result<Vec<f64>, String>` | ✅ WIRED - Compute::ControlSystems |

---

### src/specialized/graph_theory/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `shortest_path` | `(request: ShortestPathRequest) -> Result<ShortestPathResult, String>` | ✅ WIRED - Compute::Graph |
| `minimum_spanning_tree` | `(request: MSTRequest) -> Result<MSTResult, String>` | ✅ WIRED - Compute::Graph |
| `connected_components` | `(adjacency: &[Vec<bool>]) -> Vec<Vec<usize>>` | ✅ WIRED - Compute::Graph |
| `graph_properties` | `(request: GraphPropertiesRequest) -> Result<GraphPropertiesResult, String>` | ✅ WIRED - Compute::Graph |
| `topological_sort` | `(request: TopologicalSortRequest) -> Result<TopologicalSortResult, String>` | ✅ WIRED - Compute::Graph |

---

### src/specialized/information_theory/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `shannon_entropy` | `(request: EntropyRequest) -> Result<EntropyResult, String>` | ✅ WIRED - Compute::Information |

---

## src/chemistry/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_chemistry` | `(input: ChemistryInput) -> Result<ChemistryResult, String>` | ✅ WIRED - Compute::Chemistry |

---

## src/biology/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_biology` | `(input: BiologyInput) -> Result<BiologyResult, String>` | ✅ WIRED - Compute::Biology |

---

## src/thermodynamics/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_thermodynamics` | `(input: ThermodynamicsInput) -> Result<ThermodynamicsResult, String>` | ✅ WIRED - Compute::Thermodynamics |

---

## src/optics/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_optics` | `(input: OpticsInput) -> Result<OpticsResult, String>` | ✅ WIRED - Compute::Optics |

---

## src/geophysics/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_geophysics` | `(input: GeophysicsInput) -> Result<GeophysicsResult, String>` | ✅ WIRED - Compute::Geophysics |

---

## src/engineering/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_engineering` | `(input: EngineeringInput) -> Result<EngineeringResult, String>` | ✅ WIRED - Compute::Engineering |

---

## src/datetime/

#### mod.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `calculate_datetime` | `(input: DateTimeInput) -> Result<DateTimeResult, String>` | ✅ WIRED - Compute::DateTime |

---

## Summary Statistics

### Overall Counts

| Category | Count |
|----------|-------|
| **Total Functions Inventoried** | 250+ |
| **WIRED to Tool API** | ~180 functions |
| **NOT WIRED (Utilities/Helpers)** | ~70 functions |

### Wired Functions by Tool

| Tool | Function Count |
|------|----------------|
| Compute | ~95 |
| Solve | ~30 |
| Simulate | ~20 |
| Transform | ~15 |
| Differentiate | ~10 |
| Integrate | ~5 |
| Sample | ~5 |
| Analyze | ~5 |

### Not Wired Functions by Category

| Category | Count | Notes |
|----------|-------|-------|
| Helper Functions | ~40 | Internal utilities called by wired functions |
| Standalone Functions | ~30 | Potentially useful but not exposed via API |

---

## Functions Recommended for Wiring

### High Priority (Standalone with Clear Use Cases)

1. **Chaos/Fractals Module**
   - `mandelbrot`, `julia`, `burning_ship` - Fractal generation
   - `box_counting_dimension` - Fractal dimension analysis
   - `koch_snowflake` - Fractal geometry

2. **Chaos/Attractors Module**
   - `lorenz_attractor`, `rossler_attractor` - Chaotic system simulation
   - `logistic_map_bifurcation` - Bifurcation analysis

3. **Chaos/Lyapunov Module**
   - `lyapunov_1d_map`, `lyapunov_logistic_map` - Chaos quantification
   - `kaplan_yorke_dimension` - Strange attractor dimension

4. **Wormholes Module**
   - `morris_thorne_metric` - Spacetime metric computation
   - `compute_tidal_forces`, `analyze_traversal` - Traversability analysis
   - `calculate_total_energy`, `throat_energy_density` - Energy requirements

5. **Dimensional Analysis Module**
   - `dimensional_analysis` - Complete dimensional analysis
   - `convert_units` - Unit conversion
   - `are_units_compatible` - Unit compatibility checking

6. **Equation Validation Module**
   - `validate_equation` - Complete equation validation
   - `check_dimensional_consistency` - Dimensional consistency
   - `check_conservation_laws` - Physics law compliance

7. **Symbolic Regression Module**
   - `discover_equations` - Automated equation discovery from data

### Medium Priority (Useful Extensions)

1. **Game Theory Module**
   - `prisoners_dilemma`, `battle_of_sexes` - Classic games
   - Nash equilibrium computation

2. **3D Geometry**
   - `convex_hull_3d` - 3D computational geometry

3. **Symbolic Tensor Calculus**
   - `calculate_christoffel_symbols_symbolic` - Symbolic tensor operations
   - `calculate_riemann_tensor_symbolic`
   - `calculate_einstein_tensor_symbolic`

---

## Architecture Notes

### Wiring Pattern

Functions are wired through the `src/implementations/` layer:

```
User Request → ToolDispatcher → Unified Implementation → Domain Module → Result
```

### Tool Categories

1. **Solve**: Equations, systems, optimization, root finding
2. **Compute**: Calculus, transforms, field theory, sampling, matrix ops
3. **Analyze**: Series, limits, stability analysis, simplification
4. **Simulate**: Time evolution, stochastic processes, fluid dynamics
5. **Transform**: FFT, wavelets, filters, Laplace/Fourier transforms
6. **Differentiate**: Symbolic and numeric differentiation
7. **Integrate**: Definite, indefinite, contour integration
8. **FieldTheory**: EM fields, Green's functions, quantum fields
9. **Sample**: Monte Carlo, MCMC, statistical sampling
10. **Optimize**: Curve fitting, gradient descent, linear programming

### Legacy Compatibility

Legacy tool names route to primary tools:
- `differentiate` → `compute` with `operation: {differentiate: ...}`
- `integrate` → `compute` with `operation: {integrate: ...}`
- `transform` → `compute` with `operation: {transform: ...}`
- `fieldtheory` → `compute` with `operation: {field: ...}`
- `sample` → `compute` with `operation: {sample: ...}`
- `optimize` → `solve` with `equation_type: {optimize: ...}`

---

## NEWLY DISCOVERED UNWIRED MODULES

The following modules were discovered during architecture analysis and contain **REAL implementations** that are completely UNWIRED to any API.

### src/specialized/machine_learning/

**Status:** ❌ COMPLETELY UNWIRED - None of these functions are accessible via any API!

#### clustering.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `kmeans` | `(data: &[Vec<f64>], k: usize, max_iter: usize) -> Result<KMeansResult, String>` | ❌ NOT WIRED |
| `euclidean_distance_squared` | `(a: &[f64], b: &[f64]) -> f64` | ❌ NOT WIRED |
| `silhouette_score` | `(data: &[Vec<f64>], labels: &[usize]) -> Result<f64, String>` | ❌ NOT WIRED |

#### dimensionality_reduction.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `pca` | `(data: &[Vec<f64>], n_components: Option<usize>) -> Result<PCAResult, String>` | ❌ NOT WIRED |
| `transform` | `(data: &[Vec<f64>], pca_result: &PCAResult) -> Result<Vec<Vec<f64>>, String>` | ❌ NOT WIRED |
| `inverse_transform` | `(data: &[Vec<f64>], pca_result: &PCAResult) -> Result<Vec<Vec<f64>>, String>` | ❌ NOT WIRED |
| `cumulative_explained_variance` | `(pca_result: &PCAResult) -> Vec<f64>` | ❌ NOT WIRED |
| `components_for_variance` | `(pca_result: &PCAResult, target_variance: f64) -> usize` | ❌ NOT WIRED |

#### neural_network.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `activate` | `(x: &[f64], activation: ActivationFunction) -> Vec<f64>` | ❌ NOT WIRED |
| `activate_scalar` | `(x: f64, activation: ActivationFunction) -> f64` | ❌ NOT WIRED |
| `activate_softmax` | `(x: &[f64]) -> Vec<f64>` | ❌ NOT WIRED |
| `activate_derivative` | `(x: &[f64], activation: ActivationFunction) -> Vec<f64>` | ❌ NOT WIRED |
| `compute_loss` | `(y_true: &[Vec<f64>], y_pred: &[Vec<f64>], loss_fn: LossFunction) -> f64` | ❌ NOT WIRED |
| `DenseLayer::new` | Constructor for neural network layers | ❌ NOT WIRED |
| `NeuralNetwork::new` | Constructor for neural networks | ❌ NOT WIRED |

#### optimization.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `train_network` | `(network: &mut NeuralNetwork, x: &[Vec<f64>], y: &[Vec<f64>], config: &TrainingConfig) -> TrainingHistory` | ❌ NOT WIRED |
| SGD, Adam, RMSprop, Momentum optimizers | Various signatures | ❌ NOT WIRED |

#### regression.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `linear_regression` | `(x: &[Vec<f64>], y: &[f64]) -> Result<LinearRegression, String>` | ❌ NOT WIRED |
| `linear_regression_predict` | `(model: &LinearRegression, x: &[Vec<f64>]) -> Vec<f64>` | ❌ NOT WIRED |
| `ridge_regression` | `(x: &[Vec<f64>], y: &[f64], alpha: f64) -> Result<LinearRegression, String>` | ❌ NOT WIRED |
| `logistic_regression` | `(x: &[Vec<f64>], y: &[f64], config: &TrainingConfig) -> Result<LogisticRegression, String>` | ❌ NOT WIRED |
| `logistic_regression_predict` | `(model: &LogisticRegression, x: &[Vec<f64>]) -> Vec<f64>` | ❌ NOT WIRED |
| `logistic_regression_predict_proba` | `(model: &LogisticRegression, x: &[Vec<f64>]) -> Vec<f64>` | ❌ NOT WIRED |
| `mean_squared_error` | `(y_true: &[f64], y_pred: &[f64]) -> f64` | ❌ NOT WIRED |
| `mean_absolute_error` | `(y_true: &[f64], y_pred: &[f64]) -> f64` | ❌ NOT WIRED |

---

### src/mathematics/chaos/

**Status:** ❌ COMPLETELY UNWIRED - None of these functions are accessible via any API!

#### fractals.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `mandelbrot` | `(c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult` | ❌ NOT WIRED |
| `julia` | `(z0: Complex, c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult` | ❌ NOT WIRED |
| `burning_ship` | `(c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult` | ❌ NOT WIRED |
| `box_counting_dimension` | `(points: &[(f64, f64)], min_box: f64, max_box: f64, num_scales: usize) -> f64` | ❌ NOT WIRED |
| `koch_snowflake` | `(order: usize) -> Vec<(f64, f64)>` | ❌ NOT WIRED |

#### attractors.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `lorenz_attractor` | `(config: &LorenzConfig, steps: usize, dt: f64) -> AttractorTrajectory` | ❌ NOT WIRED |
| `rossler_attractor` | `(config: &RosslerConfig, steps: usize, dt: f64) -> AttractorTrajectory` | ❌ NOT WIRED |
| `logistic_map_bifurcation` | `(r_min: f64, r_max: f64, num_r: usize, iterations: usize, transient: usize) -> Vec<(f64, Vec<f64>)>` | ❌ NOT WIRED |

#### lyapunov.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `lyapunov_1d_map` | `<F>(map: F, x0: f64, iterations: usize, transient: usize) -> f64` | ❌ NOT WIRED |
| `lyapunov_logistic_map` | `(r: f64, iterations: usize, transient: usize) -> f64` | ❌ NOT WIRED |
| `lyapunov_spectrum_3d` | `<F>(...) -> [f64; 3]` | ❌ NOT WIRED |
| `kaplan_yorke_dimension` | `(exponents: &[f64]) -> f64` | ❌ NOT WIRED |

---

### src/physics/wormholes/

**Status:** ❌ COMPLETELY UNWIRED - None of these functions are accessible via any API!

#### metric.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `redshift_function` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED |
| `shape_function` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED |
| `shape_function_derivative` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED |
| `morris_thorne_metric` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> WormholeMetric` | ❌ NOT WIRED |
| `satisfies_flaring_condition` | `(r: f64, config: &WormholeConfig) -> bool` | ❌ NOT WIRED |
| `proper_distance` | `(r_start: f64, r_end: f64, config: &WormholeConfig, n_steps: usize) -> f64` | ❌ NOT WIRED |
| `embedding_surface` | `(r: f64, config: &WormholeConfig) -> f64` | ❌ NOT WIRED |

#### traversal.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_tidal_forces` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> TidalForce` | ❌ NOT WIRED |
| `analyze_traversal` | `(trajectory: &[SphericalCoordinates], config: &WormholeConfig) -> TraversalAnalysis` | ❌ NOT WIRED |
| `time_dilation_factor` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> f64` | ❌ NOT WIRED |
| `required_velocity` | `(config: &WormholeConfig, target_time: f64) -> f64` | ❌ NOT WIRED |

#### energy.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_energy_density` | `(coords: &SphericalCoordinates, config: &WormholeConfig) -> EnergyDensity` | ❌ NOT WIRED |
| `calculate_total_energy` | `(config: &WormholeConfig) -> WormholeEnergy` | ❌ NOT WIRED |
| `compare_with_black_hole` | `(config: &WormholeConfig) -> f64` | ❌ NOT WIRED |
| `throat_energy_density` | `(config: &WormholeConfig) -> f64` | ❌ NOT WIRED |

---

### src/specialized/game_theory/

**Status:** ❌ COMPLETELY UNWIRED - None of these functions are accessible via any API!

#### normal_form.rs (~270 lines)
| Function | Signature | Status |
|----------|-----------|--------|
| `prisoners_dilemma` | `() -> NormalFormGame` | ❌ NOT WIRED |
| `battle_of_sexes` | `() -> NormalFormGame` | ❌ NOT WIRED |
| `find_nash_equilibrium` | `(game: &NormalFormGame) -> Vec<Strategy>` | ❌ NOT WIRED |
| `compute_payoffs` | `(game: &NormalFormGame, strategies: &[Strategy]) -> Vec<f64>` | ❌ NOT WIRED |

#### evolutionary.rs (~83 lines)
| Function | Signature | Status |
|----------|-----------|--------|
| Evolutionary game theory functions | Various | ❌ NOT WIRED |

#### extensive_form.rs (~44 lines)
| Function | Signature | Status |
|----------|-----------|--------|
| Game tree functions | Various | ❌ NOT WIRED |

---

### src/specialized/control_theory/

**Status:** ⚠️ PARTIALLY WIRED (some functions accessible, some not)

#### pid_controller.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `ziegler_nichols_tuning` | `(ku: f64, pu: f64) -> PIDConfig` | ✅ WIRED |
| `cohen_coon_tuning` | `(k: f64, tau: f64, theta: f64) -> PIDConfig` | ✅ WIRED |

#### analysis.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `controllability` | `(a: &[Vec<f64>], b: &[Vec<f64>]) -> ControllabilityResult` | ✅ WIRED |
| `observability` | `(a: &[Vec<f64>], c: &[Vec<f64>]) -> ObservabilityResult` | ✅ WIRED |
| `eigenvalues` | `(a: &[Vec<f64>]) -> Vec<Complex64>` | ✅ WIRED |

#### lqr.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `lqr` | `(a, b, q, r) -> Result<LQRSolution, String>` | ✅ WIRED |

#### state_space.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `siso_state_space` | `(a, b, c, d, u_sequence) -> Result<Vec<f64>, String>` | ✅ WIRED |

---

### src/specialized/linear_programming/

**Status:** ⚠️ PARTIALLY WIRED

#### simplex.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `simplex` | `(lp: &LinearProgram) -> Result<LinearProgramSolution, String>` | ✅ WIRED |
| `is_feasible` | `(lp: &LinearProgram, solution: &[f64]) -> bool` | ❌ NOT WIRED |
| `objective_value` | `(lp: &LinearProgram, solution: &[f64]) -> f64` | ❌ NOT WIRED |

#### dual.rs
| Function | Signature | Status |
|----------|-----------|--------|
| `compute_dual` | `(primal: &LinearProgram) -> LinearProgram` | ❌ NOT WIRED |
| `sensitivity_analysis` | `(solution, lp) -> SensitivityAnalysis` | ❌ NOT WIRED |
| `verify_strong_duality` | `(primal_obj, dual_obj, tolerance) -> bool` | ❌ NOT WIRED |
| `interpret_shadow_price` | `(shadow_price: f64) -> String` | ❌ NOT WIRED |

---

## Updated Summary Statistics

### Overall Counts (REVISED)

| Category | Count |
|----------|-------|
| **Total Functions Discovered** | 350+ |
| **WIRED to Tool API** | ~180 functions (51%) |
| **NOT WIRED (Critical Gap)** | ~170 functions (49%) |

### Unwired by Module (Priority Order)

| Module | Functions | Lines | Priority |
|--------|-----------|-------|----------|
| `specialized/machine_learning/` | ~30 | 1,900 | **CRITICAL** |
| `mathematics/chaos/` | ~15 | 500 | **HIGH** |
| `physics/wormholes/` | ~15 | 400 | **HIGH** |
| `specialized/game_theory/` | ~10 | 440 | **MEDIUM** |
| `tools/dimensional_analysis/` | ~10 | 350 | **HIGH** |
| `tools/equation_validation/` | ~8 | 200 | **MEDIUM** |
| `specialized/linear_programming/dual` | ~4 | 170 | **LOW** |
| Various helper functions | ~78 | - | **LOW** |

### API Layer Duplication Analysis

| Layer | Location | Status | Recommendation |
|-------|----------|--------|----------------|
| Unified 10-Tool API | `src/engine/` + `src/implementations/` | INCOMPLETE | **EXPAND** |
| Module JSON API | `src/api/` | MIXED (stubs + real) | **CONSOLIDATE** |
| Legacy API | `src/api_legacy/` | DEPRECATED | **REMOVE** |

---

## Refactoring Checklist

### Phase 1: Wire Critical Missing Modules
- [ ] Wire `machine_learning/` to `Compute` or new `ML` tool
- [ ] Wire `chaos/` (fractals, attractors, Lyapunov) to `Compute` or `Simulate`
- [ ] Wire `wormholes/` to `Compute`
- [ ] Wire `dimensional_analysis/` to `Analyze`

### Phase 2: Remove Duplicated Implementations
- [ ] Remove inline math from `src/api/handlers/advanced_calculus.rs`
- [ ] Remove inline math from other API handlers
- [ ] All math must call domain modules, not be duplicated

### Phase 3: Eliminate Legacy API
- [ ] Port any unique wiring from `src/api_legacy/` to unified API
- [ ] Delete `src/api_legacy/` entirely
- [ ] Fix or delete stub handlers in `src/api/`

### Phase 4: Update Documentation
- [ ] Update CLAUDE.md with new architecture
- [ ] Update API.md with complete operation list
- [ ] Update OPERATIONS_COMPLETE_LIST.md
