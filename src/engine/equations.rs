//! Equation and Model type definitions
//!
//! This module contains all the specific equation types, models, and operations
//! that map the 180+ original operations to the 10 core tools.

use serde::{Deserialize, Serialize};

// ============================================================================
// SOLVE TOOL - Equation Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EquationType {
    // Einstein Field Equations
    Einstein(EinsteinEquation),
    // Fluid Dynamics
    Fluid(FluidEquation),
    // Differential Equations
    Differential(DifferentialEquation),
    // Electromagnetic
    Electromagnetic(EMEquation),
    // Chemical
    Chemical(ChemicalEquation),
    // Linear Systems
    LinearSystem,
    // Root Finding
    RootFinding,
    // Number Theory
    NumberTheory(NumberTheoryProblem),
    // Differential Geometry
    DifferentialGeometry(DiffGeoProblem),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EinsteinEquation {
    Vacuum,
    WithSource,
    Schwarzschild,
    KerrNewman,
    FriedmannRobertsonWalker,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FluidEquation {
    NavierStokes,
    CavityFlow,
    ChannelFlow,
    LidDrivenCavity,
    Euler,
    Bernoulli,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DifferentialEquation {
    ODE,
    PDE,
    BoundaryValue,
    InitialValue,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EMEquation {
    Maxwell,
    Wave,
    TransmissionLine,
    Waveguide,
    Helmholtz,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChemicalEquation {
    Balance,
    Thermodynamic,
    Kinetics,
    GasLaw,
    AcidBase,
    Electrochemistry,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NumberTheoryProblem {
    DiscreteLog,
    Factorization,
    PrimalityTest,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiffGeoProblem {
    Geodesic,
    ParallelTransport,
    MinimalSurface,
}

// ============================================================================
// DIFFERENTIATE TOOL - Operation Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DifferentiationOp {
    // Vector Calculus
    VectorCalc(VectorCalcOp),
    // Tensor Calculus
    TensorCalc(TensorDiffOp),
    // Variational Calculus
    Variational,
    // Differential Forms
    DifferentialForms,
    // Numeric
    Numeric,
    // Symbolic
    Symbolic,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VectorCalcOp {
    Gradient,
    Divergence,
    Curl,
    Laplacian,
    Directional,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TensorDiffOp {
    Covariant,
    Lie,
    ExteriorDerivative,
}

// ============================================================================
// INTEGRATE TOOL - Integration Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IntegrationType {
    // Line/Surface/Volume
    Geometric(GeometricIntegral),
    // Theorems
    Theorem(IntegralTheorem),
    // Complex Analysis
    ComplexAnalysis(ComplexIntegral),
    // Numeric
    Numeric(NumericIntegration),
    // Symbolic
    Symbolic,
    // Monte Carlo
    MonteCarlo,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeometricIntegral {
    Line,
    Surface,
    Volume,
    Contour,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IntegralTheorem {
    Greens,
    Stokes,
    Divergence,
    CauchyIntegral,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ComplexIntegral {
    Residue,
    Cauchy,
    Contour,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NumericIntegration {
    Trapezoidal,
    Simpson,
    GaussQuadrature,
    Adaptive,
}

// ============================================================================
// ANALYZE TOOL - Analysis Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AnalysisOp {
    // Expression Operations
    Simplify,
    Parse,
    ExtractVariables,
    // Validation
    Validate,
    CheckCorrectness,
    CheckDimensions,
    CheckPhysics,
    CheckConservation,
    CheckSymmetries,
    // Series and Limits
    PartialFraction,
    SeriesExpansion,
    LaurentSeries,
    Limit,
    // Field Analysis
    FieldAnalysis(FieldAnalysisType),
    // Dimensional Analysis
    DimensionalCheck,
    ValidateDimensions,
    InferDimensions,
    ScaleAnalysis,
    // Units
    UnitsDerive,
    UnitsAnalyze,
    // Graph Theory
    GraphComponents,
    GraphProperties,
    // Number Theory
    IsPrime,
    // Fluid Analysis
    FluidAnalysis,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FieldAnalysisType {
    Vector,
    Scalar,
    Tensor,
}

// ============================================================================
// SIMULATE TOOL - Simulation Models
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SimulationModel {
    // Stochastic Processes
    Stochastic(StochasticProcess),
    // Finance Models
    Finance(FinanceModel),
    // Fluid Dynamics
    FluidDynamics(FluidSim),
    // ODE/PDE Time-stepping
    TimeEvolution(TimeEvolutionMethod),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StochasticProcess {
    BrownianMotion,
    GeometricBrownian,
    OrnsteinUhlenbeck,
    Poisson,
    Levy,
    JumpDiffusion,
    FractionalBrownian,
    MeanReverting,
    VarianceGamma,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FinanceModel {
    Heston,
    SABR,
    StochasticVolatility,
    BlackScholes,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FluidSim {
    NavierStokes,
    Euler,
    LatticeBotzmann,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TimeEvolutionMethod {
    Euler,
    RungeKutta4,
    AdaptiveStep,
    ImplicitEuler,
}

// ============================================================================
// COMPUTE TOOL - Computation Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ComputeOp {
    // Tensor Operations
    Tensor(TensorOp),
    // Matrix Operations
    Matrix(MatrixOp),
    // Matrix Decompositions
    MatrixDecomp(MatrixDecomp),
    // Special Functions
    SpecialFunc(SpecialFunction),
    // Number Theory
    NumberTheory(NumberTheoryOp),
    // Computational Geometry
    Geometry(GeometryOp),
    // Information Theory
    Information(InformationOp),
    // Fourier Series
    FourierSeries,
    // EM Calculations
    EM(EMComputation),
    // Chemistry
    Chemistry(ChemistryOp),
    // Graph Theory
    Graph(GraphOp),
    // Scientific Formulas (2025 Expansion)
    Biology(BiologyOp),
    Thermodynamics(ThermodynamicsOp),
    Optics(OpticsOp),
    Geophysics(GeophysicsOp),
    Engineering(EngineeringOp),
    DateTime(DateTimeOp),
    // Physics (Tier 1 Wolfram Alpha expansion)
    Physics(PhysicsOp),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TensorOp {
    Christoffel,
    Riemann,
    Ricci,
    RicciScalar,
    Einstein,
    Weyl,
    Product,
    Contraction,
    ParallelTransport,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixOp {
    Norm,
    Power,
    Exp,
    Rank,
    Pseudoinverse,
    Determinant,
    Trace,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixDecomp {
    QR,
    SVD,
    Eigen,
    PCA,
    Cholesky,
    LU,
    Schur,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpecialFunction {
    Bessel,
    Gamma,
    Erf,
    Elliptic,
    OrthogonalPoly,
    Airy,
    Hypergeometric,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NumberTheoryOp {
    GeneratePrime,
    ModExp,
    ModInv,
    GCD,
    LCM,
    EulerTotient,
    CarmichaelLambda,
    ECPointAdd,
    // Cryptographic Operations
    RSAKeypair,
    RSAEncrypt,
    RSADecrypt,
    SHA256,
    SHA3_256,
    ChineseRemainder,
    DiscreteLog,
    PrimalityTest,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeometryOp {
    ConvexHull,
    Delaunay,
    Voronoi,
    PolygonArea,
    PointInPolygon,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InformationOp {
    Entropy,
    MutualInfo,
    ChannelCapacity,
    Huffman,
    Kolmogorov,
    ConditionalEntropy,
    KLDivergence,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EMComputation {
    PoyntingVector,
    SkinEffect,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChemistryOp {
    // Legacy operations
    MolarMass,
    EquilibriumConstant,
    ReactionRate,
    // New chemistry operations (2025 expansion)
    PhCalculation,
    BufferCapacity,
    Arrhenius,
    RateLaw,
    GibbsFreeEnergy,
    NernstEquation,
    BeerLambert,
    VanDerWaals,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BiologyOp {
    MichaelisMenten,
    Pharmacokinetics,
    HardyWeinberg,
    GoldmanEquation,
    AllometricScaling,
    GrowthModel,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ThermodynamicsOp {
    Conduction,
    Convection,
    Radiation,
    ThermalResistance,
    Entropy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OpticsOp {
    ThinLens,
    SnellsLaw,
    DiffractionGrating,
    FresnelEquations,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeophysicsOp {
    Seismology,
    Atmosphere,
    RadiometricDating,
    PlanetaryScience,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineeringOp {
    SoundPressureLevel,
    DopplerEffect,
    ReverberationTime,
    Stress,
    Strain,
    FractureMechanics,
    Bernoulli,
    Poiseuille,
    Drag,
    PidController,
    FirstOrderResponse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DateTimeOp {
    AddInterval,
    SubtractInterval,
    DateDifference,
    AgeCurrent,
    AgeAtDate,
    BusinessDays,
    AddBusinessDays,
    IsLeapYear,
    DaysInMonth,
    WeekNumber,
    DayOfWeek,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GraphOp {
    TopologicalSort,
    ShortestPath,
    MinimumSpanningTree,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PhysicsOp {
    // Relativity (12 operations)
    Relativity(RelativityOp),
    // Statistical Physics (10 operations)
    StatisticalPhysics(StatPhysicsOp),
    // Quantum Mechanics (15 operations)
    QuantumMechanics(QuantumMechOp),
    // Control Systems (12 operations)
    ControlSystems(ControlSystemsOp),
    // Nuclear Physics (8 operations)
    NuclearPhysics(NuclearOp),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RelativityOp {
    // Special Relativity
    LorentzTransform,
    TimeDilation,
    LengthContraction,
    RelativisticEnergy,
    VelocityAddition,
    // General Relativity
    SchwarzschildMetric,
    GravitationalTimeDilation,
    OrbitalPrecession,
    GravitationalLensing,
    // Black Holes
    BlackHoleProperties,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StatPhysicsOp {
    PartitionFunction,
    MaxwellBoltzmann,
    FermiDirac,
    BoseEinstein,
    PartitionFunctionCanonical,
    PartitionFunctionGrandCanonical,
    PhaseTransition,
    CriticalPhenomena,
    ChemicalPotential,
    FugacityCoefficient,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QuantumMechOp {
    SchrodingerEquation,
    HarmonicOscillator,
    HydrogenAtom,
    AngularMomentum,
    SpinOperators,
    PerturbationTheory,
    TunnelingProbability,
    DensityMatrix,
    EntanglementMeasure,
    QuantumEntropy,
    QuantumCoherence,
    BellInequality,
    QuantumTomography,
    QuantumGate,
    QuantumCircuit,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ControlSystemsOp {
    TransferFunction,
    PoleZeroAnalysis,
    BodePlot,
    NyquistPlot,
    RootLocus,
    StateSpace,
    Controllability,
    Observability,
    RouthHurwitz,
    GainMargin,
    PhaseMargin,
    StepResponse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NuclearOp {
    RadioactiveDecay,
    DecayChain,
    HalfLife,
    BindingEnergy,
    MassDefect,
    FissionEnergy,
    FusionEnergy,
    NuclearReaction,
}

// ============================================================================
// TRANSFORM TOOL - Transform Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TransformType {
    Fourier(FourierTransform),
    Laplace(LaplaceTransform),
    Wavelet(WaveletType),
    FFT(FFTType),
    Filter(FilterType),
    Window(WindowType),
    Conformal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FourierTransform {
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LaplaceTransform {
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WaveletType {
    Haar,
    Daubechies,
    Morlet,
    Mexican,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FFTType {
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FilterType {
    LowPass,
    HighPass,
    BandPass,
    BandStop,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum WindowType {
    Hamming,
    Hanning,
    Blackman,
    Kaiser,
}

// ============================================================================
// FIELDTHEORY TOOL - Field Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FieldType {
    EM(EMField),
    GreenFunction,
    QuantumField(QuantumFieldType),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EMField {
    Antenna,
    Waveguide,
    Scattering,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QuantumFieldType {
    ScalarField,
    DiracField,
    GaugeField,
}

// ============================================================================
// SAMPLE TOOL - Sampling Methods
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SamplingMethod {
    // Stochastic
    PathGeneration,
    Moments,
    // Monte Carlo
    MonteCarlo(MonteCarloMethod),
    // Statistical
    Stats(StatisticalMethod),
    // Signal Processing
    SignalAnalysis(SignalMethod),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MonteCarloMethod {
    Integration,
    MCMC,
    MetropolisHastings,
    Gibbs,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StatisticalMethod {
    HypothesisTest,
    ANOVA,
    Regression,
    TimeSeries,
    BasicStats,
    Correlation,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SignalMethod {
    SpectralAnalysis,
    Autocorrelation,
    CrossCorrelation,
    PowerSpectrum,
    Coherence,
    Cepstrum,
    PeakDetection,
}

// ============================================================================
// OPTIMIZE TOOL - Optimization Methods
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OptimizationMethod {
    // Curve Fitting
    Fit(FitMethod),
    // Minimization
    Minimize(MinimizationMethod),
    // Interpolation
    Interpolation(InterpolationMethod),
    // Dimensional Analysis
    DimensionalAnalysis(DimAnalysisMethod),
    // Symbolic Regression (Evolutionary Function Discovery)
    SymbolicRegression,
    // Automatic Model Selection
    Auto {
        #[serde(default)]
        criteria: SelectionCriteria,
        #[serde(default)]
        candidates: Vec<String>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SelectionCriteria {
    Aic,   // Akaike Information Criterion
    Bic,   // Bayesian Information Criterion
    Aicc,  // Corrected AIC for small samples
    RSquared, // Simple R² comparison
}

impl Default for SelectionCriteria {
    fn default() -> Self {
        SelectionCriteria::Aicc
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FitMethod {
    Polynomial,
    Exponential,
    Logarithmic,
    PowerLaw,
    Rational,
    Trigonometric,
    Custom,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MinimizationMethod {
    GradientDescent,
    NelderMead,
    ConjugateGradient,
    BFGS,
    LevenbergMarquardt,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InterpolationMethod {
    Linear,
    Polynomial,
    Spline,
    Cubic,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DimAnalysisMethod {
    BuckinghamPi,
    DimensionlessGroups,
    SimilarityAnalysis,
}
