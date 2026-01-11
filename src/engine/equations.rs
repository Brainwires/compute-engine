//! Equation and Model type definitions
//!
//! This module contains all the specific equation types, models, and operations
//! that map the 180+ original operations to the 10 core tools.

#[cfg(feature = "schemars")]
use schemars::JsonSchema;
use serde::{Deserialize, Deserializer, Serialize};
use strum::EnumIter;

// ============================================================================
// SOLVE TOOL - Equation Types
// ============================================================================

/// Custom deserializer for EquationType that accepts BOTH:
/// - Simple strings: "einstein" → Einstein(Vacuum) (uses default sub-type)
/// - Case-insensitive: "Einstein", "EINSTEIN" all work
/// - Nested objects: {"einstein": "schwarzschild"} → Einstein(Schwarzschild)
///
/// This is the UNIVERSAL solution for all 916 operations across all 10 tools
#[derive(Debug, Clone, Serialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum EquationType {
    Einstein(EinsteinEquation),
    Fluid(FluidEquation),
    Differential(DifferentialEquation),
    Electromagnetic(EMEquation),
    Chemical(ChemicalEquation),
    LinearSystem,
    RootFinding,
    NumberTheory(NumberTheoryProblem),
    DifferentialGeometry(DiffGeoProblem),
    /// Optimization problems (curve fitting, minimization, etc.)
    Optimize(OptimizationMethod),
}

// Custom deserializer implementation
impl<'de> Deserialize<'de> for EquationType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        use serde::de::{self, MapAccess, Visitor};
        use std::fmt;

        struct EquationTypeVisitor;

        impl<'de> Visitor<'de> for EquationTypeVisitor {
            type Value = EquationType;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("equation type as string or object")
            }

            fn visit_str<E>(self, value: &str) -> Result<EquationType, E>
            where
                E: de::Error,
            {
                // Accept simple strings with default sub-types
                match value.to_lowercase().as_str() {
                    "einstein" => Ok(EquationType::Einstein(EinsteinEquation::Vacuum)),
                    "fluid" => Ok(EquationType::Fluid(FluidEquation::NavierStokes)),
                    "differential" => Ok(EquationType::Differential(DifferentialEquation::ODE)),
                    "electromagnetic" => Ok(EquationType::Electromagnetic(EMEquation::Maxwell)),
                    "chemical" => Ok(EquationType::Chemical(ChemicalEquation::Balance)),
                    "linear_system" | "linearsystem" => Ok(EquationType::LinearSystem),
                    "root_finding" | "rootfinding" => Ok(EquationType::RootFinding),
                    "number_theory" | "numbertheory" => Ok(EquationType::NumberTheory(NumberTheoryProblem::PrimalityTest)),
                    "differential_geometry" | "differentialgeometry" => Ok(EquationType::DifferentialGeometry(DiffGeoProblem::Geodesic)),
                    "optimize" | "optimization" => Ok(EquationType::Optimize(OptimizationMethod::default())),
                    _ => Err(E::custom(format!("unknown equation type: {}", value))),
                }
            }

            fn visit_map<M>(self, mut map: M) -> Result<EquationType, M::Error>
            where
                M: MapAccess<'de>,
            {
                // For nested objects like {"einstein": "schwarzschild"}
                if let Some(key) = map.next_key::<String>()? {
                    let key_lower = key.to_lowercase();
                    match key_lower.as_str() {
                        "einstein" => {
                            let sub: EinsteinEquation = map.next_value()?;
                            Ok(EquationType::Einstein(sub))
                        }
                        "fluid" => {
                            let sub: FluidEquation = map.next_value()?;
                            Ok(EquationType::Fluid(sub))
                        }
                        "differential" => {
                            let sub: DifferentialEquation = map.next_value()?;
                            Ok(EquationType::Differential(sub))
                        }
                        "electromagnetic" => {
                            let sub: EMEquation = map.next_value()?;
                            Ok(EquationType::Electromagnetic(sub))
                        }
                        "chemical" => {
                            let sub: ChemicalEquation = map.next_value()?;
                            Ok(EquationType::Chemical(sub))
                        }
                        "number_theory" | "numbertheory" => {
                            let sub: NumberTheoryProblem = map.next_value()?;
                            Ok(EquationType::NumberTheory(sub))
                        }
                        "differential_geometry" | "differentialgeometry" => {
                            let sub: DiffGeoProblem = map.next_value()?;
                            Ok(EquationType::DifferentialGeometry(sub))
                        }
                        "optimize" | "optimization" => {
                            let sub: OptimizationMethod = map.next_value()?;
                            Ok(EquationType::Optimize(sub))
                        }
                        _ => Err(de::Error::custom(format!("unknown equation type: {}", key))),
                    }
                } else {
                    Err(de::Error::custom("expected equation type"))
                }
            }
        }

        deserializer.deserialize_any(EquationTypeVisitor)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum EinsteinEquation {
    #[default]
    Vacuum,
    WithSource,
    Schwarzschild,
    KerrNewman,
    FriedmannRobertsonWalker,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FluidEquation {
    #[default]
    NavierStokes,
    CavityFlow,
    ChannelFlow,
    LidDrivenCavity,
    Euler,
    Bernoulli,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DifferentialEquation {
    #[default]
    ODE,
    PDE,
    BoundaryValue,
    InitialValue,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum EMEquation {
    #[default]
    Maxwell,
    Wave,
    TransmissionLine,
    Waveguide,
    Helmholtz,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ChemicalEquation {
    #[default]
    Balance,
    Thermodynamic,
    Kinetics,
    GasLaw,
    AcidBase,
    Electrochemistry,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum NumberTheoryProblem {
    #[default]
    DiscreteLog,
    Factorization,
    PrimalityTest,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DiffGeoProblem {
    #[default]
    Geodesic,
    ParallelTransport,
    MinimalSurface,
}

// ============================================================================
// DIFFERENTIATE TOOL - Operation Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DifferentiationOp {
    // Vector Calculus
    #[serde(alias = "VectorCalc")]
    VectorCalc(VectorCalcOp),
    // Tensor Calculus
    #[serde(alias = "TensorCalc")]
    TensorCalc(TensorDiffOp),
    // Variational Calculus
    #[default]
    #[serde(alias = "Variational")]
    Variational,
    // Differential Forms
    #[serde(alias = "DifferentialForms")]
    DifferentialForms,
    // Numeric
    #[serde(alias = "Numeric")]
    Numeric,
    // Symbolic
    #[serde(alias = "Symbolic")]
    Symbolic,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum VectorCalcOp {
    #[default]
    Gradient,
    Divergence,
    Curl,
    Laplacian,
    Directional,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum TensorDiffOp {
    #[default]
    Covariant,
    Lie,
    ExteriorDerivative,
}

// ============================================================================
// INTEGRATE TOOL - Integration Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum IntegrationType {
    // Line/Surface/Volume
    #[serde(alias = "Geometric")]
    Geometric(GeometricIntegral),
    // Theorems
    #[serde(alias = "Theorem")]
    Theorem(IntegralTheorem),
    // Complex Analysis
    #[serde(alias = "ComplexAnalysis")]
    ComplexAnalysis(ComplexIntegral),
    // Numeric
    #[serde(alias = "Numeric")]
    Numeric(NumericIntegration),
    // Symbolic
    #[default]
    #[serde(alias = "Symbolic")]
    Symbolic,
    // Monte Carlo
    #[serde(alias = "MonteCarlo")]
    MonteCarlo,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum GeometricIntegral {
    #[default]
    Line,
    Surface,
    Volume,
    Contour,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum IntegralTheorem {
    #[default]
    Greens,
    Stokes,
    Divergence,
    CauchyIntegral,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ComplexIntegral {
    #[default]
    Residue,
    Cauchy,
    Contour,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum NumericIntegration {
    #[default]
    Trapezoidal,
    Simpson,
    GaussQuadrature,
    Adaptive,
}

// ============================================================================
// ANALYZE TOOL - Analysis Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum AnalysisOp {
    #[default]
    // Expression Operations
    #[serde(alias = "Simplify")]
    Simplify,
    #[serde(alias = "Parse")]
    Parse,
    #[serde(alias = "ExtractVariables")]
    ExtractVariables,
    // Validation
    #[serde(alias = "Validate")]
    Validate,
    #[serde(alias = "CheckCorrectness")]
    CheckCorrectness,
    #[serde(alias = "CheckDimensions")]
    CheckDimensions,
    #[serde(alias = "CheckPhysics")]
    CheckPhysics,
    #[serde(alias = "CheckConservation")]
    CheckConservation,
    #[serde(alias = "CheckSymmetries")]
    CheckSymmetries,
    // Series and Limits
    #[serde(alias = "PartialFraction")]
    PartialFraction,
    #[serde(alias = "SeriesExpansion")]
    SeriesExpansion,
    #[serde(alias = "LaurentSeries")]
    LaurentSeries,
    #[serde(alias = "Limit")]
    Limit,
    // Field Analysis
    #[serde(alias = "FieldAnalysis")]
    FieldAnalysis(FieldAnalysisType),
    // Dimensional Analysis
    #[serde(alias = "DimensionalCheck")]
    DimensionalCheck,
    #[serde(alias = "ValidateDimensions")]
    ValidateDimensions,
    #[serde(alias = "InferDimensions")]
    InferDimensions,
    #[serde(alias = "ScaleAnalysis")]
    ScaleAnalysis,
    // Units
    #[serde(alias = "UnitsDerive")]
    UnitsDerive,
    #[serde(alias = "UnitsAnalyze")]
    UnitsAnalyze,
    // Graph Theory
    #[serde(alias = "GraphComponents")]
    GraphComponents,
    #[serde(alias = "GraphProperties")]
    GraphProperties,
    // Number Theory
    #[serde(alias = "IsPrime")]
    IsPrime,
    // Fluid Analysis
    #[serde(alias = "FluidAnalysis")]
    FluidAnalysis,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FieldAnalysisType {
    #[default]
    Vector,
    Scalar,
    Tensor,
}

// ============================================================================
// SIMULATE TOOL - Simulation Models
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum SimulationModel {
    // Stochastic Processes
    #[serde(alias = "Stochastic")]
    Stochastic(StochasticProcess),
    // Finance Models
    #[serde(alias = "Finance")]
    Finance(FinanceModel),
    // Fluid Dynamics
    #[serde(alias = "FluidDynamics")]
    FluidDynamics(FluidSim),
    // ODE/PDE Time-stepping
    #[serde(alias = "TimeEvolution")]
    TimeEvolution(TimeEvolutionMethod),
}

impl Default for SimulationModel {
    fn default() -> Self {
        SimulationModel::TimeEvolution(TimeEvolutionMethod::default())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum StochasticProcess {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FinanceModel {
    #[default]
    Heston,
    SABR,
    StochasticVolatility,
    BlackScholes,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FluidSim {
    #[default]
    NavierStokes,
    Euler,
    LatticeBotzmann,
    /// 1D Quantum Navier-Stokes with Bohm potential correction
    QuantumNavierStokes1D,
    /// 2D Quantum Navier-Stokes with Bohm potential correction
    QuantumNavierStokes2D,
    /// 3D incompressible Navier-Stokes (Taylor-Green benchmark)
    NavierStokes3D,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum TimeEvolutionMethod {
    #[default]
    Euler,
    RungeKutta4,
    AdaptiveStep,
    ImplicitEuler,
}

// ============================================================================
// COMPUTE TOOL - Computation Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ComputeOp {
    // ========== TENSOR & MATRIX ==========
    /// Tensor operations (Christoffel symbols, Riemann, Ricci, etc.)
    #[serde(alias = "Tensor")]
    Tensor(TensorOp),
    /// Matrix operations (norm, power, etc.)
    #[serde(alias = "Matrix")]
    Matrix(MatrixOp),
    /// Matrix decompositions (LU, QR, SVD, etc.)
    #[serde(alias = "MatrixDecomp")]
    MatrixDecomp(MatrixDecomp),

    // ========== CALCULUS (absorbed from Differentiate/Integrate tools) ==========
    /// Differentiation operations
    #[serde(alias = "Differentiate")]
    Differentiate(DifferentiationOp),
    /// Integration operations
    #[serde(alias = "Integrate")]
    Integrate(IntegrationType),

    // ========== TRANSFORMS (absorbed from Transform tool) ==========
    /// Signal transforms (FFT, Laplace, wavelet, etc.)
    #[serde(alias = "Transform")]
    Transform(TransformType),

    // ========== FIELD THEORY (absorbed from FieldTheory tool) ==========
    /// Field theory operations (EM, quantum, Bohm potential, decoherence)
    #[serde(alias = "Field")]
    Field(FieldType),

    // ========== SAMPLING & STATISTICS (absorbed from Sample tool) ==========
    /// Statistical sampling and analysis
    #[serde(alias = "Sample")]
    Sample(SamplingMethod),

    // ========== SPECIAL FUNCTIONS ==========
    /// Special mathematical functions (Bessel, Gamma, etc.)
    #[serde(alias = "SpecialFunc")]
    SpecialFunc(SpecialFunction),

    // ========== NUMBER THEORY ==========
    /// Number theory operations
    #[serde(alias = "NumberTheory")]
    NumberTheory(NumberTheoryOp),

    // ========== GEOMETRY ==========
    /// Computational geometry
    #[serde(alias = "Geometry")]
    Geometry(GeometryOp),

    // ========== INFORMATION THEORY ==========
    /// Information theory (entropy, mutual information)
    #[serde(alias = "Information")]
    Information(InformationOp),

    // ========== FOURIER SERIES ==========
    /// Fourier series computation
    #[serde(alias = "FourierSeries")]
    FourierSeries,

    // ========== PHYSICS & CHEMISTRY ==========
    /// EM calculations
    #[serde(alias = "EM")]
    EM(EMComputation),
    /// Chemistry operations
    #[serde(alias = "Chemistry")]
    Chemistry(ChemistryOp),
    /// Physics formulas
    #[serde(alias = "Physics")]
    Physics(PhysicsOp),

    // ========== GRAPH THEORY ==========
    /// Graph theory operations
    #[serde(alias = "Graph")]
    Graph(GraphOp),

    // ========== SCIENTIFIC FORMULAS (2025 Expansion) ==========
    /// Biology formulas
    #[serde(alias = "Biology")]
    Biology(BiologyOp),
    /// Thermodynamics
    #[serde(alias = "Thermodynamics")]
    Thermodynamics(ThermodynamicsOp),
    /// Optics
    #[serde(alias = "Optics")]
    Optics(OpticsOp),
    /// Geophysics
    #[serde(alias = "Geophysics")]
    Geophysics(GeophysicsOp),
    /// Engineering formulas
    #[serde(alias = "Engineering")]
    Engineering(EngineeringOp),
    /// Date/time calculations
    #[serde(alias = "DateTime")]
    DateTime(DateTimeOp),
}

impl Default for ComputeOp {
    fn default() -> Self {
        ComputeOp::Tensor(TensorOp::default())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum TensorOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MatrixOp {
    #[default]
    Norm,
    Power,
    Exp,
    Rank,
    Pseudoinverse,
    Determinant,
    Trace,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MatrixDecomp {
    #[default]
    QR,
    SVD,
    Eigen,
    PCA,
    Cholesky,
    LU,
    Schur,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum SpecialFunction {
    #[default]
    Bessel,
    Gamma,
    Erf,
    Elliptic,
    OrthogonalPoly,
    Airy,
    Hypergeometric,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum NumberTheoryOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum GeometryOp {
    #[default]
    ConvexHull,
    Delaunay,
    Voronoi,
    PolygonArea,
    PointInPolygon,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum InformationOp {
    #[default]
    Entropy,
    MutualInfo,
    ChannelCapacity,
    Huffman,
    Kolmogorov,
    ConditionalEntropy,
    KLDivergence,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum EMComputation {
    #[default]
    PoyntingVector,
    SkinEffect,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ChemistryOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum BiologyOp {
    #[default]
    MichaelisMenten,
    Pharmacokinetics,
    HardyWeinberg,
    GoldmanEquation,
    AllometricScaling,
    GrowthModel,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ThermodynamicsOp {
    #[default]
    Conduction,
    Convection,
    Radiation,
    ThermalResistance,
    Entropy,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum OpticsOp {
    #[default]
    ThinLens,
    SnellsLaw,
    DiffractionGrating,
    FresnelEquations,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum GeophysicsOp {
    #[default]
    Seismology,
    Atmosphere,
    RadiometricDating,
    PlanetaryScience,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum EngineeringOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DateTimeOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum GraphOp {
    #[default]
    TopologicalSort,
    ShortestPath,
    MinimumSpanningTree,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
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

impl Default for PhysicsOp {
    fn default() -> Self {
        PhysicsOp::Relativity(RelativityOp::default())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum RelativityOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum StatPhysicsOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum QuantumMechOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ControlSystemsOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum NuclearOp {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum TransformType {
    #[serde(alias = "Fourier")]
    Fourier(FourierTransform),
    #[serde(alias = "Laplace")]
    Laplace(LaplaceTransform),
    #[serde(alias = "Wavelet")]
    Wavelet(WaveletType),
    #[serde(alias = "FFT")]
    FFT(FFTType),
    #[serde(alias = "Filter")]
    Filter(FilterType),
    #[serde(alias = "Window")]
    Window(WindowType),
    #[default]
    #[serde(alias = "Conformal")]
    Conformal,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FourierTransform {
    #[default]
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum LaplaceTransform {
    #[default]
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum WaveletType {
    #[default]
    Haar,
    Daubechies,
    Morlet,
    Mexican,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FFTType {
    #[default]
    Forward,
    Inverse,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FilterType {
    #[default]
    LowPass,
    HighPass,
    BandPass,
    BandStop,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum WindowType {
    #[default]
    Hamming,
    Hanning,
    Blackman,
    Kaiser,
}

// ============================================================================
// FIELDTHEORY TOOL - Field Types
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FieldType {
    #[serde(alias = "EM")]
    EM(EMField),
    #[default]
    #[serde(alias = "GreenFunction")]
    GreenFunction,
    #[serde(alias = "QuantumField")]
    QuantumField(QuantumFieldType),
    /// Bohm quantum potential Q = (ℏ²/2m) × ∇²√ρ / √ρ
    #[serde(alias = "bohm_potential")]
    BohmPotential,
    /// Decoherence length scale L_D = ℏ / √(2 m k_B T)
    #[serde(alias = "decoherence_scale")]
    DecoherenceScale,
    /// Quantum stress tensor for fluid dynamics
    #[serde(alias = "quantum_stress")]
    QuantumStress,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum EMField {
    #[default]
    Antenna,
    Waveguide,
    Scattering,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum QuantumFieldType {
    #[default]
    ScalarField,
    DiracField,
    GaugeField,
}

// ============================================================================
// SAMPLE TOOL - Sampling Methods
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum SamplingMethod {
    // Stochastic
    #[default]
    #[serde(alias = "PathGeneration")]
    PathGeneration,
    #[serde(alias = "Moments")]
    Moments,
    // Monte Carlo
    #[serde(alias = "MonteCarlo")]
    MonteCarlo(MonteCarloMethod),
    // Statistical
    #[serde(alias = "Stats")]
    Stats(StatisticalMethod),
    // Signal Processing
    #[serde(alias = "SignalAnalysis")]
    SignalAnalysis(SignalMethod),
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MonteCarloMethod {
    #[default]
    Integration,
    MCMC,
    MetropolisHastings,
    Gibbs,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum StatisticalMethod {
    #[default]
    HypothesisTest,
    ANOVA,
    Regression,
    TimeSeries,
    BasicStats,
    Correlation,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum SignalMethod {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum OptimizationMethod {
    // Curve Fitting
    #[serde(alias = "Fit")]
    Fit(FitMethod),
    // Minimization
    #[serde(alias = "Minimize")]
    Minimize(MinimizationMethod),
    // Interpolation
    #[serde(alias = "Interpolation")]
    Interpolation(InterpolationMethod),
    // Dimensional Analysis
    #[serde(alias = "DimensionalAnalysis")]
    DimensionalAnalysis(DimAnalysisMethod),
    // Symbolic Regression (Evolutionary Function Discovery)
    #[serde(alias = "SymbolicRegression")]
    SymbolicRegression,
    // Automatic Model Selection
    #[serde(alias = "Auto")]
    Auto {
        #[serde(default)]
        criteria: SelectionCriteria,
        #[serde(default)]
        candidates: Vec<String>,
    },
}

impl Default for OptimizationMethod {
    fn default() -> Self {
        OptimizationMethod::Fit(FitMethod::default())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum SelectionCriteria {
    Aic,      // Akaike Information Criterion
    Bic,      // Bayesian Information Criterion
    #[default]
    Aicc,     // Corrected AIC for small samples
    RSquared, // Simple R² comparison
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FitMethod {
    #[default]
    Polynomial,
    Exponential,
    Logarithmic,
    PowerLaw,
    Rational,
    Trigonometric,
    Custom,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MinimizationMethod {
    #[default]
    GradientDescent,
    NelderMead,
    ConjugateGradient,
    BFGS,
    LevenbergMarquardt,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum InterpolationMethod {
    #[default]
    Linear,
    Polynomial,
    Spline,
    Cubic,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DimAnalysisMethod {
    #[default]
    BuckinghamPi,
    DimensionlessGroups,
    SimilarityAnalysis,
}

// ============================================================================
// ML TOOL - Machine Learning Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MLOp {
    /// Clustering algorithms (K-means, DBSCAN, etc.)
    #[serde(alias = "Clustering")]
    Clustering(ClusteringMethod),
    /// Neural network operations
    #[serde(alias = "NeuralNetwork")]
    NeuralNetwork(NeuralNetworkOp),
    /// Regression methods
    #[serde(alias = "Regression")]
    Regression(RegressionMethod),
    /// Dimensionality reduction
    #[serde(alias = "DimReduction")]
    DimReduction(DimReductionMethod),
    /// Classification algorithms
    #[serde(alias = "Classification")]
    Classification(ClassificationMethod),
}

impl Default for MLOp {
    fn default() -> Self {
        MLOp::Clustering(ClusteringMethod::default())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ClusteringMethod {
    #[default]
    KMeans,
    DBSCAN,
    Hierarchical,
    GaussianMixture,
    /// Compute silhouette score for clustering quality
    SilhouetteScore,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum NeuralNetworkOp {
    #[default]
    /// Create a new layer
    CreateLayer,
    /// Forward propagation
    Forward,
    /// Backward propagation
    Backward,
    /// Train the network
    Train,
    /// Predict using trained network
    Predict,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum RegressionMethod {
    #[default]
    Linear,
    Logistic,
    Ridge,
    Lasso,
    ElasticNet,
    PolynomialRegression,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DimReductionMethod {
    #[default]
    PCA,
    TSNE,
    UMAP,
    LDA,
    /// Transform data using fitted model
    Transform,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ClassificationMethod {
    #[default]
    DecisionTree,
    RandomForest,
    SVM,
    NaiveBayes,
    KNN,
}

// ============================================================================
// CHAOS TOOL - Chaos Theory Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ChaosOp {
    /// Fractal generation (Mandelbrot, Julia, etc.)
    #[serde(alias = "Fractal")]
    Fractal(FractalType),
    /// Strange attractors (Lorenz, Rossler, etc.)
    #[serde(alias = "Attractor")]
    Attractor(AttractorType),
    /// Lyapunov exponent calculation
    #[serde(alias = "Lyapunov")]
    Lyapunov(LyapunovMethod),
    /// Bifurcation analysis
    #[serde(alias = "Bifurcation")]
    Bifurcation(BifurcationType),
    /// Fractal dimension calculation
    #[serde(alias = "Dimension")]
    Dimension(DimensionMethod),
}

impl Default for ChaosOp {
    fn default() -> Self {
        ChaosOp::Fractal(FractalType::default())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum FractalType {
    #[default]
    Mandelbrot,
    Julia,
    BurningShip,
    KochSnowflake,
    SierpinskiTriangle,
    DragonCurve,
    Cantor,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum AttractorType {
    #[default]
    Lorenz,
    Rossler,
    Henon,
    Chua,
    Thomas,
    Aizawa,
    Chen,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum LyapunovMethod {
    #[default]
    Map1D,
    Spectrum1D,
    Spectrum3D,
    MaxLyapunov,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum BifurcationType {
    #[default]
    LogisticMap,
    PeriodDoubling,
    Pitchfork,
    SaddleNode,
    Hopf,
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum DimensionMethod {
    #[default]
    BoxCounting,
    CorrelationDimension,
    HausdorffDimension,
    KaplanYorke,
}

// ============================================================================
// UNITS TOOL - Dimensional Analysis Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum UnitsOp {
    /// Convert between units
    #[default]
    #[serde(alias = "Convert")]
    Convert,
    /// Analyze dimensions of an expression
    #[serde(alias = "Analyze")]
    Analyze,
    /// Check if two units are compatible
    #[serde(alias = "CheckCompatibility")]
    CheckCompatibility,
    /// Get SI base units for a quantity
    #[serde(alias = "GetBase")]
    GetBase,
    /// Derive units for a quantity
    #[serde(alias = "Derive")]
    Derive,
    /// Parse unit string
    #[serde(alias = "Parse")]
    Parse,
    /// Simplify unit expression
    #[serde(alias = "Simplify")]
    Simplify,
}

// ============================================================================
// VALIDATE TOOL - Equation/Physics Validation Operations
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize, EnumIter, Default)]
#[cfg_attr(feature = "schemars", derive(JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum ValidateOp {
    /// Validate equation syntax and structure
    #[default]
    #[serde(alias = "Equation")]
    Equation,
    /// Check dimensional consistency
    #[serde(alias = "Dimensions")]
    Dimensions,
    /// Check conservation laws (energy, momentum, charge)
    #[serde(alias = "Conservation")]
    Conservation,
    /// Check symmetry properties
    #[serde(alias = "Symmetry")]
    Symmetry,
    /// Check physics compliance (causality, positivity, etc.)
    #[serde(alias = "Physics")]
    Physics,
    /// Validate mathematical bounds and constraints
    #[serde(alias = "Bounds")]
    Bounds,
    /// Check for singularities
    #[serde(alias = "Singularities")]
    Singularities,
}
