//! Unified Analyzer Module
//!
//! Routes analysis requests to equation validation, dimensional analysis, and simplification modules.
//!
//! This module provides the ANALYZE tool implementation for:
//! - Symbolic operations (parse, simplify, extract variables, partial fractions)
//! - Series expansions (Taylor, Laurent, limits)
//! - Validation (equations, dimensions, physics, conservation, symmetries)
//! - Field analysis (vector, scalar, tensor fields)
//! - Graph analysis (components, properties)
//! - Fluid analysis (flow regime)

pub mod series;
pub mod stability;
pub mod symbolic;
pub mod validation;

use crate::engine::*;

pub struct UnifiedAnalyzer;

impl UnifiedAnalyzer {
    pub fn new() -> Self {
        Self
    }
}

impl Default for UnifiedAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

impl Analyze for UnifiedAnalyzer {
    fn analyze(&self, input: &AnalyzeInput) -> ToolResult<AnalyzeOutput> {
        match &input.operation {
            // Symbolic operations
            AnalysisOp::Parse => symbolic::analyze_parse(input),
            AnalysisOp::Simplify => symbolic::analyze_simplify(input),
            AnalysisOp::ExtractVariables => symbolic::analyze_extract_variables(input),
            AnalysisOp::PartialFraction => symbolic::analyze_partial_fraction(input),

            // Series operations
            AnalysisOp::SeriesExpansion => series::analyze_series_expansion(input),
            AnalysisOp::LaurentSeries => series::analyze_laurent_series(input),
            AnalysisOp::Limit => series::analyze_limit(input),

            // Validation operations
            AnalysisOp::Validate => validation::analyze_validate(input),
            AnalysisOp::CheckCorrectness => validation::analyze_validate(input),
            AnalysisOp::CheckDimensions => validation::analyze_dimensions(input),
            AnalysisOp::CheckPhysics => validation::analyze_physics(input),
            AnalysisOp::CheckConservation => validation::analyze_conservation(input),
            AnalysisOp::CheckSymmetries => validation::analyze_symmetries(input),
            AnalysisOp::DimensionalCheck => validation::analyze_dimensions(input),
            AnalysisOp::ValidateDimensions => validation::analyze_dimensions(input),
            AnalysisOp::InferDimensions => validation::analyze_infer_dimensions(input),
            AnalysisOp::ScaleAnalysis => validation::analyze_scale(input),
            AnalysisOp::UnitsDerive => validation::analyze_units_derive(input),
            AnalysisOp::UnitsAnalyze => validation::analyze_units(input),

            // Stability/field analysis operations
            AnalysisOp::FieldAnalysis(field_type) => stability::analyze_field(field_type, input),
            AnalysisOp::GraphComponents => stability::analyze_graph_components(input),
            AnalysisOp::GraphProperties => stability::analyze_graph_properties(input),
            AnalysisOp::FluidAnalysis => stability::analyze_fluid(input),

            // Number theory
            AnalysisOp::IsPrime => {
                let n: u64 = input
                    .expression
                    .parse()
                    .map_err(|_| "Expression must be a positive integer")?;

                let is_prime = stability::is_prime_simple(n);

                Ok(AnalyzeOutput {
                    result: serde_json::json!(is_prime),
                    latex: None,
                    validation: None,
                    details: Some(serde_json::json!({
                        "number": n,
                        "is_prime": is_prime
                    })),
                })
            }
        }
    }
}
