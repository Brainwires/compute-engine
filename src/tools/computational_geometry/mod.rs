//! Computational geometry module
//!
//! Provides geometric algorithms and spatial computations:
//! - Convex hull algorithms (2D and 3D)
//! - Voronoi diagrams
//! - Triangulation
//! - Spatial indexing (KD-trees)

mod advanced;
mod lib;
pub mod spatial_3d;

pub use advanced::*;
pub use spatial_3d::*;
