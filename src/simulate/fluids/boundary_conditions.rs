use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BoundaryType {
    NoSlip,
    FreeSlip,
    VelocityInlet { u: f64, v: f64 },
    PressureOutlet { pressure: f64 },
    Wall { u: f64, v: f64 },
    Symmetry,
    Periodic,
}

#[derive(Debug, Clone)]
pub struct BoundaryConditions {
    pub left: BoundaryType,
    pub right: BoundaryType,
    pub bottom: BoundaryType,
    pub top: BoundaryType,
}

impl BoundaryConditions {
    pub fn new(
        left: BoundaryType,
        right: BoundaryType,
        bottom: BoundaryType,
        top: BoundaryType,
    ) -> Self {
        Self {
            left,
            right,
            bottom,
            top,
        }
    }

    pub fn all_no_slip() -> Self {
        Self {
            left: BoundaryType::NoSlip,
            right: BoundaryType::NoSlip,
            bottom: BoundaryType::NoSlip,
            top: BoundaryType::NoSlip,
        }
    }

    pub fn lid_driven_cavity(lid_velocity: f64) -> Self {
        Self {
            left: BoundaryType::NoSlip,
            right: BoundaryType::NoSlip,
            bottom: BoundaryType::NoSlip,
            top: BoundaryType::Wall {
                u: lid_velocity,
                v: 0.0,
            },
        }
    }

    pub fn channel_flow(inlet_velocity: f64) -> Self {
        Self {
            left: BoundaryType::VelocityInlet {
                u: inlet_velocity,
                v: 0.0,
            },
            right: BoundaryType::PressureOutlet { pressure: 0.0 },
            bottom: BoundaryType::NoSlip,
            top: BoundaryType::NoSlip,
        }
    }

    pub fn apply_velocity_boundary(
        &self,
        u: &mut ndarray::Array2<f64>,
        v: &mut ndarray::Array2<f64>,
    ) {
        let (nx, ny) = u.dim();

        // Left boundary (i = 0)
        for j in 0..ny {
            match &self.left {
                BoundaryType::NoSlip => {
                    u[[0, j]] = 0.0;
                    v[[0, j]] = 0.0;
                }
                BoundaryType::FreeSlip => {
                    u[[0, j]] = 0.0; // Normal velocity = 0
                    // v[[0, j]] remains unchanged (free slip in tangential direction)
                }
                BoundaryType::VelocityInlet {
                    u: u_inlet,
                    v: v_inlet,
                } => {
                    u[[0, j]] = *u_inlet;
                    v[[0, j]] = *v_inlet;
                }
                BoundaryType::Wall {
                    u: u_wall,
                    v: v_wall,
                } => {
                    u[[0, j]] = *u_wall;
                    v[[0, j]] = *v_wall;
                }
                BoundaryType::Symmetry => {
                    u[[0, j]] = 0.0; // Normal velocity = 0
                    // Tangential velocity gradient = 0 (handled in solver)
                }
                _ => {} // Other types handled elsewhere
            }
        }

        // Right boundary (i = nx-1)
        for j in 0..ny {
            match &self.right {
                BoundaryType::NoSlip => {
                    u[[nx - 1, j]] = 0.0;
                    v[[nx - 1, j]] = 0.0;
                }
                BoundaryType::FreeSlip => {
                    u[[nx - 1, j]] = 0.0;
                }
                BoundaryType::Wall {
                    u: u_wall,
                    v: v_wall,
                } => {
                    u[[nx - 1, j]] = *u_wall;
                    v[[nx - 1, j]] = *v_wall;
                }
                BoundaryType::PressureOutlet { .. } => {
                    // Zero gradient condition (handled in solver)
                    u[[nx - 1, j]] = u[[nx - 2, j]];
                    v[[nx - 1, j]] = v[[nx - 2, j]];
                }
                _ => {}
            }
        }

        // Bottom boundary (j = 0)
        for i in 0..nx {
            match &self.bottom {
                BoundaryType::NoSlip => {
                    u[[i, 0]] = 0.0;
                    v[[i, 0]] = 0.0;
                }
                BoundaryType::FreeSlip => {
                    v[[i, 0]] = 0.0;
                }
                BoundaryType::Wall {
                    u: u_wall,
                    v: v_wall,
                } => {
                    u[[i, 0]] = *u_wall;
                    v[[i, 0]] = *v_wall;
                }
                BoundaryType::Symmetry => {
                    v[[i, 0]] = 0.0;
                }
                _ => {}
            }
        }

        // Top boundary (j = ny-1)
        for i in 0..nx {
            match &self.top {
                BoundaryType::NoSlip => {
                    u[[i, ny - 1]] = 0.0;
                    v[[i, ny - 1]] = 0.0;
                }
                BoundaryType::FreeSlip => {
                    v[[i, ny - 1]] = 0.0;
                }
                BoundaryType::Wall {
                    u: u_wall,
                    v: v_wall,
                } => {
                    u[[i, ny - 1]] = *u_wall;
                    v[[i, ny - 1]] = *v_wall;
                }
                BoundaryType::Symmetry => {
                    v[[i, ny - 1]] = 0.0;
                }
                _ => {}
            }
        }
    }

    pub fn apply_pressure_boundary(&self, p: &mut ndarray::Array2<f64>) {
        let (nx, ny) = p.dim();

        // Left boundary
        for j in 0..ny {
            match &self.left {
                BoundaryType::PressureOutlet { pressure } => {
                    p[[0, j]] = *pressure;
                }
                BoundaryType::VelocityInlet { .. }
                | BoundaryType::NoSlip
                | BoundaryType::FreeSlip
                | BoundaryType::Wall { .. } => {
                    // Neumann condition (zero gradient)
                    p[[0, j]] = p[[1, j]];
                }
                _ => {}
            }
        }

        // Right boundary
        for j in 0..ny {
            match &self.right {
                BoundaryType::PressureOutlet { pressure } => {
                    p[[nx - 1, j]] = *pressure;
                }
                _ => {
                    // Neumann condition (zero gradient)
                    p[[nx - 1, j]] = p[[nx - 2, j]];
                }
            }
        }

        // Bottom boundary
        for i in 0..nx {
            match &self.bottom {
                BoundaryType::PressureOutlet { pressure } => {
                    p[[i, 0]] = *pressure;
                }
                _ => {
                    // Neumann condition (zero gradient)
                    p[[i, 0]] = p[[i, 1]];
                }
            }
        }

        // Top boundary
        for i in 0..nx {
            match &self.top {
                BoundaryType::PressureOutlet { pressure } => {
                    p[[i, ny - 1]] = *pressure;
                }
                _ => {
                    // Neumann condition (zero gradient)
                    p[[i, ny - 1]] = p[[i, ny - 2]];
                }
            }
        }
    }
}

impl Default for BoundaryConditions {
    fn default() -> Self {
        Self::all_no_slip()
    }
}
