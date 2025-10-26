use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Grid2D {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

impl Grid2D {
    pub fn new(nx: usize, ny: usize, x_min: f64, x_max: f64, y_min: f64, y_max: f64) -> Self {
        let dx = (x_max - x_min) / (nx - 1) as f64;
        let dy = (y_max - y_min) / (ny - 1) as f64;

        Self {
            nx,
            ny,
            dx,
            dy,
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    pub fn uniform(nx: usize, ny: usize, length_x: f64, length_y: f64) -> Self {
        Self::new(nx, ny, 0.0, length_x, 0.0, length_y)
    }

    pub fn x_coord(&self, i: usize) -> f64 {
        self.x_min + (i as f64) * self.dx
    }

    pub fn y_coord(&self, j: usize) -> f64 {
        self.y_min + (j as f64) * self.dy
    }

    pub fn x_coordinates(&self) -> Vec<f64> {
        (0..self.nx).map(|i| self.x_coord(i)).collect()
    }

    pub fn y_coordinates(&self) -> Vec<f64> {
        (0..self.ny).map(|j| self.y_coord(j)).collect()
    }

    pub fn is_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || i == self.nx - 1 || j == 0 || j == self.ny - 1
    }

    pub fn is_interior(&self, i: usize, j: usize) -> bool {
        !self.is_boundary(i, j)
    }

    pub fn cfl_condition(&self, u_max: f64, v_max: f64) -> f64 {
        let dt_convective = 0.5 / (u_max / self.dx + v_max / self.dy);
        dt_convective
    }

    pub fn diffusion_condition(&self, viscosity: f64) -> f64 {
        let dt_diffusive =
            0.25 / (viscosity * (1.0 / (self.dx * self.dx) + 1.0 / (self.dy * self.dy)));
        dt_diffusive
    }

    pub fn stable_timestep(&self, u_max: f64, v_max: f64, viscosity: f64) -> f64 {
        let dt_cfl = self.cfl_condition(u_max, v_max);
        let dt_diffusion = self.diffusion_condition(viscosity);
        dt_cfl.min(dt_diffusion) * 0.8 // Safety factor
    }
}

#[derive(Debug, Clone)]
pub struct StaggeredGrid2D {
    pub grid: Grid2D,
}

impl StaggeredGrid2D {
    pub fn new(grid: Grid2D) -> Self {
        Self { grid }
    }

    /// x-coordinate for u-velocity (staggered in x-direction)
    pub fn x_u(&self, i: usize) -> f64 {
        self.grid.x_min + (i as f64 + 0.5) * self.grid.dx
    }

    /// y-coordinate for v-velocity (staggered in y-direction)  
    pub fn y_v(&self, j: usize) -> f64 {
        self.grid.y_min + (j as f64 + 0.5) * self.grid.dy
    }

    /// Check if u-velocity point is on boundary
    pub fn is_u_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || i == self.grid.nx || j == 0 || j == self.grid.ny - 1
    }

    /// Check if v-velocity point is on boundary
    pub fn is_v_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || i == self.grid.nx - 1 || j == 0 || j == self.grid.ny
    }
}

