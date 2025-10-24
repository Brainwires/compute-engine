//! Barnes-Hut Algorithm for Fast N-Body Simulation
//!
//! Reduces complexity from O(N²) to O(N log N) using octree spatial decomposition.
//!
//! **Algorithm:**
//! 1. Build octree of space containing all bodies
//! 2. For each body, traverse tree and:
//!    - If node is far enough (θ < s/d), treat as single mass
//!    - Otherwise, recurse into children
//!
//! **Opening angle θ:** Tradeoff between speed and accuracy (typical: 0.5)

use super::{Body, Vec3};
use serde::{Deserialize, Serialize};

/// Octree node for Barnes-Hut algorithm
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OctreeNode {
    pub center: Vec3,
    pub half_width: f64,
    pub total_mass: f64,
    pub center_of_mass: Vec3,
    pub children: Option<Box<[OctreeNode; 8]>>,
    pub body_index: Option<usize>,
}

impl OctreeNode {
    /// Create empty octree node
    pub fn new(center: Vec3, half_width: f64) -> Self {
        Self {
            center,
            half_width,
            total_mass: 0.0,
            center_of_mass: Vec3::zero(),
            children: None,
            body_index: None,
        }
    }

    /// Insert body into octree
    pub fn insert(&mut self, body: &Body, index: usize) {
        // If node is empty, place body here
        if self.total_mass == 0.0 {
            self.total_mass = body.mass;
            self.center_of_mass = body.position;
            self.body_index = Some(index);
            return;
        }

        // If node contains a single body, subdivide
        if self.children.is_none() && self.body_index.is_some() {
            self.subdivide();
        }

        // Update center of mass
        let new_total_mass = self.total_mass + body.mass;
        self.center_of_mass = (self.center_of_mass * self.total_mass + body.position * body.mass) / new_total_mass;
        self.total_mass = new_total_mass;

        // Insert into appropriate child
        let octant = self.get_octant(&body.position);
        if let Some(ref mut children) = self.children {
            children[octant].insert(body, index);
        }
    }

    /// Subdivide node into 8 children
    fn subdivide(&mut self) {
        let quarter = self.half_width / 2.0;
        let offsets = [
            Vec3::new(-quarter, -quarter, -quarter),
            Vec3::new( quarter, -quarter, -quarter),
            Vec3::new(-quarter,  quarter, -quarter),
            Vec3::new( quarter,  quarter, -quarter),
            Vec3::new(-quarter, -quarter,  quarter),
            Vec3::new( quarter, -quarter,  quarter),
            Vec3::new(-quarter,  quarter,  quarter),
            Vec3::new( quarter,  quarter,  quarter),
        ];

        let children: Vec<OctreeNode> = offsets.iter()
            .map(|&offset| OctreeNode::new(self.center + offset, quarter))
            .collect();

        self.children = Some(Box::new([
            children[0].clone(), children[1].clone(), children[2].clone(), children[3].clone(),
            children[4].clone(), children[5].clone(), children[6].clone(), children[7].clone(),
        ]));
    }

    /// Determine which octant a position belongs to
    fn get_octant(&self, pos: &Vec3) -> usize {
        let mut octant = 0;
        if pos.x > self.center.x { octant |= 1; }
        if pos.y > self.center.y { octant |= 2; }
        if pos.z > self.center.z { octant |= 4; }
        octant
    }

    /// Calculate force using Barnes-Hut approximation
    pub fn calculate_force(&self, body: &Body, theta: f64, softening: f64) -> Vec3 {
        // Skip if this is the body itself
        if let Some(idx) = self.body_index {
            if body.name == format!("body_{}", idx) {
                return Vec3::zero();
            }
        }

        let r = self.center_of_mass - body.position;
        let dist = r.mag();

        // Check opening criterion: s/d < θ
        if self.children.is_none() || (2.0 * self.half_width / dist) < theta {
            // Treat as single mass
            if dist < 1e-10 {
                return Vec3::zero();
            }
            let dist_sq = dist * dist + softening * softening;
            let f_mag = super::G * body.mass * self.total_mass / dist_sq;
            return r.normalize() * f_mag;
        }

        // Recurse into children
        let mut force = Vec3::zero();
        if let Some(ref children) = self.children {
            for child in children.iter() {
                if child.total_mass > 0.0 {
                    force = force + child.calculate_force(body, theta, softening);
                }
            }
        }
        force
    }
}

/// Build octree from bodies
pub fn build_octree(bodies: &[Body]) -> OctreeNode {
    // Find bounding box
    let mut min = bodies[0].position;
    let mut max = bodies[0].position;

    for body in bodies {
        min.x = min.x.min(body.position.x);
        min.y = min.y.min(body.position.y);
        min.z = min.z.min(body.position.z);
        max.x = max.x.max(body.position.x);
        max.y = max.y.max(body.position.y);
        max.z = max.z.max(body.position.z);
    }

    let center = (min + max) * 0.5;
    let half_width = ((max.x - min.x).max(max.y - min.y).max(max.z - min.z)) / 2.0 * 1.1; // Slight margin

    let mut root = OctreeNode::new(center, half_width);

    for (i, body) in bodies.iter().enumerate() {
        root.insert(body, i);
    }

    root
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_octree_build() {
        let bodies = vec![
            Body::new(1.0, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "b1"),
            Body::new(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::zero(), "b2"),
            Body::new(1.0, Vec3::new(0.0, 1.0, 0.0), Vec3::zero(), "b3"),
        ];

        let tree = build_octree(&bodies);

        assert_eq!(tree.total_mass, 3.0);
        assert!(tree.center_of_mass.x > 0.0);
    }

    #[test]
    fn test_octree_force() {
        let bodies = vec![
            Body::new(1e30, Vec3::new(0.0, 0.0, 0.0), Vec3::zero(), "sun"),
            Body::new(1e24, Vec3::new(1e11, 0.0, 0.0), Vec3::zero(), "earth"),
        ];

        let tree = build_octree(&bodies);
        let test_body = Body::new(1.0, Vec3::new(5e10, 0.0, 0.0), Vec3::zero(), "probe");

        let force = tree.calculate_force(&test_body, 0.5, 0.0);

        // Force should point away from origin (towards sun's mass concentration)
        assert!(force.mag() > 0.0);
    }
}
