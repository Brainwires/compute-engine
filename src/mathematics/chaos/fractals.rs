//! Fractal Generation

use super::Complex;
use serde::{Deserialize, Serialize};

/// Mandelbrot set iteration result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MandelbrotResult {
    pub iterations: usize,
    pub escaped: bool,
    pub final_magnitude: f64,
}

/// Calculate iterations for Mandelbrot set: z_(n+1) = z_n² + c
pub fn mandelbrot(c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult {
    let mut z = Complex::zero();
    let escape_sq = escape_radius * escape_radius;

    for i in 0..max_iter {
        if z.mag_sq() > escape_sq {
            return MandelbrotResult {
                iterations: i,
                escaped: true,
                final_magnitude: z.mag(),
            };
        }
        z = z * z + c;
    }

    MandelbrotResult {
        iterations: max_iter,
        escaped: false,
        final_magnitude: z.mag(),
    }
}

/// Julia set: z_(n+1) = z_n² + c, but c is fixed and z_0 varies
pub fn julia(z0: Complex, c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult {
    let mut z = z0;
    let escape_sq = escape_radius * escape_radius;

    for i in 0..max_iter {
        if z.mag_sq() > escape_sq {
            return MandelbrotResult {
                iterations: i,
                escaped: true,
                final_magnitude: z.mag(),
            };
        }
        z = z * z + c;
    }

    MandelbrotResult {
        iterations: max_iter,
        escaped: false,
        final_magnitude: z.mag(),
    }
}

/// Burning Ship fractal: z_(n+1) = (|Re(z_n)| + i|Im(z_n)|)² + c
pub fn burning_ship(c: Complex, max_iter: usize, escape_radius: f64) -> MandelbrotResult {
    let mut z = Complex::zero();
    let escape_sq = escape_radius * escape_radius;

    for i in 0..max_iter {
        if z.mag_sq() > escape_sq {
            return MandelbrotResult {
                iterations: i,
                escaped: true,
                final_magnitude: z.mag(),
            };
        }
        // Apply absolute values before squaring
        let z_abs = Complex::new(z.re.abs(), z.im.abs());
        z = z_abs * z_abs + c;
    }

    MandelbrotResult {
        iterations: max_iter,
        escaped: false,
        final_magnitude: z.mag(),
    }
}

/// Box-counting fractal dimension
pub fn box_counting_dimension(points: &[(f64, f64)], min_box: f64, max_box: f64, num_scales: usize) -> f64 {
    if points.is_empty() {
        return 0.0;
    }

    let mut log_sizes = vec![];
    let mut log_counts = vec![];

    let log_min = min_box.ln();
    let log_max = max_box.ln();
    let step = (log_max - log_min) / (num_scales as f64 - 1.0);

    for i in 0..num_scales {
        let box_size = (log_min + i as f64 * step).exp();
        let count = count_boxes(points, box_size);

        if count > 0 {
            log_sizes.push(box_size.ln());
            log_counts.push((count as f64).ln());
        }
    }

    // Linear regression: log(N) ~ -D * log(ε)
    if log_sizes.len() < 2 {
        return 0.0;
    }

    let n = log_sizes.len() as f64;
    let sum_x: f64 = log_sizes.iter().sum();
    let sum_y: f64 = log_counts.iter().sum();
    let sum_xx: f64 = log_sizes.iter().map(|x| x * x).sum();
    let sum_xy: f64 = log_sizes.iter().zip(&log_counts).map(|(x, y)| x * y).sum();

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    -slope // Dimension is negative of slope
}

fn count_boxes(points: &[(f64, f64)], box_size: f64) -> usize {
    use std::collections::HashSet;

    let mut occupied_boxes = HashSet::new();

    for &(x, y) in points {
        let box_x = (x / box_size).floor() as i64;
        let box_y = (y / box_size).floor() as i64;
        occupied_boxes.insert((box_x, box_y));
    }

    occupied_boxes.len()
}

/// Koch snowflake iteration
pub fn koch_snowflake(order: usize) -> Vec<(f64, f64)> {
    let mut points = vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (0.5, (3.0_f64).sqrt() / 2.0),
        (0.0, 0.0),
    ];

    for _ in 0..order {
        points = koch_iterate(&points);
    }

    points
}

fn koch_iterate(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut result = vec![];

    for i in 0..(points.len() - 1) {
        let p1 = points[i];
        let p2 = points[i + 1];

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;

        let p3 = (p1.0 + dx / 3.0, p1.1 + dy / 3.0);
        let p5 = (p1.0 + 2.0 * dx / 3.0, p1.1 + 2.0 * dy / 3.0);

        // Peak of triangle (60° rotation)
        let p4 = (
            p3.0 + (dx / 3.0) * 0.5 - (dy / 3.0) * (3.0_f64).sqrt() / 2.0,
            p3.1 + (dy / 3.0) * 0.5 + (dx / 3.0) * (3.0_f64).sqrt() / 2.0,
        );

        result.push(p1);
        result.push(p3);
        result.push(p4);
        result.push(p5);
    }

    result.push(*points.last().unwrap());
    result
}

