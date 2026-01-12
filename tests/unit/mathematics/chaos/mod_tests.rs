// Unit tests for mathematics::chaos::mod
use computational_engine::chaos::mod::*;

use super::*;

    #[test]
    fn test_complex_operations() {
        let c1 = Complex::new(1.0, 2.0);
        let c2 = Complex::new(3.0, 4.0);

        let sum = c1 + c2;
        assert_eq!(sum.re, 4.0);
        assert_eq!(sum.im, 6.0);

        let prod = c1 * c2;
        assert_eq!(prod.re, -5.0); // (1*3 - 2*4)
        assert_eq!(prod.im, 10.0);  // (1*4 + 2*3)
    }

    #[test]
    fn test_complex_magnitude() {
        let c = Complex::new(3.0, 4.0);
        assert_eq!(c.mag(), 5.0);
        assert_eq!(c.mag_sq(), 25.0);
    }

    #[test]
    fn test_point3d_distance() {
        let p1 = Point3D::new(0.0, 0.0, 0.0);
        let p2 = Point3D::new(3.0, 4.0, 0.0);
        assert_eq!(p1.distance(&p2), 5.0);
    }
