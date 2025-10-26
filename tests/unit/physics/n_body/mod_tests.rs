// Unit tests for physics::n_body::mod
use computational_engine::physics::n_body::mod::*;

use super::*;

    #[test]
    fn test_vec3_operations() {
        let v1 = Vec3::new(1.0, 2.0, 3.0);
        let v2 = Vec3::new(4.0, 5.0, 6.0);

        let sum = v1 + v2;
        assert_eq!(sum.x, 5.0);
        assert_eq!(sum.y, 7.0);
        assert_eq!(sum.z, 9.0);

        let diff = v2 - v1;
        assert_eq!(diff.x, 3.0);

        let mag = Vec3::new(3.0, 4.0, 0.0).mag();
        assert!((mag - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_body_energy() {
        let body = Body::new(
            1.0,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(10.0, 0.0, 0.0),
            "test"
        );

        let ke = body.kinetic_energy();
        assert!((ke - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_body_system() {
        // Earth-Sun system (simplified)
        let sun = Body::new(
            1.989e30,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            "Sun"
        );

        let earth = Body::new(
            5.972e24,
            Vec3::new(1.496e11, 0.0, 0.0), // 1 AU
            Vec3::new(0.0, 29780.0, 0.0),   // Orbital velocity
            "Earth"
        );

        let config = NBodyConfig {
            bodies: vec![sun, earth],
            dt: 3600.0, // 1 hour
            steps: 100,
            method: IntegrationMethod::Verlet,
            softening: 0.0,
        };

        let force = config.total_force(1);
        assert!(force.mag() > 0.0);

        let pe = config.potential_energy();
        assert!(pe < 0.0); // Bound system
    }

    #[test]
    fn test_system_properties() {
        let bodies = vec![
            Body::new(1.0, Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), "b1"),
            Body::new(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0), "b2"),
        ];

        let config = NBodyConfig {
            bodies,
            dt: 0.01,
            steps: 10,
            method: IntegrationMethod::Verlet,
            softening: 0.01,
        };

        let props = config.system_properties();

        // Equal masses with opposite velocities → zero total momentum
        assert!(props.total_momentum.mag() < 1e-10);

        // Center of mass should be at (0.5, 0, 0)
        assert!((props.center_of_mass.x - 0.5).abs() < 1e-10);
    }
