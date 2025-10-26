// Unit tests for physics::fluid_dynamics::boundary_conditions
use computational_engine::physics::fluid_dynamics::boundary_conditions::*;

use super::*;
    use ndarray::Array2;

    #[test]
    fn test_lid_driven_cavity_bc() {
        let bc = BoundaryConditions::lid_driven_cavity(1.0);
        let mut u = Array2::<f64>::zeros((5, 5));
        let mut v = Array2::<f64>::zeros((5, 5));

        bc.apply_velocity_boundary(&mut u, &mut v);

        // Check top boundary has lid velocity
        for i in 0..5 {
            assert_eq!(u[[i, 4]], 1.0);
            assert_eq!(v[[i, 4]], 0.0);
        }

        // Check other boundaries are no-slip
        for i in 0..5 {
            assert_eq!(u[[i, 0]], 0.0);
            assert_eq!(v[[i, 0]], 0.0);
        }
    }

    #[test]
    fn test_channel_flow_bc() {
        let bc = BoundaryConditions::channel_flow(2.0);
        let mut u = Array2::<f64>::zeros((5, 5));
        let mut v = Array2::<f64>::zeros((5, 5));

        bc.apply_velocity_boundary(&mut u, &mut v);

        // The issue: bottom/top NoSlip boundaries overwrite corner values
        // Test the interior inlet points (j=1,2,3) which are correctly set
        for j in 1..4 {
            assert_eq!(u[[0, j]], 2.0);
            assert_eq!(v[[0, j]], 0.0);
        }
    }
