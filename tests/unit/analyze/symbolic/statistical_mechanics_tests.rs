// Unit tests for mathematics::symbolic_cas::statistical_mechanics
use computational_engine::analyze::symbolic::statistical_mechanics::*;

use super::*;

    #[test]
    fn test_boltzmann_distribution() {
        let p = boltzmann_distribution("E", "T", None);
        println!("Boltzmann distribution: {}", p);
    }

    #[test]
    fn test_partition_functions() {
        let z_discrete = partition_function_discrete();
        let z_gas = partition_function_ideal_gas();

        println!("Z (discrete): {}", z_discrete);
        println!("Z (ideal gas): {}", z_gas);
    }

    #[test]
    fn test_free_energies() {
        let f = helmholtz_free_energy("T", None);
        let g = gibbs_free_energy();

        println!("Helmholtz free energy: {}", f);
        println!("Gibbs free energy: {}", g);
    }

    #[test]
    fn test_maxwell_boltzmann() {
        let f = maxwell_boltzmann_speed_distribution();
        println!("Maxwell-Boltzmann: {}", f);
    }

    #[test]
    fn test_quantum_distributions() {
        let fd = fermi_dirac_distribution("E", "μ", "T");
        let be = bose_einstein_distribution("E", "μ", "T");

        println!("Fermi-Dirac: {}", fd);
        println!("Bose-Einstein: {}", be);
    }

    #[test]
    fn test_planck_distribution() {
        let u = planck_distribution();
        println!("Planck distribution: {}", u);
    }

    #[test]
    fn test_ideal_gas_law() {
        let pv = ideal_gas_law();
        println!("Ideal gas law: {}", pv);
    }

    #[test]
    fn test_van_der_waals() {
        let vdw = van_der_waals_equation();
        println!("Van der Waals: {}", vdw);
    }

    #[test]
    fn test_carnot_efficiency() {
        let eta = carnot_efficiency("T_c", "T_h");
        println!("Carnot efficiency: {}", eta);
    }
