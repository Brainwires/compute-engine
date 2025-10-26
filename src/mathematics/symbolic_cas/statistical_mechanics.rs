//! Statistical mechanics and thermodynamics
//!
//! Partition functions, entropy, free energy, and statistical distributions

use super::Expr;
use super::simplify::simplify;

/// Boltzmann distribution: P(E) = exp(-E/kT) / Z
///
/// Returns the probability expression
pub fn boltzmann_distribution(energy: &str, temperature: &str, k_b: Option<&str>) -> Expr {
    let e = Expr::sym(energy);
    let t = Expr::sym(temperature);
    let k = Expr::sym(k_b.unwrap_or("k_B"));

    // exp(-E/kT)
    let exponent = Expr::mul(
        Expr::num(-1),
        Expr::mul(e, Expr::pow(Expr::mul(k, t), Expr::num(-1))),
    );

    let numerator = Expr::func("exp", vec![exponent]);

    // Divide by partition function Z
    simplify(&Expr::mul(
        numerator,
        Expr::pow(Expr::sym("Z"), Expr::num(-1)),
    ))
}

/// Partition function for discrete states: Z = Σ exp(-E_i/kT)
///
/// Returns symbolic representation
pub fn partition_function_discrete() -> Expr {
    Expr::func(
        "Σ",
        vec![Expr::func(
            "exp",
            vec![Expr::mul(
                Expr::num(-1),
                Expr::mul(
                    Expr::sym("E_i"),
                    Expr::pow(Expr::mul(Expr::sym("k_B"), Expr::sym("T")), Expr::num(-1)),
                ),
            )],
        )],
    )
}

/// Canonical partition function for classical ideal gas
///
/// Z = (V/λ³)^N / N!  where λ is thermal wavelength
pub fn partition_function_ideal_gas() -> Expr {
    let v_over_lambda_cubed = Expr::pow(
        Expr::mul(Expr::sym("V"), Expr::pow(Expr::sym("λ"), Expr::num(-3))),
        Expr::sym("N"),
    );

    // Divide by N!
    Expr::mul(
        v_over_lambda_cubed,
        Expr::pow(Expr::func("factorial", vec![Expr::sym("N")]), Expr::num(-1)),
    )
}

/// Helmholtz free energy: F = -kT ln(Z)
pub fn helmholtz_free_energy(temperature: &str, k_b: Option<&str>) -> Expr {
    let t = Expr::sym(temperature);
    let k = Expr::sym(k_b.unwrap_or("k_B"));

    // F = -kT ln(Z)
    simplify(&Expr::mul(
        Expr::mul(Expr::num(-1), Expr::mul(k, t)),
        Expr::func("ln", vec![Expr::sym("Z")]),
    ))
}

/// Gibbs free energy: G = H - TS = U + PV - TS
pub fn gibbs_free_energy() -> Expr {
    let u = Expr::sym("U");
    let p = Expr::sym("P");
    let v = Expr::sym("V");
    let t = Expr::sym("T");
    let s = Expr::sym("S");

    // G = U + PV - TS
    simplify(&Expr::add(
        Expr::add(u, Expr::mul(p, v)),
        Expr::mul(Expr::num(-1), Expr::mul(t, s)),
    ))
}

/// Entropy from partition function: S = k(ln(Z) + T ∂ln(Z)/∂T)
///
/// Symbolic representation
pub fn entropy_from_partition_function() -> Expr {
    let k = Expr::sym("k_B");
    let ln_z = Expr::func("ln", vec![Expr::sym("Z")]);

    // T ∂ln(Z)/∂T
    let derivative_term = Expr::mul(Expr::sym("T"), Expr::func("∂/∂T", vec![ln_z.clone()]));

    // S = k(ln(Z) + T ∂ln(Z)/∂T)
    simplify(&Expr::mul(k, Expr::add(ln_z, derivative_term)))
}

/// Maxwell-Boltzmann speed distribution
///
/// f(v) = 4π(m/2πkT)^(3/2) v² exp(-mv²/2kT)
pub fn maxwell_boltzmann_speed_distribution() -> Expr {
    let m = Expr::sym("m");
    let k = Expr::sym("k_B");
    let t = Expr::sym("T");
    let v = Expr::sym("v");

    // Prefactor: 4π(m/2πkT)^(3/2)
    let prefactor = Expr::mul(
        Expr::mul(Expr::num(4), Expr::sym("π")),
        Expr::pow(
            Expr::mul(
                m.clone(),
                Expr::pow(
                    Expr::mul(
                        Expr::mul(Expr::num(2), Expr::sym("π")),
                        Expr::mul(k.clone(), t.clone()),
                    ),
                    Expr::num(-1),
                ),
            ),
            Expr::rational_unchecked(3, 2),
        ),
    );

    // v² term
    let v_squared = Expr::pow(v.clone(), Expr::num(2));

    // Exponential: exp(-mv²/2kT)
    let exponent = Expr::mul(
        Expr::num(-1),
        Expr::mul(
            Expr::mul(m, v_squared.clone()),
            Expr::pow(Expr::mul(Expr::mul(Expr::num(2), k), t), Expr::num(-1)),
        ),
    );

    let exponential = Expr::func("exp", vec![exponent]);

    // Combine
    simplify(&Expr::mul(Expr::mul(prefactor, v_squared), exponential))
}

/// Fermi-Dirac distribution (fermions)
///
/// f(E) = 1 / (exp((E-μ)/kT) + 1)
pub fn fermi_dirac_distribution(energy: &str, chemical_potential: &str, temperature: &str) -> Expr {
    let e = Expr::sym(energy);
    let mu = Expr::sym(chemical_potential);
    let k = Expr::sym("k_B");
    let t = Expr::sym(temperature);

    // (E - μ) / kT
    let exponent = Expr::mul(
        Expr::add(e, Expr::mul(Expr::num(-1), mu)),
        Expr::pow(Expr::mul(k, t), Expr::num(-1)),
    );

    // exp((E-μ)/kT) + 1
    let denominator = Expr::add(Expr::func("exp", vec![exponent]), Expr::num(1));

    simplify(&Expr::pow(denominator, Expr::num(-1)))
}

/// Bose-Einstein distribution (bosons)
///
/// f(E) = 1 / (exp((E-μ)/kT) - 1)
pub fn bose_einstein_distribution(
    energy: &str,
    chemical_potential: &str,
    temperature: &str,
) -> Expr {
    let e = Expr::sym(energy);
    let mu = Expr::sym(chemical_potential);
    let k = Expr::sym("k_B");
    let t = Expr::sym(temperature);

    // (E - μ) / kT
    let exponent = Expr::mul(
        Expr::add(e, Expr::mul(Expr::num(-1), mu)),
        Expr::pow(Expr::mul(k, t), Expr::num(-1)),
    );

    // exp((E-μ)/kT) - 1
    let denominator = Expr::add(Expr::func("exp", vec![exponent]), Expr::num(-1));

    simplify(&Expr::pow(denominator, Expr::num(-1)))
}

/// Planck distribution (blackbody radiation)
///
/// u(ν) = (8πhν³/c³) / (exp(hν/kT) - 1)
pub fn planck_distribution() -> Expr {
    let h = Expr::sym("h");
    let nu = Expr::sym("ν");
    let c = Expr::sym("c");
    let k = Expr::sym("k_B");
    let t = Expr::sym("T");

    // 8πhν³/c³
    let prefactor = Expr::mul(
        Expr::mul(
            Expr::mul(Expr::num(8), Expr::sym("π")),
            Expr::mul(h.clone(), Expr::pow(nu.clone(), Expr::num(3))),
        ),
        Expr::pow(Expr::pow(c, Expr::num(3)), Expr::num(-1)),
    );

    // exp(hν/kT) - 1
    let exponent = Expr::mul(Expr::mul(h, nu), Expr::pow(Expr::mul(k, t), Expr::num(-1)));

    let denominator = Expr::add(Expr::func("exp", vec![exponent]), Expr::num(-1));

    simplify(&Expr::mul(prefactor, Expr::pow(denominator, Expr::num(-1))))
}

/// Stefan-Boltzmann law for blackbody radiation
///
/// j = σT⁴
pub fn stefan_boltzmann_law(temperature: &str) -> Expr {
    let sigma = Expr::sym("σ");
    let t = Expr::sym(temperature);

    simplify(&Expr::mul(sigma, Expr::pow(t, Expr::num(4))))
}

/// Ideal gas law: PV = NkT
pub fn ideal_gas_law() -> Expr {
    let p = Expr::sym("P");
    let v = Expr::sym("V");
    let n = Expr::sym("N");
    let k = Expr::sym("k_B");
    let t = Expr::sym("T");

    // PV - NkT = 0
    simplify(&Expr::add(
        Expr::mul(p, v),
        Expr::mul(Expr::num(-1), Expr::mul(Expr::mul(n, k), t)),
    ))
}

/// Van der Waals equation: (P + a/V²)(V - b) = NkT
///
/// Returns symbolic representation
pub fn van_der_waals_equation() -> Expr {
    let p = Expr::sym("P");
    let v = Expr::sym("V");
    let a = Expr::sym("a");
    let b = Expr::sym("b");
    let n = Expr::sym("N");
    let k = Expr::sym("k_B");
    let t = Expr::sym("T");

    // (P + a/V²)
    let pressure_term = Expr::add(p, Expr::mul(a, Expr::pow(v.clone(), Expr::num(-2))));

    // (V - b)
    let volume_term = Expr::add(v, Expr::mul(Expr::num(-1), b));

    // LHS = (P + a/V²)(V - b)
    let lhs = Expr::mul(pressure_term, volume_term);

    // RHS = NkT
    let rhs = Expr::mul(Expr::mul(n, k), t);

    // Equation: LHS - RHS = 0
    simplify(&Expr::add(lhs, Expr::mul(Expr::num(-1), rhs)))
}

/// Heat capacity at constant volume: C_v = (∂U/∂T)_V
///
/// Symbolic representation
pub fn heat_capacity_constant_volume() -> Expr {
    Expr::func("(∂U/∂T)_V", vec![])
}

/// Heat capacity at constant pressure: C_p = (∂H/∂T)_P
///
/// Symbolic representation
pub fn heat_capacity_constant_pressure() -> Expr {
    Expr::func("(∂H/∂T)_P", vec![])
}

/// Carnot efficiency: η = 1 - T_c/T_h
pub fn carnot_efficiency(cold_temp: &str, hot_temp: &str) -> Expr {
    let t_c = Expr::sym(cold_temp);
    let t_h = Expr::sym(hot_temp);

    // η = 1 - T_c/T_h
    simplify(&Expr::add(
        Expr::num(1),
        Expr::mul(Expr::num(-1), Expr::mul(t_c, Expr::pow(t_h, Expr::num(-1)))),
    ))
}

