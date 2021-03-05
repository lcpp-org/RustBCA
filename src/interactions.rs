use super::*;

/// Analytic solutions to outermost root of the interaction potential.
pub fn crossing_point_doca(interaction_potential: InteractionPotential) -> f64 {

    match interaction_potential {
        InteractionPotential::LENNARD_JONES_12_6{sigma, ..} | InteractionPotential::LENNARD_JONES_65_6{sigma, ..} => sigma,
        InteractionPotential::MORSE{D, alpha, r0} => (alpha*r0 - (2.0_f64).ln())/alpha,
        InteractionPotential::WW => 50.*ANGSTROM,
        _ => 10.*ANGSTROM,
        //_ => panic!("Input error: potential never crosses zero for r > 0. Consider using the Newton rootfinder.")
    }

}

/// Interaction potential between two particles a and b at a distance `r`.
pub fn interaction_potential(r: f64, a: f64, Za: f64, Zb: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::ZBL | InteractionPotential::TRIDYN => {
            screened_coulomb(r, a, Za, Zb, interaction_potential)
        },
        InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => {
            lennard_jones(r, sigma, epsilon)
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => {
            lennard_jones_65_6(r, sigma, epsilon)
        }
        InteractionPotential::MORSE{D, alpha, r0} => {
            morse(r, D, alpha, r0)
        }
        InteractionPotential::WW => {
            tungsten_tungsten_cubic_spline(r)
        }
        InteractionPotential::COULOMB{Za, Zb} => {
            coulomb(r, Za, Zb)
        }
    }
}

/// Analytic energy threshold above which the distance of closest approach function is guaranteed to have only one root.
pub fn energy_threshold_single_root(interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential{
        InteractionPotential::LENNARD_JONES_12_6{..} | InteractionPotential::LENNARD_JONES_65_6{..} => f64::INFINITY,
        InteractionPotential::MORSE{..} => f64::INFINITY,
        InteractionPotential::WW => f64::INFINITY,
        InteractionPotential::COULOMB{..} => f64::INFINITY,
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::ZBL | InteractionPotential::TRIDYN => 0.,
    }
}

/// Distance of closest approach function for screened coulomb potentials.
fn doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: InteractionPotential) -> f64 {
    //Nonlinear equation to determine distance of closest approach
    return x0 - interactions::phi(x0, interaction_potential)/reduced_energy - beta*beta/x0;
}

/// First derivative w.r.t. `r` of the distance of closest approach function for screened coulomb potentials.
fn diff_doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: InteractionPotential) -> f64 {
    //First differential of distance of closest approach function for N-R solver
    return beta*beta/x0/x0 - interactions::dphi(x0, interaction_potential)/reduced_energy + 1.
}

/// Singularity-free version of the screened-coulomb distance of closest approach function.
fn doca_function_transformed(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: InteractionPotential) -> f64 {
    //Singularity free version of doca function
    return x0*x0 - x0*interactions::phi(x0, interaction_potential)/reduced_energy - beta*beta;
}

/// First derivative w.r.t. `r` of the singularity-free version of the screened-coulomb distance of closest approach function.
fn diff_doca_function_transformed(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: InteractionPotential) -> f64 {
    //First differential of distance of closest approach function for N-R solver
    return 2.*x0 - interactions::phi(x0, interaction_potential)/reduced_energy
}

/// Distance of closest approach function. The outermost root of this function is the distance of closest approach, or classical turning point, which is the lower bound to the scattering integral.
pub fn distance_of_closest_approach_function(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::ZBL | InteractionPotential::TRIDYN => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            doca_function(r/a, beta, reduced_energy, interaction_potential)

        },
        InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => {
            doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => {
            doca_lennard_jones_65_6(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::MORSE{D, alpha, r0} => {
            doca_morse(r, impact_parameter, relative_energy, D, alpha, r0)
        },
        InteractionPotential::WW => {
            doca_tungsten_tungsten_cubic_spline(r, impact_parameter, relative_energy)
        },
        InteractionPotential::COULOMB{..} => panic!("Coulombic potential cannot be used with rootfinder.")
    }
}

/// Singularity-free distance of closest approach function. The outermost root of this function is the distance of closest approach, or classical turning point, which is the lower bound to the scattering integral.
pub fn distance_of_closest_approach_function_singularity_free(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: InteractionPotential) -> f64 {
    if r.is_nan() {
        panic!("Numerical error: r is NaN in distance of closest approach function. Check Rootfinder.")
    }
    match interaction_potential {
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::ZBL | InteractionPotential::TRIDYN => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            doca_function_transformed(r/a, beta, reduced_energy, interaction_potential)
        },
        InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => {
            doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => {
            doca_lennard_jones_65_6(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::MORSE{D, alpha, r0} => {
            doca_morse(r, impact_parameter, relative_energy, D, alpha, r0)
        },
        InteractionPotential::WW => {
            doca_tungsten_tungsten_cubic_spline(r, impact_parameter, relative_energy)
        },
        InteractionPotential::COULOMB{..} => panic!("Coulombic potential cannot be used with rootfinder.")
    }
}

/// Scaling function used to keep the distance of closest approach function ~1 for the CPR rootfinder.
pub fn scaling_function(r: f64, a: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::ZBL | InteractionPotential::TRIDYN => {
            1./(1. + (r/a).powi(2))
        },
        InteractionPotential::LENNARD_JONES_12_6{sigma, ..} => {
            let n = 11.;
            1./(1. + (r/sigma).powf(n))
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, ..} => {
            let n = 6.;
            1./(1. + (r/sigma).powf(n))
        },
        InteractionPotential::MORSE{D, alpha, r0} => {
            1./(1. + (r*alpha).powi(2))
        }
        InteractionPotential::WW => {
            1.
        },
        InteractionPotential::COULOMB{..} => panic!("Coulombic potential cannot be used with rootfinder.")
    }
}

/// First derivative w.r.t. `r` of the distance of closest approach function. The outermost root of this function is the distance of closest approach, or classical turning point, which is the lower bound to the scattering integral.
pub fn diff_distance_of_closest_approach_function(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN |InteractionPotential::ZBL | InteractionPotential::TRIDYN => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            //let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            diff_doca_function(r/a, beta, reduced_energy, interaction_potential)
        },
        InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => {
            diff_doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => {
            diff_doca_lennard_jones_65_6(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::MORSE{D, alpha, r0} => {
            diff_doca_morse(r, impact_parameter, relative_energy, D, alpha, r0)
        },
        _ => panic!("Input error: {} does not have an implemented derivative. Try using the derivative free CPR-rootfinder.", interaction_potential)
    }
}

/// First derivative w.r.t. `r` of the singularity-free distance of closest approach function. The outermost root of this function is the distance of closest approach, or classical turning point, which is the lower bound to the scattering integral.
pub fn diff_distance_of_closest_approach_function_singularity_free(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::ZBL | InteractionPotential::TRIDYN => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            //let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            diff_doca_function_transformed(r/a, beta, reduced_energy, interaction_potential)
        },
        InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => {
            diff_doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => {
            diff_doca_lennard_jones_65_6(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        InteractionPotential::MORSE{D, alpha, r0} => {
            diff_doca_morse(r, impact_parameter, relative_energy, D, alpha, r0)
        },
        _ => panic!("Input error: {} does not have an implemented derivative. Try using the derivative free CPR-rootfinder.", interaction_potential)
    }
}

/// Screened coulomb interaction potential.
pub fn screened_coulomb(r: f64, a: f64, Za: f64, Zb: f64, interaction_potential: InteractionPotential) -> f64 {
    Za*Zb*Q*Q/4./PI/EPS0/r*phi(r/a, interaction_potential)
}

/// Coulombic interaction potential.
pub fn coulomb(r: f64, Za: f64, Zb: f64) -> f64 {
    return Za*Zb*Q*Q/4./PI/EPS0/r;
}

/// Screening functions for screened-coulomb interaction potentials.
pub fn phi(xi: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE => moliere(xi),
        InteractionPotential::KR_C => kr_c(xi),
        InteractionPotential::ZBL => zbl(xi),
        InteractionPotential::LENZ_JENSEN => lenz_jensen(xi),
        InteractionPotential::TRIDYN => kr_c(xi),
        _ => panic!("Input error: Screened Coulomb cannot be used with {}", interaction_potential)
    }
}

/// First derivative w.r.t. `r` of the screening functions for screened-coulomb interaction potentials.
pub fn dphi(xi: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE => diff_moliere(xi),
        InteractionPotential::KR_C => diff_kr_c(xi),
        InteractionPotential::ZBL => diff_zbl(xi),
        InteractionPotential::LENZ_JENSEN => diff_lenz_jensen(xi),
        InteractionPotential::TRIDYN => diff_kr_c(xi),
        _ => panic!("Input error: Screened Coulomb cannot be used with {}", interaction_potential)
    }
}

/// Screening length of interatomic potentials.
pub fn screening_length(Za: f64, Zb: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        //ZBL screening length, Eckstein (4.1.8)
        InteractionPotential::ZBL => 0.88534*A0/(Za.powf(0.23) + Zb.powf(0.23)),
        //Lindhard/Firsov screening length, Eckstein (4.1.5)
        InteractionPotential::MOLIERE | InteractionPotential::KR_C | InteractionPotential::LENZ_JENSEN | InteractionPotential::TRIDYN | InteractionPotential::WW => 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.),
        InteractionPotential::LENNARD_JONES_12_6{..} | InteractionPotential::LENNARD_JONES_65_6{..} => 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.),
        InteractionPotential::MORSE{D, alpha, r0} => alpha,
        InteractionPotential::COULOMB{Za: Z1, Zb: Z2} => 0.88534*A0/(Z1.powf(0.23) + Z2.powf(0.23)),
    }
}

/// Coefficients of inverse-polynomial interaction potentials.
pub fn polynomial_coefficients(relative_energy: f64, impact_parameter: f64, interaction_potential: InteractionPotential) -> Vec<f64> {
    match interaction_potential {
        InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => {
            vec![1., 0., -impact_parameter.powi(2), 0., 0., 0., 4.*epsilon*sigma.powf(6.)/relative_energy, 0., 0., 0., 0., 0., -4.*epsilon*sigma.powf(12.)/relative_energy]
        },
        InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => {
            vec![1., 0., 0., 0., -impact_parameter.powi(2), 0., 0., 0., 0., 0., 0., 0., 4.*epsilon*sigma.powf(6.)/relative_energy, -4.*epsilon*sigma.powf(6.5)/relative_energy]
        },
        _ => panic!("Input error: non-polynomial interaction potential used with polynomial root-finder.")
    }
}

/// Inverse-polynomial interaction potentials transformation to remove singularities for the CPR root-finder.
pub fn inverse_transform(x: f64, interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::LENNARD_JONES_12_6{..} => {
            x
        },
        InteractionPotential::LENNARD_JONES_65_6{..} => {
            x*x
        },
        _ => panic!("Input error: non-polynomial interaction potential used with polynomial root-finder transformation.")
    }
}

/// Lennard-Jones 12-6
pub fn lennard_jones(r: f64, sigma: f64, epsilon: f64) -> f64 {
    4.*epsilon*((sigma/r).powf(12.) - (sigma/r).powf(6.))
}

/// Lennard-Jones 6.5-6
pub fn lennard_jones_65_6(r: f64, sigma: f64, epsilon: f64) -> f64 {
    4.*epsilon*((sigma/r).powf(6.5) - (sigma/r).powf(6.))
}

/// Morse potential
pub fn morse(r: f64, D: f64, alpha: f64, r0: f64) -> f64 {
    D*((-2.*alpha*(r - r0)).exp() - 2.*(-alpha*(r - r0)).exp())
}

/// Distance of closest approach function for Morse potential.
pub fn doca_morse(r: f64, impact_parameter: f64, relative_energy: f64, D: f64, alpha: f64, r0: f64) -> f64 {
    (r*alpha).powi(2) - (r*alpha).powi(2)*D/relative_energy*((-2.*alpha*(r - r0)).exp() - 2.*(-alpha*(r - r0)).exp()) - (impact_parameter*alpha).powi(2)
}

/// First derivative w.r.t. `r` of the distance of closest approach function for Morse potential.
pub fn diff_doca_morse(r: f64, impact_parameter: f64, relative_energy: f64, D: f64, alpha: f64, r0: f64) -> f64 {
    2.*alpha.powi(2)*r - 2.*alpha.powi(2)*D*r*(-2.*alpha*(r - r0) - 1.).exp()*(alpha*r*(alpha*(r - r0)).exp() - 2.*(alpha*(r - r0)).exp() - r*alpha + 1.)
}

/// Distance of closest approach function for LJ 6.5-6 potential.
pub fn doca_lennard_jones_65_6(r: f64, p: f64, relative_energy: f64, sigma: f64, epsilon: f64) -> f64 {
    (r/sigma).powf(6.5) - 4.*epsilon/relative_energy*(1. - (r/sigma).powf(0.5)) - (p/sigma).powi(2)*(r/sigma).powf(4.5)
}

/// Distance of closest approach function for LJ 12-6 potential.
pub fn doca_lennard_jones(r: f64, p: f64, relative_energy: f64, sigma: f64, epsilon: f64) -> f64 {
    (r/sigma).powf(12.) - 4.*epsilon/relative_energy*(1. - (r/sigma).powf(6.)) - p.powi(2)*r.powf(10.)/sigma.powf(12.)
}

/// First derivative w.r.t. `r` of the distance of closest approach function for LJ 12-6 potential.
pub fn diff_doca_lennard_jones(r: f64, p: f64, relative_energy: f64, sigma: f64, epsilon: f64) -> f64 {
    12.*(r/sigma).powf(11.)/sigma + 4.*epsilon/relative_energy*6.*(r/sigma).powf(5.)/sigma - 10.*p.powi(2)*r.powf(9.)/sigma.powf(12.)
}

/// First derivative w.r.t. `r` of the distance of closest approach function for LJ 6.5-6 potential.
pub fn diff_doca_lennard_jones_65_6(r: f64, p: f64, relative_energy: f64, sigma: f64, epsilon: f64) -> f64 {
    6.5*(r/sigma).powf(5.5)/sigma + 4.*epsilon/relative_energy*0.5*(sigma*r).powf(-0.5) - (p/sigma).powi(2)*4.5*(r/sigma).powf(3.5)/sigma
}

/// W-W cublic spline potential from Bjorkas et al.
pub fn tungsten_tungsten_cubic_spline(r: f64) -> f64 {

    let x = r/ANGSTROM;
    let x1 = 1.10002200044;
    let x2 = 2.10004200084;

    if x <= x1 {

        let a = screening_length(74., 74., InteractionPotential::ZBL);
        screened_coulomb(r, a, 74., 74., InteractionPotential::ZBL)

    } else if x <= x2 {

        let a = vec![
            1.389653276380862E4,
            -3.596912431628216E4,
            3.739206756369099E4,
            -1.933748081656593E4,
            0.495516793802426E4,
            -0.050264585985867E4
        ];

        (a[0] + a[1]*x + a[2]*x.powi(2) + a[3]*x.powi(3) + a[4]*x.powf(4.) + a[5]*x.powf(5.))*EV

    } else {

        let a = vec![
            -0.1036435865158945,
            -0.2912948318493851,
            -2.096765499656263,
            19.16045452701010,
            -41.01619862085917,
            46.05205617244703,
            26.42203930654883,
            15.35211507804088,
            14.12806259323987,
        ];

        let delta = vec![
            4.268900000000000,
            3.985680000000000,
            3.702460000000000,
            3.419240000000000,
            3.136020000000000,
            2.852800000000000,
            2.741100000000000,
            2.604045000000000,
            2.466990000000000,
        ];

        a.iter().zip(delta).map(|(&a_i, delta_i)| EV*a_i*(delta_i - x).powi(3)*heaviside(delta_i - x)).sum::<f64>()
    }
}

/// Distance of closest approach function for the W-W cublic spline potential from Bjorkas et al.
pub fn doca_tungsten_tungsten_cubic_spline(r: f64, p: f64, relative_energy: f64) -> f64 {

    let x = r/ANGSTROM;
    let x1 = 1.10002200044;
    let x2 = 2.10004200084;

    if x <= x1 {
        let a = screening_length(74., 74., InteractionPotential::ZBL);
        distance_of_closest_approach_function_singularity_free(r, a, 74., 74., relative_energy, p, InteractionPotential::ZBL)
    } else {
        (r/ANGSTROM).powi(2) - (r/ANGSTROM).powi(2)*tungsten_tungsten_cubic_spline(r)/relative_energy - p.powi(2)/ANGSTROM.powi(2)
    }
}

fn heaviside(x: f64) -> f64 {
    if x >= 0. {
        1.
    } else {
        0.
    }
}

fn moliere(xi: f64) -> f64 {
    0.35*(-0.3*xi).exp() + 0.55*(-1.2*xi).exp() + 0.10*(-6.0*xi).exp()
}

fn kr_c(xi: f64) -> f64 {
    0.190945*(-0.278544*xi).exp() + 0.473674*(-0.637174*xi).exp() + 0.335381*(-1.919249*xi).exp()
}

fn zbl(xi: f64) -> f64 {
    0.02817*(-0.20162*xi).exp() + 0.28022*(-0.40290*xi).exp() + 0.50986*(-0.94229*xi).exp() + 0.18175*(-3.1998*xi).exp()
}

fn lenz_jensen(xi: f64) -> f64 {
    0.01018*(-0.206*xi).exp() + 0.24330*(-0.3876*xi).exp() + 0.7466*(-1.038*xi).exp()
}

fn diff_moliere(xi: f64) -> f64 {
    -0.35*0.3*(-0.3*xi).exp() -0.55*1.2*(-1.2*xi).exp() -0.10*6.0*(-6.0*xi).exp()
}

fn diff_kr_c(xi: f64) -> f64 {
    -0.278544*0.190945*(-0.278544*xi).exp() - 0.637174*0.473674*(-0.637174*xi).exp() - 0.335381*1.919249*(-1.919249*xi).exp()
}

fn diff_zbl(xi: f64) -> f64 {
    -0.20162*0.02817*(-0.20162*xi).exp() -0.40290*0.28022*(-0.40290*xi).exp() -0.94229*0.50986*(-0.94229*xi).exp() -3.1998*0.18175*(-3.1998*xi).exp()
}

fn diff_lenz_jensen(xi: f64) -> f64 {
     -0.206*0.01018*(-0.206*xi).exp() -0.3876*0.24330*(-0.3876*xi).exp() -1.038*0.7466*(-1.038*xi).exp()
}

/// Normalized first screening radius. Used in Oen-Robinson electronic stopping.
pub fn first_screening_radius(interaction_potential: InteractionPotential) -> f64 {
    match interaction_potential {
        InteractionPotential::MOLIERE => 0.3,
        InteractionPotential::KR_C => 0.278544,
        InteractionPotential::ZBL => 0.20162,
        InteractionPotential::LENZ_JENSEN => 0.206,
        InteractionPotential::TRIDYN => 0.278544,
        InteractionPotential::LENNARD_JONES_12_6{..} => 1.,
        InteractionPotential::LENNARD_JONES_65_6{..} => 1.,
        InteractionPotential::MORSE{..} => 1.,
        InteractionPotential::WW => 1.,
        InteractionPotential::COULOMB{..} => 0.3,
    }
}
