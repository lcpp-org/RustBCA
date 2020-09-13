use super::*;

const LENNARD_JONES_EPSILON: f64 = 0.343*EV;
const LENNARD_JONES_SIGMA: f64 = 1.*ANGSTROM;

pub fn interaction_potential(r: f64, a: f64, Za: f64, Zb: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE | KR_C | LENZ_JENSEN | ZBL => {
            screened_coulomb(r, a, Za, Zb, interaction_potential)
        },
        LENNARD_JONES_12_6 => {
            let epsilon = LENNARD_JONES_EPSILON;
            let sigma = LENNARD_JONES_SIGMA;
            lennard_jones(r, sigma, epsilon)
        },
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

pub fn energy_threshold_single_root(interaction_potential: i32) -> f64 {
    match interaction_potential{
        LENNARD_JONES_12_6 => 4./5.*LENNARD_JONES_EPSILON*1E9,
        MOLIERE | KR_C | LENZ_JENSEN | ZBL => 0.,
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

fn doca_function_transformed(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    return x0*x0 - x0*interactions::phi(x0, interaction_potential)/reduced_energy - beta*beta;
}

fn doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //Transcendental function to _ distance of closest approach
    return x0 - interactions::phi(x0, interaction_potential)/reduced_energy - beta*beta/x0;
}

fn diff_doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //First differential of distance of closest approach function for N-R solver
    return beta*beta/x0/x0 - interactions::dphi(x0, interaction_potential)/reduced_energy + 1.
}

fn diff_doca_function_transformed(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //First differential of distance of closest approach function for N-R solver
    return 2.*x0 - interactions::phi(x0, interaction_potential)/reduced_energy
}

pub fn distance_of_closest_approach_function(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: i32) -> f64 {
    if r.is_nan() {
        panic!("r is nan")
    }
    match interaction_potential {
        MOLIERE | KR_C | LENZ_JENSEN | ZBL => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            doca_function(r/a, beta, reduced_energy, interaction_potential)

        },
        LENNARD_JONES_12_6 => {
            let epsilon = LENNARD_JONES_EPSILON;
            let sigma = LENNARD_JONES_SIGMA;
            doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        }
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

pub fn distance_of_closest_approach_function_singularity_free(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: i32) -> f64 {
    if r.is_nan() {
        panic!("r is nan")
    }
    match interaction_potential {
        MOLIERE | KR_C | LENZ_JENSEN | ZBL => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            doca_function_transformed(r/a, beta, reduced_energy, interaction_potential)
        },
        LENNARD_JONES_12_6 => {
            let epsilon = LENNARD_JONES_EPSILON;
            let sigma = LENNARD_JONES_SIGMA;
            doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        }
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

pub fn scaling_function(r: f64, a: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE | KR_C | LENZ_JENSEN | ZBL => {
            1./(a + r)
        },
        LENNARD_JONES_12_6 => {
            let n = 10.;
            1./(a.powf(n) + r.powf(n))
        },
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

pub fn diff_distance_of_closest_approach_function(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE | KR_C | LENZ_JENSEN |ZBL => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            //let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            diff_doca_function(r/a, beta, reduced_energy, interaction_potential)
        },
        LENNARD_JONES_12_6 => {
            let epsilon = LENNARD_JONES_EPSILON;
            let sigma = LENNARD_JONES_SIGMA;
            diff_doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

pub fn diff_distance_of_closest_approach_function_singularity_free(r: f64, a: f64, Za: f64, Zb: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE | KR_C | LENZ_JENSEN |ZBL => {
            let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
            //let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a/Za/Zb*relative_energy;
            let beta: f64 = impact_parameter/a;
            diff_doca_function_transformed(r/a, beta, reduced_energy, interaction_potential)
        },
        LENNARD_JONES_12_6 => {
            let epsilon = LENNARD_JONES_EPSILON;
            let sigma = LENNARD_JONES_SIGMA;
            diff_doca_lennard_jones(r, impact_parameter, relative_energy, sigma, epsilon)
        },
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}

pub fn screened_coulomb(r: f64, a: f64, Za: f64, Zb: f64, interaction_potential: i32) -> f64 {
    Za*Zb*Q*Q/4./PI/EPS0/r*phi(r/a, interaction_potential)
}

pub fn phi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => moliere(xi),
        KR_C => kr_c(xi),
        ZBL => zbl(xi),
        LENZ_JENSEN => lenz_jensen(xi),
        TRIDYN => kr_c(xi),
        _ => panic!("Input error: unimplemented screened coulomb potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN", interaction_potential)
    }
}

pub fn dphi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => diff_moliere(xi),
        KR_C => diff_kr_c(xi),
        ZBL => diff_zbl(xi),
        LENZ_JENSEN => diff_lenz_jensen(xi),
        TRIDYN => diff_kr_c(xi),
        _ => panic!("Input error: unimplemented screened coulomb potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN", interaction_potential)
    }
}

pub fn screening_length(Za: f64, Zb: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        //ZBL screening length, Eckstein (4.1.8)
        ZBL => 0.88534*A0/(Za.powf(0.23) + Zb.powf(0.23)),
        //Lindhard/Firsov screening length, Eckstein (4.1.5)
        _ => 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.)
    }
}

pub fn lennard_jones(r: f64, sigma: f64, epsilon: f64) -> f64 {
    4.*epsilon*((sigma/r).powf(12.) - (sigma/r).powf(6.))
}

pub fn doca_lennard_jones(r: f64, p: f64, relative_energy: f64, sigma: f64, epsilon: f64) -> f64 {
    r.powf(12.) - 4.*epsilon/relative_energy*(sigma.powf(12.) - sigma.powf(6.)*r.powf(6.)) - p*p*r.powf(10.)
}

pub fn diff_doca_lennard_jones(r: f64, p: f64, relative_energy: f64, sigma: f64, epsilon: f64) -> f64 {
    12.*r.powf(11.) + 4.*epsilon*sigma.powf(6.)*6.*r.powf(5.) - p*p*10.*r.powf(9.)
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

pub fn first_screening_radius(interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => 0.3,
        KR_C => 0.278544,
        ZBL => 0.20162,
        LENZ_JENSEN => 0.206,
        TRIDYN => 0.278544,
        LENNARD_JONES_12_6 => LENNARD_JONES_SIGMA,
        _ => panic!("Input error: unimplemented interaction potential: {}, Use 0: MOLIERE 1: KR_C 2: ZBL 3: LENZ_JENSEN 4: LENNARD_JONES_12_6", interaction_potential)
    }
}
