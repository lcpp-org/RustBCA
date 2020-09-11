use super::*;

pub const LJ: &dyn Fn(f64) -> f64 = &lennard_jones;

pub const PHI_MOLIERE: &dyn Fn(f64) -> f64 = &moliere;
pub const PHI_KR_C: &dyn Fn(f64) -> f64 = &kr_c;
pub const PHI_LENZ_JENSEN: &dyn Fn(f64) -> f64 = &lenz_jensen;
pub const PHI_ZBL: &dyn Fn(f64) -> f64 = &zbl;

pub const DPHI_MOLIERE: &dyn Fn(f64) -> f64 = &diff_moliere;
pub const DPHI_KR_C: &dyn Fn(f64) -> f64 = &diff_kr_c;
pub const DPHI_LENZ_JENSEN: &dyn Fn(f64) -> f64 = &diff_lenz_jensen;
pub const DPHI_ZBL: &dyn Fn(f64) -> f64 = &diff_zbl;

pub fn screened_coulomb(r: f64, a: f64, Za: f64, Zb: f64, screening_function: &dyn Fn(f64) -> f64) -> f64 {
    Za*Zb*Q*Q/4./PI/EPS0/r*screening_function(r/a)
}

pub fn doca_function(r: f64, relative_energy: f64, impact_parameter: f64, interaction_potential: i32) -> f64 {

}

pub fn phi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => moliere(xi),
        KR_C => kr_c(xi),
        ZBL => zbl(xi),
        LENZ_JENSEN => lenz_jensen(xi),
        TRIDYN => kr_c(xi),
        _ => panic!("Unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
    }
}

pub fn dphi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => diff_moliere(xi),
        KR_C => diff_kr_c(xi),
        ZBL => diff_zbl(xi),
        LENZ_JENSEN => diff_lenz_jensen(xi),
        TRIDYN => diff_kr_c(xi),
        _ => panic!("Unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
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

pub fn lennard_jones(r: f64) -> f64 {
    let epsilon = 0.343*EV;
    let sigma = 2.*ANGSTROM;

    4.*epsilon*((sigma/r).powf(12.) - (sigma/r).powf(6.))
}

pub fn doca_lennard_jones(r: f64, p: f64, relative_energy: f64) -> f64 {
    let epsilon = 0.343*EV;
    let sigma = 2.*ANGSTROM;

    r.powf(12.) - 4.*epsilon/relative_energy*(sigma.powf(12.) - sigma.powf(6.)*r.powf(6.)) - p*p*r.powf(10.)
}

pub fn diff_doca_lennard_jones(r: f64, p: f64, relative_energy: f64) -> f64 {

    let epsilon = 0.343*EV;
    let sigma = 2.*ANGSTROM;

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
        _ => panic!("Input error: unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
    }
}
