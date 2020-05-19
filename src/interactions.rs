use super::*;

pub fn phi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => 0.35*(-0.3*xi).exp() + 0.55*(-1.2*xi).exp() + 0.10*(-6.0*xi).exp(),
        KR_C => 0.190945*(-0.278544*xi).exp() + 0.473674*(-0.637174*xi).exp() + 0.335381*(-1.919249*xi).exp(),
        ZBL => 0.02817*(-0.20162*xi).exp() + 0.28022*(-0.40290*xi).exp() + 0.50986*(-0.94229*xi).exp() + 0.18175*(-3.1998*xi).exp(),
        LENZ_JENSEN => 0.01018*(-0.206*xi).exp() + 0.24330*(-0.3876*xi).exp() + 0.7466*(-1.038*xi).exp(),
        TRIDYN => 0.190945*(-0.278544*xi).exp() + 0.473674*(-0.637174*xi).exp() + 0.335381*(-1.919249*xi).exp(),
        _ => panic!("Unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
    }
}

pub fn dphi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => -0.35*0.3*(-0.3*xi).exp() + -0.55*1.2*(-1.2*xi).exp() + -0.10*6.0*(-6.0*xi).exp(),
        KR_C => -0.278544*0.190945*(-0.278544*xi).exp() - 0.637174*0.473674*(-0.637174*xi).exp() - 0.335381*1.919249*(-1.919249*xi).exp(),
        ZBL => -0.20162*0.02817*(-0.20162*xi).exp() + -0.40290*0.28022*(-0.40290*xi).exp() + -0.94229*0.50986*(-0.94229*xi).exp() + -3.1998*0.18175*(-3.1998*xi).exp(),
        LENZ_JENSEN => -0.206*0.01018*(-0.206*xi).exp() + -0.3876*0.24330*(-0.3876*xi).exp() + -1.038*0.7466*(-1.038*xi).exp(),
        TRIDYN => -0.278544*0.190945*(-0.278544*xi).exp() - 0.637174*0.473674*(-0.637174*xi).exp() - 0.335381*1.919249*(-1.919249*xi).exp(),
        _ => panic!("Unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
    }
}

pub fn screening_length(Za: f64, Zb: f64, interaction_potential: i32) -> f64 {
    println!("{} {}", 0.88534*A0/(Za.powf(0.23) + Zb.powf(0.23))/ANGSTROM, 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.)/ANGSTROM);
    match interaction_potential {
        //ZBL screening length, Eckstein (4.1.8)
        ZBL => 0.88534*A0/(Za.powf(0.23) + Zb.powf(0.23)),
        //Lindhard/Firsov screening length, Eckstein (4.1.5)
        _ => 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.)
    }
}
