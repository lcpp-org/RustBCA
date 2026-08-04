use super::*;
/// This function finds an orthonormal basis from a unit vector n
/// it should avoid all numerical / consistency issues
/// Duff et al., JCGT 2017
/// http://jcgt.org/published/0006/01/01/
pub fn duff_orthonormal_basis(n: Vector) -> (Vector, Vector) {
    let sign = (1.0_f64).copysign(n.z);
    let a = -1.0_f64 / (sign + n.z);
    let b = n.x*n.y*a;
    let b1 = Vector::new(1.0 + sign*n.x*n.x*a, sign*b, -sign*n.x);
    let b2 = Vector::new(b, sign + n.y*n.y*a, -n.y);
    (b1, b2)
}

pub fn triangular_index(i: &mut usize, j: &mut usize) -> usize {
    if i < j {
        std::mem::swap(i, j);
    }
    (*i*(*i + 1)/2) + *j
}