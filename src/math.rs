use super::*;
/// This function finds an orthonormal basis from a unit vector n
/// it should avoid all numerical / consistency issues
/// Duff et al., JCGT 2017
/// http://jcgt.org/published/0006/01/01/
pub fn duff_orthonormal_basis(n: Vector) -> (Vector, Vector) {
    if n.z < 0. {
        let a = 1.0 / (1.0 - n.z);
        let b = n.x * n.y * a;
        let b1 = Vector::new(1.0 - n.x*n.x*a, -b, n.x);
        let b2 = Vector::new(b, n.y * n.y*a- 1.0, -n.y);
        (b1, b2)
    } else {
        let a = 1.0 / (1.0 + n.z);
        let b = -n.x * n.y * a;
        let b1 = Vector::new(1.0 - n.x*n.x*a, b, -n.x);
        let b2 = Vector::new(b, 1.0 - n.y*n.y*a, -n.y);
        (b1, b2)
    }
}