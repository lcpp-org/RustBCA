/// 3D vector.
#[derive(Clone)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Vector {
        Vector {
            x,
            y,
            z
        }
    }

    /// Calculates vector magnitude.
    pub fn magnitude(&self) -> f64 {
        (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    /// Assigns vector values from another vector.
    pub fn assign(&mut self, other: &Vector) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }

    /// Normalizes vector components to magnitude 1.
    pub fn normalize(&mut self) {
        let magnitude = self.magnitude();
        self.x /= magnitude;
        self.y /= magnitude;
        self.z /= magnitude;
    }

    /// Add this vector and another and return a new vector.
    pub fn add(&self, other: &Vector) -> Vector {
        Vector::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    pub fn dot(&self, other: &Vector) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }
}

/// Vector4 is a trajectory-tracking object that includes x, y, z, and the current energy.
#[derive(Clone)]
pub struct Vector4 {
    pub E: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector4 {
    pub fn new(E: f64, x: f64, y: f64, z: f64) -> Vector4 {
        Vector4 {
            E,
            x,
            y,
            z
        }
    }
}

/// Energy loss is an output tracker that tracks the separate nuclear and electronic energy losses.
#[derive(Clone)]
pub struct EnergyLoss {
    pub En: f64,
    pub Ee: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl EnergyLoss {
    pub fn new(Ee: f64, En: f64, x: f64, y: f64, z: f64) -> EnergyLoss {
        EnergyLoss {
            En,
            Ee,
            x,
            y,
            z
        }
    }
}
