//! This crate contains approximations for a set of mathematical functions 
//! commonly referred to as "special functions". 

pub mod cephes64;

pub struct EllipJOutput<T> {
    pub sn: T,
    pub cn: T,
    pub dn: T,
    pub phi: T
}

/// Implementations of elliptic integrals as a trait
pub trait Ellip: Sized {
    /// Jacobian elliptic functions
    fn ellip_j(&self, m: Self) -> EllipJOutput<Self>;

    /// Complete elliptic integral of the second kind
    fn ellip_e(&self) -> Self;

    /// Complete elliptic integral of the first kind
    fn ellip_k(&self) -> Self;

    /// Incomplete elliptic integral of the second kind
    fn ellip_e_inc(&self, m: Self) -> Self;

    /// Incomplete elliptic integral of the first kind
    fn ellip_k_inc(&self, m: Self) -> Self;
}

impl Ellip for f64 {
    fn ellip_j(&self, m: f64) -> EllipJOutput<f64> {
        //! Uses [`cephes64::ellpj`](crate::cephes64::ellpj)
        let res = crate::cephes64::ellpj(*self, m);
        EllipJOutput::<f64> {sn: res.0, cn: res.1, dn: res.2, phi: res.3}
    }
    fn ellip_e(&self) -> f64 {
        //! Uses [`cephes64::ellpe`](crate::cephes64::ellpe) evaluated at `1.0 - self`
        crate::cephes64::ellpe(1.0 - *self)
    }
    fn ellip_k(&self) -> f64 {
        //! Uses [`cephes64::ellpk`](crate::cephes64::ellpk) evaluated at `1.0 - self`
        crate::cephes64::ellpk(1.0 - *self)
    }
    fn ellip_e_inc(&self, m: f64) -> f64 {
        //! Uses [`cephes64::ellie`](crate::cephes64::ellie)
        crate::cephes64::ellie(*self, m)
    }
    fn ellip_k_inc(&self, m: f64) -> f64 {
        //! Uses [`cephes64::ellik`](crate::cephes64::ellik)
        crate::cephes64::ellik(*self, m)
    }
}

/// Implementations of Bessel functions as a trait
pub trait Bessel {
    /// Bessel function of the first kind, order zero
    fn bessel_j0(&self) -> Self;

    /// Bessel function of the second kind, order zero
    fn bessel_y0(&self) -> Self;

    /// Bessel function of the first kind, order one
    fn bessel_j1(&self) -> Self;

    /// Bessel function of the second kind, order one
    fn bessel_y1(&self) -> Self;
}

impl Bessel for f64 {
    fn bessel_j0(&self) -> f64 {
        //! Uses [`cephes64::j0`](crate::cephes64::j0)
        crate::cephes64::j0(*self)
    }
    fn bessel_y0(&self) -> f64 {
        //! Uses [`cephes64::y0`](crate::cephes64::y0)
        crate::cephes64::y0(*self)
    }
    fn bessel_j1(&self) -> f64 {
        //! Uses [`cephes64::j1`](crate::cephes64::j1)
        crate::cephes64::j1(*self)
    }
    fn bessel_y1(&self) -> f64 {
        //! Uses [`cephes64::y1`](crate::cephes64::y1)
        crate::cephes64::y1(*self)
    }
}