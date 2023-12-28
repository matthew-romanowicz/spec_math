//! This crate contains approximations for a set of mathematical functions 
//! commonly referred to as "special functions". The goal of this crate is
//! to eventually contain a full Rust re-implementation of the CEPHES 
//! library with minimal changes in modules such as the [`cephes64`](crate::cephes64)
//! module. The special functions will also be available as traits for easier
//! use in Rust, potentially with modifications to correct bugs or improve performance
//! from the original CEPHES implementation.
//!
//! Currently, the error functions, gamma functions, elliptic integrals, and bessel functions are 
//! implemented.

pub mod cephes64;

pub struct EllipJOutput<T> {
    pub sn: T,
    pub cn: T,
    pub dn: T,
    pub phi: T
}

pub struct AiryOutput<T> {
    pub ai: T,
    pub aip: T,
    pub bi: T,
    pub bip: T
}

/// Implementations of error functions as a trait
pub trait Erf {
    /// Error function
    fn erf(&self) -> Self;

    /// Complimentary error function
    fn erfc(&self) -> Self;
}

impl Erf for f64 {
    fn erf(&self) -> f64 {
        //! Uses [`cephes64::erf`](crate::cephes64::erf)
        crate::cephes64::erf(*self)
    }
    fn erfc(&self) -> f64 {
        //! Uses [`cephes64::erfc`](crate::cephes64::erfc) 
        crate::cephes64::erfc(*self)
    }
}

/// Implementations of gamma and beta functions as a trait
pub trait Gamma {
    /// Gamma function
    fn gamma(&self) -> Self;

    /// Natural logarithm of the absolute value of the gamma function
    fn lgamma(&self) -> Self;

    /// Regularized lower incomplete gamma function
    fn igamma(&self, x: Self) -> Self;

    /// Regularized upper incomplete gamma function
    fn igammac(&self, x: Self) -> Self;

    /// Beta function
    fn beta(&self, other: Self) -> Self;

    /// Natural logarithm of the absolute value of the beta function
    fn lbeta(&self, other: Self) -> Self;
}

impl Gamma for f64 {
    fn gamma(&self) -> f64 {
        //! Uses [`cephes64::gamma`](crate::cephes64::gamma)
        crate::cephes64::gamma(*self)
    }
    fn lgamma(&self) -> f64 {
        //! Uses [`cephes64::lgam`](crate::cephes64::lgam) 
        crate::cephes64::lgam(*self)
    }
    fn igamma(&self, x: f64) -> f64 {
        //! Uses [`cephes64::igam`](crate::cephes64::igam) 
        crate::cephes64::igam(*self, x)
    }
    fn igammac(&self, x: f64) -> f64 {
        //! Uses [`cephes64::igamc`](crate::cephes64::igamc) 
        crate::cephes64::igamc(*self, x)
    }
    fn beta(&self, other: f64) -> f64 {
        //! Uses [`cephes64::beta`](crate::cephes64::beta)
        crate::cephes64::beta(*self, other)
    }
    fn lbeta(&self, other: f64) -> f64 {
        //! Uses [`cephes64::lbeta`](crate::cephes64::lbeta) 
        crate::cephes64::lbeta(*self, other)
    }
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
pub trait Bessel: Sized {
    /// Airy functions
    fn airy(&self) -> AiryOutput<Self>;

    /// Bessel function of the first kind, order zero
    fn bessel_j0(&self) -> Self;

    /// Bessel function of the second kind, order zero
    fn bessel_y0(&self) -> Self;

    /// Modified Bessel function of order zero
    fn bessel_i0(&self) -> Self;

    /// Modified Bessel function of order zero, exponentially scaled
    fn bessel_i0e(&self) -> Self;

    /// Bessel function of the first kind, order one
    fn bessel_j1(&self) -> Self;

    /// Bessel function of the second kind, order one
    fn bessel_y1(&self) -> Self;

    /// Modified Bessel function of order one
    fn bessel_i1(&self) -> Self;

    /// Modified Bessel function of order one, exponentially scaled
    fn bessel_i1e(&self) -> Self;

    /// Modified Bessel function of the third kind, order zero
    fn bessel_k0(&self) -> Self;

    /// Modified Bessel function of the third kind, order zero exponentially scaled
    fn bessel_k0e(&self) -> Self;

    /// Modified Bessel function of the third kind, order one
    fn bessel_k1(&self) -> Self;

    /// Modified Bessel function of the third kind, order one exponentially scaled
    fn bessel_k1e(&self) -> Self;

    /// Bessel function of second kind of integer order
    fn bessel_yn(&self, n: isize) -> Self;

    /// Modified Bessel function, third kind, integer order
    fn bessel_kn(&self, n: isize) -> Self;

    /// Bessel function of first kind of real order
    fn bessel_jv(&self, v: f64) -> Self;

    /// Bessel function of second kind of real order
    fn bessel_yv(&self, v: f64) -> Self;
}

impl Bessel for f64 {
    fn airy(&self) -> AiryOutput<f64> {
        //! Uses [`cephes64::airy`](crate::cephes64::airy)
        let res = crate::cephes64::airy(*self);
        AiryOutput::<f64> {ai: res.0, aip: res.1, bi: res.2, bip: res.3}
    }
    fn bessel_j0(&self) -> f64 {
        //! Uses [`cephes64::j0`](crate::cephes64::j0)
        crate::cephes64::j0(*self)
    }
    fn bessel_y0(&self) -> f64 {
        //! Uses [`cephes64::y0`](crate::cephes64::y0)
        crate::cephes64::y0(*self)
    }
    fn bessel_i0(&self) -> f64 {
        //! Uses [`cephes64::i0`](crate::cephes64::i0)
        crate::cephes64::i0(*self)
    }
    fn bessel_i0e(&self) -> f64 {
        //! Uses [`cephes64::i0e`](crate::cephes64::i0e)
        crate::cephes64::i0e(*self)
    }
    fn bessel_j1(&self) -> f64 {
        //! Uses [`cephes64::j1`](crate::cephes64::j1)
        crate::cephes64::j1(*self)
    }
    fn bessel_y1(&self) -> f64 {
        //! Uses [`cephes64::y1`](crate::cephes64::y1)
        crate::cephes64::y1(*self)
    }
    fn bessel_i1(&self) -> f64 {
        //! Uses [`cephes64::i1`](crate::cephes64::i1)
        crate::cephes64::i1(*self)
    }
    fn bessel_i1e(&self) -> f64 {
        //! Uses [`cephes64::i1e`](crate::cephes64::i1e)
        crate::cephes64::i1e(*self)
    }
    fn bessel_k0(&self) -> f64 {
        //! Uses [`cephes64::k0`](crate::cephes64::k0)
        crate::cephes64::k0(*self)
    }
    fn bessel_k0e(&self) -> f64 {
        //! Uses [`cephes64::k0e`](crate::cephes64::k0e)
        crate::cephes64::k0e(*self)
    }
    fn bessel_k1(&self) -> f64 {
        //! Uses [`cephes64::k1`](crate::cephes64::k1)
        crate::cephes64::k1(*self)
    }
    fn bessel_k1e(&self) -> f64 {
        //! Uses [`cephes64::k1e`](crate::cephes64::k1e)
        crate::cephes64::k1e(*self)
    }
    fn bessel_yn(&self, n: isize) -> f64 {
        //! Uses [`cephes64::yn`](crate::cephes64::yn)
        crate::cephes64::yn(n, *self)
    }
    fn bessel_kn(&self, n: isize) -> f64 {
        //! Uses [`cephes64::kn`](crate::cephes64::kn)
        crate::cephes64::kn(n, *self)
    }
    fn bessel_jv(&self, v: f64) -> f64 {
        //! Uses [`cephes64::jv`](crate::cephes64::jv)
        crate::cephes64::jv(v, *self)
    }
    fn bessel_yv(&self, v: f64) -> f64 {
        //! Uses [`cephes64::yv`](crate::cephes64::yv)
        crate::cephes64::yv(v, *self)
    }
}