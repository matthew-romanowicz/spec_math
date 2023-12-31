//! This crate contains approximations for a set of mathematical functions 
//! commonly referred to as "special functions". The goal of this crate is
//! to eventually contain a full Rust re-implementation of the CEPHES 
//! library with minimal changes in modules such as the [`cephes64`](crate::cephes64)
//! module. The special functions will also be available as traits for easier
//! use in Rust, potentially with modifications to correct bugs or improve performance
//! from the original CEPHES implementation.
//!
//! Currently, the error functions, gamma functions, beta functions, fresnel integrals, 
//! sine and cosine integrals, elliptic integrals, and bessel functions are implemented.

pub mod utils;
pub mod cephes64;
pub mod misc;

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

pub struct FresnelOutput<T> {
    pub s: T,
    pub c: T
}

pub struct SiCiOutput<T> {
    pub si: T,
    pub ci: T
}

/// Implementations of error functions as a trait
pub trait Erf {
    /// Error function
    fn erf(&self) -> Self;

    /// Complimentary error function
    fn erfc(&self) -> Self;

    /// Inverse error function
    fn erf_inv(&self) -> Self;

    /// Inverse complimentary error function
    fn erfc_inv(&self) -> Self;
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
    fn erf_inv(&self) -> f64 {
        //! Uses [`misc::erfinv`](crate::misc::erfinv)
        crate::misc::erfinv(*self)
    }
    fn erfc_inv(&self) -> f64 {
        //! Uses [`misc::erfcinv`](crate::misc::erfcinv)
        crate::misc::erfcinv(*self)
    }
}

/// Implementations of Fresnel integrals as a trait
pub trait Fresnel: Sized {
    /// Fresnel integrals
    fn fresnel(&self) -> FresnelOutput<Self>;
}

impl Fresnel for f64 {
    fn fresnel(&self) -> FresnelOutput<f64> {
        //! Uses [`cephes64::fresnl`](crate::cephes64::fresnl)
        let res = crate::cephes64::fresnl(*self);
        FresnelOutput::<f64> {s: res.0, c: res.1}
    }
}

/// Implementations of sine and cosine integrals as a trait
pub trait SiCi: Sized {
    /// Fresnel integrals
    fn sici(&self) -> SiCiOutput<Self>;
}

impl SiCi for f64 {
    fn sici(&self) -> SiCiOutput<f64> {
        //! Uses [`cephes64::sici`](crate::cephes64::sici)
        let res = crate::cephes64::sici(*self);
        SiCiOutput::<f64> {si: res.0, ci: res.1}
    }
}

/// Implementations of Dawson integral as a trait
pub trait Dawson {
    /// Dawson integral
    fn dawson(&self) -> Self;
}

impl Dawson for f64 {
    fn dawson(&self) -> f64 {
        //! Uses [`cephes64::dawsn`](crate::cephes64::dawsn)
        crate::cephes64::dawsn(*self)
    }
}

/// Implementations of Normal (Gaussian) distribution as a trait
pub trait NormDist {
    /// Normal probability density function
    fn norm_pdf(&self) -> Self;

    /// Normal cumulative distribution function
    fn norm_cdf(&self) -> Self;

    /// Inverse of normal cumulative distribution function
    fn norm_cdf_inv(&self) -> Self;
}

impl NormDist for f64 {
    fn norm_pdf(&self) -> f64 {
        //! Uses [`misc::norm_pdf`](crate::misc::norm_pdf)
        crate::misc::norm_pdf(*self)
    }
    fn norm_cdf(&self) -> f64 {
        //! Uses [`cephes64::ndtr`](crate::cephes64::ndtr)
        crate::cephes64::ndtr(*self)
    }
    fn norm_cdf_inv(&self) -> f64 {
        //! Uses [`cephes64::ndtri`](crate::cephes64::ndtri)
        crate::cephes64::ndtri(*self)
    }
}

/// Implementations of gamma functions as a trait
pub trait Gamma {
    /// Gamma function
    fn gamma(&self) -> Self;

    /// Natural logarithm of the absolute value of the gamma function
    fn lgamma(&self) -> Self;

    /// Regularized lower incomplete gamma function
    fn igamma(&self, x: Self) -> Self;

    /// Regularized upper incomplete gamma function
    fn igammac(&self, x: Self) -> Self;

    /// Inverse of the regularized lower incomplete gamma function
    fn igamma_inv(&self, x: Self) -> Self;

    /// Inverse of the regularized upper incomplete gamma function
    fn igammac_inv(&self, x: Self) -> Self;

    /// Digamma function
    fn digamma(&self) -> Self;

    /// Reciprocal of the gamma function
    fn rgamma(&self) -> Self;
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
    fn igamma_inv(&self, x: f64) -> f64 {
        //! Uses [`cephes64::igami`](crate::cephes64::igami) 
        crate::cephes64::igami(*self, x)
    }
    fn igammac_inv(&self, x: f64) -> f64 {
        //! Uses [`cephes64::igamci`](crate::cephes64::igamci) 
        crate::cephes64::igamci(*self, x)
    }
    fn digamma(&self) -> f64 {
        //! Uses [`cephes64::psi`](crate::cephes64::psi)
        crate::cephes64::psi(*self)
    }
    fn rgamma(&self) -> f64 {
        //! Uses [`cephes64::rgamma`](crate::cephes64::rgamma)
        crate::cephes64::rgamma(*self)
    }
}

/// Implementations of beta functions as a trait
pub trait Beta {
    /// Beta function
    fn beta(&self, other: Self) -> Self;

    /// Natural logarithm of the absolute value of the beta function
    fn lbeta(&self, other: Self) -> Self;

    /// Regularized incomplete beta function
    fn ibeta(&self, a: Self, b: Self) -> Self;

    /// Inverse of the regularized incomplete beta function
    ///
    /// Finds `x` such that self is equal to `x.ibeta(a, b)`
    fn ibeta_inv(&self, a: Self, b: Self) -> Self;
}

impl Beta for f64 {
    fn beta(&self, other: f64) -> f64 {
        //! Uses [`cephes64::beta`](crate::cephes64::beta)
        crate::cephes64::beta(*self, other)
    }
    fn lbeta(&self, other: f64) -> f64 {
        //! Uses [`cephes64::lbeta`](crate::cephes64::lbeta) 
        crate::cephes64::lbeta(*self, other)
    }

    fn ibeta(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::incbet`](crate::cephes64::incbet) 
        crate::cephes64::incbet(a, b, *self)
    }

    fn ibeta_inv(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::incbi`](crate::cephes64::incbi) 
        crate::cephes64::incbi(a, b, *self)
    }
}

/// Implementations of zeta functions as a trait
pub trait Zeta {
    /// Riemann zeta function
    fn zeta(&self) -> Self;

    /// Riemann zeta function minus 1
    fn zetac(&self) -> Self;

    /// Hurwitz zeta function
    fn hzeta(&self, q: Self) -> Self;
}

impl Zeta for f64 {
    fn zeta(&self) -> f64 {
        //! Uses [`cephes64::reimann_zeta`](crate::cephes64::reimann_zeta)
        crate::cephes64::riemann_zeta(*self)
    }
    fn zetac(&self) -> f64 {
        //! Uses [`cephes64::zetac`](crate::cephes64::zetac)
        crate::cephes64::zetac(*self)
    }
    fn hzeta(&self, q: f64) -> f64 {
        //! Uses [`cephes64::zeta`](crate::cephes64::zeta)
        crate::cephes64::zeta(*self, q)
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

    /// Weighted integral of the Bessel function of the first kind
    fn bessel_poly(&self, lambda: Self, nu: Self) -> Self;
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
    fn bessel_poly(&self, lambda: f64, nu: f64) -> f64 {
        //! Uses [`misc::besselpoly`](crate::misc::besselpoly)
        crate::misc::besselpoly(*self, lambda, nu)
    }
}