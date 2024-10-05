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

pub struct ShiChiOutput<T> {
    pub shi: T,
    pub chi: T
}

/// Implementations of error functions as a trait
pub trait Erf {

    /// Error function
    ///
    /// # Description
    ///
    /// Approximation of the error function
    ///
    /// # Domain
    ///
    /// Function is defined for all real numbers.
    ///
    /// # Examples
    ///
    /// Particular Values
    ///
    /// ```
    /// use spec_math::Erf;
    /// 
    /// assert_eq!(0.0_f64.erf(), 0.0);
    /// assert_eq!(f64::INFINITY.erf(), 1.0);
    /// assert_eq!((-f64::INFINITY).erf(), -1.0);
    /// ```
    ///
    /// Symmetry
    ///
    /// ```
    /// use spec_math::Erf;
    /// 
    /// // 0 to 5 in increments of 0.05
    /// for x in (0..=100).map(|i| i as f64 * 0.05) {
    ///     assert_eq!((-x).erf(), -x.erf());
    /// }
    /// ```
    ///
    /// Relationship with incomplete gamma function ([DLMF 7.11.1](https://dlmf.nist.gov/7.11))
    ///
    /// ```
    /// use spec_math::{Erf, Gamma};
    /// 
    /// let sqrt_pi = std::f64::consts::PI.sqrt();
    ///
    /// // 0 to 5 in increments of 0.05
    /// for x in (0..=100).map(|i| i as f64 * 0.05) {
    ///     let a = sqrt_pi * x.erf();
    ///     let b = (0.5_f64).igamma(x * x) * (0.5_f64).gamma();
    ///     if a != 0.0 {
    ///         let rel_err = (a - b) / a;
    ///         assert!(rel_err.abs() < 5e-15); 
    ///     } else {
    ///         let abs_err = a - b;
    ///         assert!(abs_err.abs() < 5e-15); 
    ///     }
    /// }
    /// ```
    fn erf(&self) -> Self;

    /// Compliment of the error function
    ///
    /// # Description
    ///
    /// Approximation of the compliment of the error function. Useful for 
    /// precise evaluation at large inputs.
    ///
    /// # Domain
    ///
    /// Function is defined for all real numbers.
    ///
    /// # Examples
    ///
    /// Relationship with erf
    ///
    /// ```
    /// use spec_math::Erf;
    /// 
    /// // -5 to 5 in increments of 0.05
    /// for x in (-100..=100).map(|i| i as f64 * 0.05) {
    ///     let a = x.erf();
    ///     let b = x.erfc();
    ///     let abs_err = (a + b - 1.0).abs();
    ///     assert!(abs_err < 1e-15);
    /// }
    /// ```
    ///
    /// Relationship with incomplete gamma function ([DLMF 7.11.2](https://dlmf.nist.gov/7.11))
    ///
    /// ```
    /// use spec_math::{Erf, Gamma};
    /// 
    /// let sqrt_pi = std::f64::consts::PI.sqrt();
    ///
    /// // 0 to 5 in increments of 0.05
    /// for x in (0..=100).map(|i| i as f64 * 0.05) {
    ///     let a = sqrt_pi * x.erfc();
    ///     let b = (0.5_f64).igammac(x * x) * (0.5_f64).gamma();
    ///     if a != 0.0 {
    ///         let rel_err = (a - b) / a;
    ///         assert!(rel_err.abs() < 3e-14); 
    ///     } else {
    ///         let abs_err = a - b;
    ///         assert!(abs_err.abs() < 5e-15); 
    ///     }
    /// }
    /// ```
    fn erfc(&self) -> Self;

    /// Inverse error function
    ///
    /// # Description
    ///
    /// Approximation of the inverse of the error function.
    ///
    /// # Domain
    ///
    /// Returns `NAN` for inputs less than -1.0 and greater than 1.0.
    /// Returns infinity for 1.0 and negative infinity for -1.0. Note
    /// that the function becomes very steep near -1.0 and 1.0 and so
    /// the floating point precision can heavily influence the output
    /// in those regions.
    ///
    /// # Examples
    ///
    /// Domain
    ///
    /// ```
    /// use spec_math::Erf;
    ///
    /// assert!((-1.1).erf_inv().is_nan());
    /// assert!((1.1).erf_inv().is_nan());
    ///
    /// assert_eq!((-1.0).erf_inv(), -f64::INFINITY);
    /// assert_eq!((1.0).erf_inv(), f64::INFINITY);
    /// ```
    ///
    /// Relationship with erf
    ///
    /// ```
    /// use spec_math::Erf;
    /// 
    /// // -2 to 2 in increments of 0.02
    /// for x in (-100..=100).map(|i| i as f64 * 0.02) {
    ///     let y = x.erf();
    ///     let x2 = y.erf_inv();
    ///     if x != 0.0 {
    ///         let rel_err = (x - x2) / x;
    ///         assert!(rel_err.abs() < 1e-14); 
    ///     } else {
    ///         let abs_err = x - x2;
    ///         assert!(abs_err.abs() < 1e-15); 
    ///     }
    /// }
    /// ```
    fn erf_inv(&self) -> Self;

    /// Inverse complimentary error function
    ///
    /// # Description
    ///
    /// Approximation of the inverse of the complimentary error function.
    ///
    /// # Domain
    ///
    /// Returns `NAN` for inputs less than 0.0 and greater than 2.0.
    /// Returns infinity for 0.0 and negative infinity for 2.0. Note
    /// that the function becomes very steep near 2.0 and 0.0 and so
    /// the floating point precision can heavily influence the output
    /// in those regions. This is less significant in the 0.0 region than
    /// the 2.0 region.
    ///
    /// # Examples
    ///
    /// Domain
    ///
    /// ```
    /// use spec_math::Erf;
    ///
    /// assert!((-0.1).erfc_inv().is_nan());
    /// assert!((2.1).erfc_inv().is_nan());
    ///
    /// assert_eq!((0.0).erfc_inv(), f64::INFINITY);
    /// assert_eq!((2.0).erfc_inv(), -f64::INFINITY);
    /// ```
    ///
    /// Relationship with erfc
    ///
    /// ```
    /// use spec_math::Erf;
    /// 
    /// // -2 to 2 in increments of 0.02
    /// for x in (-100..=100).map(|i| i as f64 * 0.02) {
    ///     let y = x.erfc();
    ///     let x2 = y.erfc_inv();
    ///     if x != 0.0 {
    ///         let rel_err = (x - x2) / x;
    ///         assert!(rel_err.abs() < 1e-14); 
    ///     } else {
    ///         let abs_err = x - x2;
    ///         assert!(abs_err.abs() < 1e-15); 
    ///     }
    /// }
    /// ```
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
    /// Sine and cosine integrals
    fn sici(&self) -> SiCiOutput<Self>;
}

impl SiCi for f64 {
    fn sici(&self) -> SiCiOutput<f64> {
        //! Uses [`cephes64::sici`](crate::cephes64::sici)
        let res = crate::cephes64::sici(*self);
        SiCiOutput::<f64> {si: res.0, ci: res.1}
    }
}

/// Implementations of hyperbolic sine and cosine integrals as a trait
pub trait ShiChi: Sized {
    /// Hyperbolic sine and cosine integrals
    fn shichi(&self) -> ShiChiOutput<Self>;
}

impl ShiChi for f64 {
    fn shichi(&self) -> ShiChiOutput<f64> {
        //! Uses [`cephes64::shichi`](crate::cephes64::shichi)
        let res = crate::cephes64::shichi(*self);
        ShiChiOutput::<f64> {shi: res.0, chi: res.1}
    }
}

/// Implementations of the generalized exponential integral as a trait
pub trait Expn {
    /// Generalized exponential integral
    fn expn(&self, n: i32) -> Self;
}

impl Expn for f64 {
    fn expn(&self, n: i32) -> f64 {
        //! Uses [`cephes64::expn`](crate::cephes64::expn)
        crate::cephes64::expn(n, *self)
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

/// Implementations of polylogarithms as a trait
pub trait Polylog {
    /// Dilogarithm integral
    fn li2(&self) -> Self;
}

impl Polylog for f64 {
    fn li2(&self) -> f64 {
        //! Uses [`cephes64::spence`](crate::cephes64::spence)
        crate::cephes64::spence(*self)
    }
}

/// Implementations of Normal (Gaussian) distribution as a trait
pub trait NormDist {
    /// Normal probability density function
    ///
    /// # Description:
    ///
    /// Returns the probability density of the normal distribution 
    /// function at `self`.
    ///
    /// # Domain:
    ///
    /// Function is defined for all real numbers. Returns `NAN` if
    /// input is `NAN`.
    ///
    /// # Examples:
    ///
    /// Particular values:
    ///
    /// ```
    /// use spec_math::NormDist;
    /// 
    /// // Normal PDF approaches 0 at infinity and negative infinity
    /// assert_eq!((f64::INFINITY).norm_pdf(), 0.0);
    /// assert_eq!((-f64::INFINITY).norm_pdf(), 0.0);
    /// 
    /// // Normal PDF is equal to 1/sqrt(2 * pi) at 0
    /// assert_eq!((0.0_f64).norm_pdf(), 1.0 / (2.0 * std::f64::consts::PI).sqrt());
    /// ```
    ///
    /// Symmetry about 0.0:
    ///
    /// ```
    /// use spec_math::NormDist;
    /// 
    /// // 0.0 to 10.0 at 0.1 increments
    /// for x in (0..=100).map(|i| i as f64 * 0.1) {
    ///     assert_eq!(x.norm_pdf(), (-x).norm_pdf());
    /// }    
    /// ```
    fn norm_pdf(&self) -> Self;

    /// Normal cumulative distribution function
    ///
    /// # Description:
    ///
    /// Returns the area under the normal distribution curve from
    /// negative infinity to `self`.
    ///
    /// # Domain:
    ///
    /// Function is defined for all real numbers. Returns `NAN` if
    /// input is `NAN`.
    ///
    /// # Examples:
    ///
    /// Particular values:
    ///
    /// ```
    /// use spec_math::NormDist;
    /// 
    /// // Normal CDF approaches 1 at infinity and 0 at negative infinity
    /// assert_eq!((f64::INFINITY).norm_cdf(), 1.0);
    /// assert_eq!((-f64::INFINITY).norm_cdf(), 0.0);
    /// 
    /// // One-tailed z-value for 95% confidence is approximately 1.645
    /// let p_95 = (1.645_f64).norm_cdf();
    /// assert!((p_95 - 0.95).abs() < 1e-3);
    ///
    /// // One-tailed z-value for 99% confidence is approximately 2.33
    /// let p_99 = (2.33_f64).norm_cdf();
    /// assert!((p_99 - 0.99).abs() < 1e-3);
    /// ```
    ///
    /// Association with the error function ([DLMF 7.20(iii)](https://dlmf.nist.gov/7.20#iii))
    ///
    /// ```
    /// use spec_math::{Erf, NormDist};
    /// 
    /// // 0.0 to 10.0 at 0.1 increments
    /// for x in (0..=100).map(|i| i as f64 * 0.1) {
    ///     let a = x.norm_cdf();
    ///     let b = 0.5 * (-x / (2.0_f64).sqrt()).erfc();
    ///     let abs_err = a - b;
    ///     assert!(abs_err.abs() < 1e-14);
    /// }    
    /// ```
    fn norm_cdf(&self) -> Self;

    /// Inverse of normal cumulative distribution function
    ///
    /// # Description:
    ///
    /// Returns `x` such that `x.norm_cdf()` is equal to `self`.
    ///
    /// # Domain:
    ///
    /// Returns `NAN` for input values outside of [0.0, 1.0]. Returns negative
    /// infinity for 0.0 and positive infinity for 1.0.
    ///
    /// # Examples:
    ///
    /// Particular values:
    ///
    /// ```
    /// use spec_math::NormDist;
    /// 
    /// // Normal CDF approaches 1 at infinity and 0 at negative infinity
    /// assert_eq!((0.0_f64).norm_cdf_inv(), -f64::INFINITY);
    /// assert_eq!((1.0_f64).norm_cdf_inv(), f64::INFINITY);
    /// 
    /// // One-tailed z-value for 95% confidence is approximately 1.645
    /// let z_95 = (0.95_f64).norm_cdf_inv();
    /// assert!((z_95 - 1.645).abs() < 1e-3);
    ///
    /// // One-tailed z-value for 99% confidence is approximately 2.33
    /// let z_99 = (0.99_f64).norm_cdf_inv();
    /// assert!((z_99 - 2.33).abs() < 5e-3);
    /// ```
    ///
    /// Association with normal CDF:
    ///
    /// ```
    /// use spec_math::NormDist;
    /// 
    /// // -5.0 to 5.0 at 0.1 increments
    /// for x in (-50..=50).map(|i| i as f64 * 0.1) {
    ///     let p = x.norm_cdf();
    ///     let x2 = p.norm_cdf_inv();
    ///     let abs_err = x - x2;
    ///     assert!(abs_err.abs() < 1e-10);
    /// }    
    /// ```
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

/// Implementations of the student's t distribution as a trait
pub trait TDist {
    /// Student's t probability density function
    fn t_pdf(&self, df: Self) -> Self;

    /// Student's t cumulative distribution function
    fn t_cdf(&self, df: isize) -> Self;

    /// Inverse of student's t cumulative distribution function
    fn t_cdf_inv(&self, df: isize) -> Self;
}

impl TDist for f64 {
    fn t_pdf(&self, df: f64) -> f64 {
        //! Uses [`misc::t_pdf`](crate::misc::t_pdf)
        crate::misc::t_pdf(df, *self)
    }
    fn t_cdf(&self, df: isize) -> f64 {
        //! Uses [`cephes64::stdtr`](crate::cephes64::stdtr)
        crate::cephes64::stdtr(df, *self)
    }
    fn t_cdf_inv(&self, df: isize) -> f64 {
        //! Uses [`cephes64::stdtri`](crate::cephes64::stdtri)
        crate::cephes64::stdtri(df, *self)
    }
}

/// Implementations of binomial distribution as a trait
pub trait BinomDist {
    /// Binomial probability mass function
    fn binom_pmf(&self, k: f64, n: isize) -> Self;

    /// Binomial cumulative distribution function
    fn binom_cdf(&self, k: Self, n: isize) -> Self;

    /// Compliment of the binomial cumulative distribution function
    fn binom_cdfc(&self, k: Self, n: isize) -> Self;

    /// Inverse of the binomial cumulative distribution function
    fn binom_cdf_inv(&self, k: Self, n: isize) -> Self;
}

impl BinomDist for f64 {
    fn binom_pmf(&self, k: f64, n: isize) -> f64 {
        //! Uses [`misc::binom_pmf`](crate::misc::binom_pmf)
        crate::misc::binom_pmf(k, n as i32, *self)
    }
    fn binom_cdf(&self, k: f64, n: isize) -> f64 {
        //! Uses [`cephes64::bdtr`](crate::cephes64::bdtr)
        crate::cephes64::bdtr(k, n, *self)
    }
    fn binom_cdfc(&self, k: f64, n: isize) -> f64 {
        //! Uses [`cephes64::bdtrc`](crate::cephes64::bdtrc)
        crate::cephes64::bdtrc(k, n, *self)
    }
    fn binom_cdf_inv(&self, k: Self, n: isize) -> f64 {
        //! Uses [`cephes64::bdtri`](crate::cephes64::bdtri)
        crate::cephes64::bdtri(k, n, *self)
    }
}

/// Implementations of negative binomial distribution as a trait
pub trait NBinomDist {
    /// Negative binomial probability mass function
    fn nbinom_pmf(&self, k: isize, n: isize) -> Self;

    /// Negative binomial cumulative distribution function
    fn nbinom_cdf(&self, k: isize, n: isize) -> Self;

    /// Compliment of the negative binomial cumulative distribution function
    fn nbinom_cdfc(&self, k: isize, n: isize) -> Self;

    /// Inverse of the negative binomial cumulative distribution function
    fn nbinom_cdf_inv(&self, k: isize, n: isize) -> Self;
}

impl NBinomDist for f64 {
    fn nbinom_pmf(&self, k: isize, n: isize) -> f64 {
        //! Uses [`misc::nbinom_pmf`](crate::misc::nbinom_pmf)
        crate::misc::nbinom_pmf(k as i32, n as i32, *self)
    }
    fn nbinom_cdf(&self, k: isize, n: isize) -> f64 {
        //! Uses [`cephes64::nbdtr`](crate::cephes64::nbdtr)
        crate::cephes64::nbdtr(k, n, *self)
    }
    fn nbinom_cdfc(&self, k: isize, n: isize) -> f64 {
        //! Uses [`cephes64::nbdtrc`](crate::cephes64::nbdtrc)
        crate::cephes64::nbdtrc(k, n, *self)
    }
    fn nbinom_cdf_inv(&self, k: isize, n: isize) -> f64 {
        //! Uses [`cephes64::nbdtri`](crate::cephes64::nbdtri)
        crate::cephes64::nbdtri(k, n, *self)
    }
}

/// Implementations of gamma functions as a trait
pub trait Gamma {

    /// Gamma function
    ///
    /// # Description
    ///
    /// Approximation of the gamma function
    ///
    /// # Domain
    ///
    /// Returns NAN for the following inputs: NAN, 0.0, negative 
    /// integers, and negative infinity. Returns positivy infinity 
    /// when a positive overflow condition is encountered.
    ///
    /// # Examples
    ///
    /// Gamma of a positive integer `n` is equal to the factorial of `n - 1`:
    ///
    /// ```
    /// use spec_math::Gamma;
    /// 
    /// let mut fac = 1.0;
    /// for n in (1..20).map(|i| i as f64) {
    ///     assert_eq!(fac, n.gamma());
    ///     fac *= n;
    /// }
    /// ```
    ///
    /// Gamma of a nonpositive integer is NAN:
    ///
    /// ```
    /// use spec_math::Gamma;
    /// 
    /// for n in (-100..=0).map(|i| i as f64) {
    ///     assert!(n.gamma().is_nan());
    /// }
    /// ```
    ///
    /// Gamma of a real number `x + 1` is equal to `x * x.gamma()`
    ///
    /// ```
    /// use spec_math::Gamma;
    ///    
    /// // Positive numbers at increments of 0.5 (0.5, 1.0, 1.5, ...)
    /// for x in (1..100).map(|i| 0.5 * (i as f64)) {
    ///     let a = (x + 1.0).gamma();
    ///     let b = x * x.gamma();
    ///     let rel_err = (a - b) / b;
    ///     assert!(rel_err < 1e-15);
    /// }
    ///
    /// // Negative numbers at increments of 0.5 offset by 0.25
    /// // (-0.25, -0.75, -1.25, ...)
    /// for x in (0..100).map(|i| -0.5 * (i as f64) - 0.25) {
    ///     let a = (x + 1.0).gamma();
    ///     let b = x * x.gamma();
    ///     let rel_err = (a - b) / b;
    ///     assert!(rel_err < 1e-15);        
    /// }
    /// ```
    ///
    /// Particular values (from [wikipedia](https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function)):
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// let a = (1.0 / 2.0).gamma();
    /// let b = std::f64::consts::PI.sqrt();
    /// let rel_err = (a - b) / b;
    /// assert!(rel_err.abs() < 1e-15);
    /// ```
    fn gamma(&self) -> Self;

    /// Natural logarithm of the absolute value of the gamma function
    ///
    /// # Description
    ///
    /// Approximation of the natural logarithm of the absolute value of the 
    /// gamma function. Has much slower growth characteristics than the 
    /// standard gamma function.
    ///
    /// # Domain
    ///
    /// Returns NAN for negative infinity. Returns positivy infinity 
    /// for 0.0, negative integers, and when a positive overflow 
    /// condition is encountered.
    ///
    /// # Examples
    ///
    /// Lgamma is the natural logarithm of gamma:
    ///
    /// ```
    /// use spec_math::Gamma;
    ///    
    /// // Positive numbers at increments of 0.5 (0.5, 1.0, 1.5, ...)
    /// for x in (1..100).map(|i| 0.5 * (i as f64)) {
    ///     let a = x.lgamma();
    ///     let b = x.gamma().ln();
    ///     let abs_err = a - b;
    ///     assert!(abs_err.abs() < 1e-13);
    /// }
    ///
    /// // Negative numbers at increments of 0.5 offset by 0.25
    /// // (-0.25, -0.75, -1.25, ...)
    /// for x in (0..100).map(|i| -0.5 * (i as f64) - 0.25) {
    ///     let a = x.lgamma();
    ///     let b = x.gamma().abs().ln();
    ///     let rel_err = (a - b) / b;
    ///     assert!(rel_err.abs() < 1e-13);        
    /// }
    /// ```
    ///
    /// Lgamma of a nonpositive integer is infinity:
    ///
    /// ```
    /// use spec_math::Gamma;
    /// 
    /// for n in (-100..=0).map(|i| i as f64) {
    ///     assert_eq!(n.lgamma(), f64::INFINITY);
    /// }
    /// ```
    fn lgamma(&self) -> Self;

    /// Regularized lower incomplete gamma function
    ///
    /// # Examples
    ///
    /// Upper and lower incomplete gamma functions sum to 1:
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// for s in (1..100).map(|i| i as f64) {
    ///     for x in (0..100).map(|i| i as f64) {
    ///         let lower = s.igamma(x);
    ///         let upper = s.igammac(x);
    ///         assert!(lower + upper - 1.0 < 1e-14);
    ///     }
    /// }
    /// ```
    ///
    /// Lower gamma recurrence relation (from [DLMF 8.8.1](https://dlmf.nist.gov/8.8))
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// for s in (1..100).map(|i| i as f64) {
    ///     for x in (1..100).map(|i| i as f64) {
    ///         // Results must be multiplied by s.gamma() to de-regularize them
    ///         let a = (s + 1.0).igamma(x) * (s + 1.0).gamma();
    ///         let b = s * s.igamma(x) * s.gamma() - x.powf(s) * (-x).exp();
    ///         if a != 0.0 {
    ///             let rel_err = (a - b) / a;
    ///             assert!(rel_err.abs() < 8e-12); 
    ///         } else {
    ///             let abs_err = a - b;
    ///             assert!(abs_err.abs() < 1e-14); 
    ///         }
    ///     }
    /// }
    /// ```
    fn igamma(&self, x: Self) -> Self;

    /// Regularized upper incomplete gamma function
    ///
    /// # Examples
    ///
    /// Upper gamma recurrence relation (from [DLMF 8.8.2](https://dlmf.nist.gov/8.8))
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// for s in (1..100).map(|i| i as f64) {
    ///     for x in (1..100).map(|i| i as f64) {
    ///         // Results must be multiplied by s.gamma() to de-regularize them
    ///         let a = (s + 1.0).igammac(x) * (s + 1.0).gamma();
    ///         let b = s * s.igammac(x) * s.gamma() + x.powf(s) * (-x).exp();
    ///         let rel_err = (a - b) / b;
    ///         if a != 0.0 {
    ///             let rel_err = (a - b) / a;
    ///             assert!(rel_err.abs() < 1e-13); 
    ///         } else {
    ///             let abs_err = a - b;
    ///             assert!(abs_err.abs() < 1e-14); 
    ///         }
    ///     }
    /// }
    /// ```
    ///
    /// See also [`igamma`](crate::Gamma::igamma).
    fn igammac(&self, x: Self) -> Self;

    /// Inverse of the regularized lower incomplete gamma function
    ///
    /// # Examples
    ///
    /// `igamma_inv` is the inverse of `igamma`:
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// for s in (1..100).map(|i| i as f64) {
    ///     for x in (0..20).map(|i| i as f64) {
    ///         let y = s.igamma(x);
    ///         let x2 = s.igamma_inv(y);
    ///         if x != 0.0 {
    ///             let rel_err = (x - x2) / x;
    ///             assert!(rel_err.abs() < 1e-9); 
    ///         } else {
    ///             let abs_err = x - x2;
    ///             assert!(abs_err.abs() < 1e-14); 
    ///         }
    ///     }
    /// }
    /// ```    
    fn igamma_inv(&self, x: Self) -> Self;

    /// Inverse of the regularized upper incomplete gamma function
    fn igammac_inv(&self, x: Self) -> Self;

    /// Digamma function
    ///
    /// Approximates the derivative of the gamma function.
    ///
    /// # Examples
    ///
    /// Digamma of a real number `x + 1` is equal to `x.digamma() + 1.0 / x`
    ///
    /// ```
    /// use spec_math::Gamma;
    ///    
    /// // Positive numbers at increments of 0.5 (0.5, 1.0, 1.5, ...)
    /// for x in (1..100).map(|i| 0.5 * (i as f64)) {
    ///     let a = (x + 1.0).digamma();
    ///     let b = x.digamma() + 1.0 / x;
    ///     let rel_err = (a - b) / b;
    ///     assert!(rel_err < 1e-14);
    /// }
    ///
    /// // Negative numbers at increments of 0.5 offset by 0.25
    /// // (-0.25, -0.75, -1.25, ...)
    /// for x in (0..100).map(|i| -0.5 * (i as f64) - 0.25) {
    ///     let a = (x + 1.0).digamma();
    ///     let b = x.digamma() + 1.0 / x;
    ///     let rel_err = (a - b) / b;
    ///     assert!(rel_err < 7e-14);        
    /// }
    /// ```
    ///
    /// Digamma is the derivative of lgamma:
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// const DX: f64 = 1e-10;
    ///    
    /// // Positive numbers at increments of 0.5 (0.5, 1.0, 1.5, ...)
    /// for x in (1..100).map(|i| 0.5 * (i as f64)) {
    ///     let digam = x.digamma();
    ///     let dydx = ((x + 0.5 * DX).lgamma() - (x - 0.5 * DX).lgamma()) / DX;
    ///     let rel_err = (digam - dydx) / dydx;
    ///     assert!(rel_err.abs() < 1e-3); // dydx approximation is imprecise
    /// }
    ///
    /// // Negative numbers at increments of 0.5 offset by 0.25
    /// // (-0.25, -0.75, -1.25, ...)
    /// for x in (0..100).map(|i| -0.5 * (i as f64) - 0.25) {
    ///     let digam = x.digamma();
    ///     let dydx = ((x + 0.5 * DX).lgamma() - (x - 0.5 * DX).lgamma()) / DX;
    ///     let rel_err = (digam - dydx) / dydx;
    ///     assert!(rel_err.abs() < 1e-2); // dydx approximation is imprecise      
    /// }
    /// ```
    fn digamma(&self) -> Self;

    /// Reciprocal of the gamma function
    ///
    /// # Examples
    ///
    /// Gamma and reciprocal gamma functions multiply to 1:
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// for x in (1..100).map(|i| i as f64) {
    ///     let gam = x.gamma();
    ///     let rgam = x.rgamma();
    ///     assert!(gam * rgam - 1.0 < 1e-13);
    /// }
    /// ```
    fn rgamma(&self) -> Self;

    /// Pochhammer Symbol
    ///
    /// # Description
    ///
    /// Approximates the value of the [Pochhammer symbol](https://mathworld.wolfram.com/PochhammerSymbol.html)
    ///
    /// # Examples
    ///
    /// Association with gamma:
    ///
    /// ```
    /// use spec_math::Gamma;
    ///
    /// // 0.5 to 15 in increments of 0.5
    /// for x in (1..=30).map(|i| i as f64 * 0.5) {
    ///     // 0 to 5 in increments of 1.0
    ///     for n in (0..5).map(|i| i as f64) {
    ///         let a = x.poch(n);
    ///         let b = (x + n).gamma() / x.gamma();
    ///         let err = (a - b) / a;
    ///         assert!(err < 1e-15);
    ///     }
    /// }
    /// ```
    fn poch(&self, n: Self) -> Self;
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
    fn poch(&self, n: f64) -> f64 {
        crate::misc::poch(*self, n)
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

/// Implementations of beta distribution as a trait
pub trait BetaDist {
    /// Beta probability density function
    fn beta_pdf(&self, a: Self, b: Self) -> Self;

    /// Beta cumulative distribution function
    fn beta_cdf(&self, a: Self, b: Self) -> Self;

    /// Inverse of the beta cumulative distribution function
    fn beta_cdf_inv(&self, a: Self, b: Self) -> Self;
}

impl BetaDist for f64 {
    fn beta_pdf(&self, a: f64, b: f64) -> f64 {
        //! Uses [`misc::beta_pdf`](crate::misc::beta_pdf)
        crate::misc::beta_pdf(a, b, *self)
    }
    fn beta_cdf(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::btdtr`](crate::cephes64::btdtr)
        crate::cephes64::btdtr(a,b, *self)
    }
    fn beta_cdf_inv(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::incbi`](crate::cephes64::incbi) 
        crate::cephes64::incbi(a, b, *self)
    }
}

/// Implementations of Chi squared distribution as a trait
pub trait Chi2Dist {
    /// Chi squared probability density function
    fn chi2_pdf(&self, df: Self) -> Self;

    /// Chi squared cumulative distribution function
    fn chi2_cdf(&self, df: Self) -> Self;

    /// Compliment of the chi squared cumulative distribution function
    fn chi2_cdfc(&self, df: Self) -> Self;

    /// Inverse of the chi squared cumulative distribution function
    fn chi2_cdf_inv(&self, df: Self) -> Self;
}

impl Chi2Dist for f64 {
    fn chi2_pdf(&self, df: f64) -> f64 {
        //! Uses [`misc::chi2_pdf`](crate::misc::chi2_pdf)
        crate::misc::chi2_pdf(df, *self)
    }
    fn chi2_cdf(&self, df: f64) -> f64 {
        //! Uses [`cephes64::chdtr`](crate::cephes64::chdtr)
        crate::cephes64::chdtr(df, *self)
    }
    fn chi2_cdfc(&self, df: f64) -> f64 {
        //! Uses [`cephes64::chdtrc`](crate::cephes64::chdtrc)
        crate::cephes64::chdtrc(df, *self)
    }
    fn chi2_cdf_inv(&self, df: f64) -> f64 {
        //! Uses [`cephes64::chdtri`](crate::cephes64::chdtri) 
        crate::cephes64::chdtri(df, *self)
    }
}

/// Implementations of the Poisson distribution as a trait
pub trait PoissonDist {
    /// Poisson probability mass function
    fn pois_pmf(&self, k: i32) -> Self;

    /// Poisson cumulative distribution function
    fn pois_cdf(&self, k: Self) -> Self;

    /// Compliment of the Poisson cumulative distribution function
    fn pois_cdfc(&self, k: Self) -> Self;

    /// Inverse of the Poisson cumulative distribution function
    fn pois_cdf_inv(&self, k: isize) -> Self;
}

impl PoissonDist for f64 {
    fn pois_pmf(&self, k: i32) -> f64 {
        //! Uses [`misc::pois_pmf`](crate::misc::pois_pmf)
        crate::misc::pois_pmf(k, *self)
    }
    fn pois_cdf(&self, k: f64) -> f64 {
        //! Uses [`cephes64::pdtr`](crate::cephes64::pdtr)
        crate::cephes64::pdtr(k, *self)
    }
    fn pois_cdfc(&self, k: f64) -> f64 {
        //! Uses [`cephes64::pdtrc`](crate::cephes64::pdtrc)
        crate::cephes64::pdtrc(k, *self)
    }
    fn pois_cdf_inv(&self, k: isize) -> f64 {
        //! Uses [`cephes64::pdtri`](crate::cephes64::pdtri) 
        crate::cephes64::pdtri(k, *self)
    }
}

/// Implementations of f distribution as a trait
pub trait FDist {
    /// F probability density function
    fn f_pdf(&self, a: Self, b: Self) -> Self;

    /// F cumulative distribution function
    fn f_cdf(&self, a: Self, b: Self) -> Self;

    /// Compliment of the f cumulative distribution function
    fn f_cdfc(&self, a: Self, b: Self) -> Self;

    /// Inverse of the f cumulative distribution function
    fn f_cdf_inv(&self, a: Self, b: Self) -> Self;
}

impl FDist for f64 {
    fn f_pdf(&self, a: f64, b: f64) -> f64 {
        //! Uses [`misc::f_pdf`](crate::misc::f_pdf)
        crate::misc::f_pdf(a, b, *self)
    }
    fn f_cdf(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::fdtr`](crate::cephes64::fdtr)
        crate::cephes64::fdtr(a,b, *self)
    }
    fn f_cdfc(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::fdtrc`](crate::cephes64::fdtrc)
        crate::cephes64::fdtrc(a,b, *self)
    }
    fn f_cdf_inv(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::fdtri`](crate::cephes64::fdtri) 
        crate::cephes64::fdtri(a, b, *self)
    }
}

/// Implementations of gamma distribution as a trait
pub trait GammaDist {
    /// Gamma probability density function
    fn gamma_pdf(&self, a: Self, b: Self) -> Self;

    /// Gamma cumulative distribution function
    fn gamma_cdf(&self, a: Self, b: Self) -> Self;

    /// Compliment of the gamma cumulative distribution function
    fn gamma_cdfc(&self, a: Self, b: Self) -> Self;

    /// Inverse of the gamma cumulative distribution function
    fn gamma_cdf_inv(&self, a: Self, b: Self) -> Self;
}

impl GammaDist for f64 {
    fn gamma_pdf(&self, a: f64, b: f64) -> f64 {
        //! Uses [`misc::gamma_pdf`](crate::misc::gamma_pdf)
        crate::misc::gamma_pdf(a, b, *self)
    }
    fn gamma_cdf(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::gdtr`](crate::cephes64::gdtr)
        crate::cephes64::gdtr(a,b, *self)
    }
    fn gamma_cdfc(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::gdtrc`](crate::cephes64::gdtrc)
        crate::cephes64::gdtrc(a,b, *self)
    }
    fn gamma_cdf_inv(&self, a: f64, b: f64) -> f64 {
        //! Uses [`cephes64::gdtri`](crate::cephes64::gdtri) 
        crate::cephes64::gdtri(a, b, *self)
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
    ///
    /// # Description:
    ///
    /// Approximates the complete elliptic integral of the second kind 
    /// with `self` as the k squared term.
    ///
    /// [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_integral#Complete_elliptic_integral_of_the_second_kind)
    ///
    /// # Domain:
    ///
    /// Returns `NAN` for values greater than 1.0
    ///
    /// # Examples:
    ///
    /// Particular values:
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // E(0) = pi / 2
    /// assert_eq!((0.0_f64).ellip_e(), 0.5 * std::f64::consts::PI);
    ///
    /// // E(1) = 1
    /// assert_eq!((1.0_f64).ellip_e(), 1.0);
    /// ```
    ///
    /// Calculation of the polar circumference of Earth from WGS84
    /// constants:
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // Equatorial radius of Earth (km)
    /// let a: f64 = 6378.1370;
    ///
    /// // Flattening of Earth
    /// let f: f64 = 1.0 / 298.257223563;
    ///
    /// // Eccentricity of Earth squared
    /// let ecc2 = f * (2.0 - f);
    ///
    /// // Polar circumference of Earth (km)
    /// let pc: f64 = 40007.863;
    ///
    /// let pc_calc = 4.0 * a * ecc2.ellip_e();
    ///
    /// assert!((pc - pc_calc).abs() < 1e-3);
    /// ```
    ///
    /// Association with the incomplete elliptic integral of the 
    /// second kind ([DLMF 19.2.8](https://dlmf.nist.gov/19.2#E8)):
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.0 to 1.0 in 0.01 increments
    /// for x in (0..=100).map(|i| i as f64 * 0.01) {
    ///     let p_e = x.ellip_e();
    ///     let i_e = x.ellip_e_inc(std::f64::consts::FRAC_PI_2);
    ///     let rel_err = (i_e - p_e) / p_e;
    ///     assert!(rel_err.abs() < 1e-15);
    /// }
    ///
    /// // -100.0 to -0.5 in 0.5 increments
    /// for x in (-200..0).map(|i| i as f64 * 0.5) {
    ///     let p_e = x.ellip_e();
    ///     let i_e = x.ellip_e_inc(std::f64::consts::FRAC_PI_2);
    ///     let rel_err = (i_e - p_e) / p_e;
    ///     assert!(rel_err.abs() < 1e-15);
    /// }
    /// ```
    /// 
    /// Association with Carlson's symmetric elliptic integral of the 
    /// second kind ([DLMF 19.25.1](https://dlmf.nist.gov/19.25#E1)):
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.0 to 1.0 in 0.01 increments
    /// for x in (0..=100).map(|i| i as f64 * 0.01) {
    ///     let p_e = x.ellip_e();
    ///     let rg = (1.0 - x).carlson_rg(0.0, 1.0) * 2.0;
    ///     let rel_err = (rg - p_e) / p_e;
    ///     assert!(rel_err.abs() < 1e-14);
    /// }
    ///
    /// // -100.0 to -0.5 in 0.5 increments
    /// for x in (-200..0).map(|i| i as f64 * 0.5) {
    ///     let p_e = x.ellip_e();
    ///     let rg = (1.0 - x).carlson_rg(0.0, 1.0) * 2.0;
    ///     let rel_err = (rg - p_e) / p_e;
    ///     assert!(rel_err.abs() < 1e-14);
    /// }
    /// ```  
    fn ellip_e(&self) -> Self;

    /// Complete elliptic integral of the first kind
    ///
    /// # Description:
    ///
    /// Approximates the complete elliptic integral of the first kind 
    /// with `self` as the k squared term.
    ///
    /// [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_integral#Complete_elliptic_integral_of_the_first_kind)
    ///
    /// # Domain:
    ///
    /// Returns `NAN` for values greater than 1.0. Returns infinity for 1.0.
    ///
    /// # Examples:
    ///
    /// Particular values:
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // K(0) = pi / 2
    /// assert_eq!((0.0_f64).ellip_k(), 0.5 * std::f64::consts::PI);
    ///
    /// // K(x) approaches 0 as x approaches negative infinity
    /// assert_eq!((-f64::INFINITY).ellip_k(), 0.0);
    ///
    /// // K(x) approaches infinity as x approaches 1.0
    /// assert_eq!((1.0_f64).ellip_k(), f64::INFINITY);
    /// ```
    ///
    /// Association with the incomplete elliptic integral of the 
    /// first kind ([DLMF 19.2.8](https://dlmf.nist.gov/19.2#E8)):
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.0 to 0.99 in 0.01 increments
    /// for x in (0..100).map(|i| i as f64 * 0.01) {
    ///     let p_k = x.ellip_k();
    ///     let i_k = x.ellip_k_inc(std::f64::consts::FRAC_PI_2);
    ///     let rel_err = (i_k - p_k) / p_k;
    ///     assert!(rel_err.abs() < 1e-15);
    /// }
    ///
    /// // -100.0 to -0.5 in 0.5 increments
    /// for x in (-200..0).map(|i| i as f64 * 0.5) {
    ///     let p_k = x.ellip_k();
    ///     let i_k = x.ellip_k_inc(std::f64::consts::FRAC_PI_2);
    ///     let rel_err = (i_k - p_k) / p_k;
    ///     assert!(rel_err.abs() < 1e-15);
    /// }
    /// ```    
    ///
    /// Association with Carlson's symmetric elliptic integral of the 
    /// first kind ([DLMF 19.25.1](https://dlmf.nist.gov/19.25#E1))
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.0 to 0.99 in 0.01 increments
    /// for x in (0..100).map(|i| i as f64 * 0.01) {
    ///     let p_k = x.ellip_k();
    ///     let rf = (1.0 - x).carlson_rf(0.0, 1.0);
    ///     let rel_err = (rf - p_k) / p_k;
    ///     assert!(rel_err.abs() < 1e-15);
    /// }
    ///
    /// // -100.0 to -0.5 in 0.5 increments
    /// for x in (-200..0).map(|i| i as f64 * 0.5) {
    ///     let p_k = x.ellip_k();
    ///     let rf = (1.0 - x).carlson_rf(0.0, 1.0);
    ///     let rel_err = (rf - p_k) / p_k;
    ///     assert!(rel_err.abs() < 1e-15);
    /// }
    /// ```      
    fn ellip_k(&self) -> Self;

    /// Complete elliptic integral of the third kind
    ///
    /// # Description:
    ///
    /// Approximates the complete elliptic integral of the third kind 
    /// with `self` as the k squared term.
    ///
    /// [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_integral#Complete_elliptic_integral_of_the_third_kind)
    ///
    /// # Domain:
    ///
    /// Function returns `NAN` if `n == 1.0` or `self >= 1.0`
    ///
    /// # Examples:
    fn ellip_pi(&self, n: Self) -> Self;

    /// Incomplete elliptic integral of the second kind
    ///
    /// # Description:
    ///
    /// Approximates the incomplete elliptic integral of the second kind 
    /// with `self` as the k squared term.
    ///
    /// [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_integral#Incomplete_elliptic_integral_of_the_second_kind)
    ///
    /// # Examples:
    ///
    /// Association with Carlson's symmetric elliptic integral of the 
    /// first kind and Carlson's elliptic integral symmetric in only 
    /// two variables ([DLMF 19.25.9](https://dlmf.nist.gov/19.25#E9)):
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.1 to 1.5 in 0.1 increments
    /// for phi in (1..=15).map(|i| i as f64 * 0.1) {
    ///     let c = phi.sin().powi(-2);
    ///     // 0.0 to 0.99 in 0.01 increments
    ///     for x in (0..100).map(|i| i as f64 * 0.01) {
    ///         let p_k = x.ellip_e_inc(phi);
    ///         let rf = (c - x).carlson_rf(c, c - 1.0);
    ///         let rd = (c - 1.0).carlson_rd(c - x, c);
    ///         let carl = rf - x * rd / 3.0;
    ///         let rel_err = (carl - p_k) / p_k;
    ///         assert!(rel_err.abs() < 2e-15);
    ///     }
    ///
    ///     // -100.0 to -0.5 in 0.5 increments
    ///     for x in (-200..0).map(|i| i as f64 * 0.5) {
    ///         let p_k = x.ellip_e_inc(phi);
    ///         let rf = (c - x).carlson_rf(c, c - 1.0);
    ///         let rd = (c - 1.0).carlson_rd(c - x, c);
    ///         let carl = rf - x * rd / 3.0;
    ///         let rel_err = (carl - p_k) / p_k;
    ///         assert!(rel_err.abs() < 2e-15);
    ///     }
    /// }
    /// ``` 
    fn ellip_e_inc(&self, phi: Self) -> Self;

    /// Incomplete elliptic integral of the first kind
    ///
    /// # Description:
    ///
    /// Approximates the incomplete elliptic integral of the first kind 
    /// with `self` as the k squared term.
    ///
    /// [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_integral#Incomplete_elliptic_integral_of_the_first_kind)
    ///
    /// # Domain:
    ///
    /// Returns `NAN` if either paramter is `NAN`, if `self > 1.0`, or 
    /// if both paramters are infinite. Returns infinity if `self == 1.0`
    /// and phi is outside of `[-pi / 2, pi / 2]`
    ///
    /// # Examples:
    ///
    /// Association with Carlson's symmetric elliptic integral of the 
    /// first kind ([DLMF 19.25.5](https://dlmf.nist.gov/19.25#E5)):
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.1 to 1.5 in 0.1 increments
    /// for phi in (1..=15).map(|i| i as f64 * 0.1) {
    ///     let c = phi.sin().powi(-2);
    ///     // 0.0 to 0.99 in 0.01 increments
    ///     for x in (0..100).map(|i| i as f64 * 0.01) {
    ///         let p_k = x.ellip_k_inc(phi);
    ///         let rf = (c - x).carlson_rf(c, c - 1.0);
    ///         let rel_err = (rf - p_k) / p_k;
    ///         assert!(rel_err.abs() < 1e-14);
    ///     }
    ///
    ///     // -100.0 to -0.5 in 0.5 increments
    ///     for x in (-200..0).map(|i| i as f64 * 0.5) {
    ///         let p_k = x.ellip_k_inc(phi);
    ///         let rf = (c - x).carlson_rf(c, c - 1.0);
    ///         let rel_err = (rf - p_k) / p_k;
    ///         assert!(rel_err.abs() < 1e-14);
    ///     }
    /// }
    /// ``` 
    fn ellip_k_inc(&self, phi: Self) -> Self;

    /// Incomplete elliptic integral of the third kind
    ///
    /// # Description:
    ///
    /// Approximates the incomplete elliptic integral of the third kind 
    /// with `self` as the k squared term.
    ///
    /// [Wikipedia](https://en.wikipedia.org/wiki/Elliptic_integral#Incomplete_elliptic_integral_of_the_third_kind)
    /// ``` 
    fn ellip_pi_inc(&self, phi: Self, n: Self) -> Self;

    /// Degenerate Carlson symmetric elliptic integral of the first kind
    ///
    /// # Domain:
    ///
    /// Returns `NAN` if `x < 0.0` and infinity if `y == 0.0`.
    ///
    /// # Examples
    ///
    /// Specific values:
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // RC(0, 1/4) = pi
    /// assert_eq!((0.0_f64).carlson_rc(0.25), std::f64::consts::PI);
    ///
    /// // RC(9/4, 2) == ln(2)
    /// let a = (9.0 / 4.0).carlson_rc(2.0);
    /// let b = (2.0_f64).ln();
    /// assert!((a - b).abs() < 1e-15);
    ///
    /// // RC(1/4, -2) = ln(2) / 3
    /// let a = (0.25_f64).carlson_rc(-2.0);
    /// let b = (2.0_f64).ln() / 3.0;
    /// assert!((a - b).abs() < 1e-15);
    /// ```
    fn carlson_rc(&self, y: Self) -> Self;

    /// Carlson elliptic integral symmetric in only two variables
    ///
    /// # Domain
    ///
    /// Returns `NAN` if any paramters are less than zero and infinity if
    /// two or more of the paramters are equal to zero.
    fn carlson_rd(&self, y: Self, z: Self) -> Self;

    /// Carlson symmetric elliptic integral of the first kind
    ///
    /// # Domain
    ///
    /// Returns `NAN` if any paramters are less than zero and infinity if
    /// two or more of the paramters are equal to zero.
    /// 
    /// # Examples:
    ///
    /// Special values:
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // First lemniscate constant
    /// let a: f64 = 1.3110287771461;
    /// let a_calc = (0.0_f64).carlson_rf(1.0, 2.0);
    /// let err = (a - a_calc) / a;
    ///
    /// assert!(err < 1e-13);
    /// ```
    ///
    /// Degeneracy ([DLMF 19.20(i)](https://dlmf.nist.gov/19.20#i)):
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // 0.5 to 50 in increments of 0.5
    /// for y in (1..=100).map(|i| i as f64 * 0.5) {
    ///     // RF(0, y, y) = pi/2 / sqrt(y)
    ///     let a1 = (0.0_f64).carlson_rf(y, y);
    ///     let a2 = y.carlson_rf(y, 0.0);
    ///     let a3 = y.carlson_rf(0.0, y);
    ///     let b = std::f64::consts::FRAC_PI_2 / y.sqrt();
    ///     assert!((b - a1) / b < 1e-15);
    ///     assert!((b - a2) / b < 1e-15);
    ///     assert!((b - a3) / b < 1e-15);
    /// }
    /// ```
    fn carlson_rf(&self, y: Self, z: Self) -> Self;

    /// Carlson symmetric elliptic integral of the second kind
    ///
    /// # Domain
    ///
    /// Returns `NAN` if any paramters are less than zero and infinity if
    /// two or more of the paramters are equal to zero.
    ///
    /// # Examples
    ///
    /// Specific values:
    ///
    /// ```
    /// use spec_math::Ellip;
    ///
    /// // RG(0, 16, 16) = pi
    /// let a = (0.0_f64).carlson_rg(16.0, 16.0);
    /// let b = std::f64::consts::PI;
    /// assert!((a - b).abs() < 1e15);
    /// ```
    fn carlson_rg(&self, y: Self, z: Self) -> Self;

    /// Carlson symmetric elliptic integral of the third kind
    ///
    /// # Domain
    ///
    /// Returns `NAN` if `self`, `y`, or `z` are less than zero and infinity if
    /// two or more of `self`, `y`, and `z` are equal to zero or `p` is equal 
    /// to zero.
    fn carlson_rj(&self, y: Self, z: Self, p: Self) -> Self;
}

impl Ellip for f64 {
    fn ellip_j(&self, m: f64) -> EllipJOutput<f64> {
        //! Uses [`cephes64::ellpj`](crate::cephes64::ellpj)
        let res = crate::cephes64::ellpj(*self, m);
        EllipJOutput::<f64> {sn: res.0, cn: res.1, dn: res.2, phi: res.3}
    }
    fn ellip_e(&self) -> f64 {
        //! Uses [`cephes64::ellpe`](crate::cephes64::ellpe)
        crate::cephes64::ellpe(*self)
    }
    fn ellip_k(&self) -> f64 {
        //! Uses [`cephes64::ellpk`](crate::cephes64::ellpk) evaluated at `1.0 - self`
        crate::cephes64::ellpk(1.0 - *self)
    }
    fn ellip_pi(&self, n: f64) -> f64 {
        crate::misc::ellippi(n, *self)
    }
    fn ellip_e_inc(&self, phi: f64) -> f64 {
        //! Uses [`cephes64::ellie`](crate::cephes64::ellie)
        crate::cephes64::ellie(phi, *self)
    }
    fn ellip_k_inc(&self, phi: f64) -> f64 {
        //! Uses [`cephes64::ellik`](crate::cephes64::ellik)
        crate::cephes64::ellik(phi, *self)
    }
    fn ellip_pi_inc(&self, phi: f64, n: f64) -> f64 {
        crate::misc::ellippi_inc(phi, *self, n)
    }
    fn carlson_rc(&self, y: f64) -> f64 {
        crate::misc::elliprc(*self, y)
    }
    fn carlson_rd(&self, y: f64, z: f64) -> f64 {
        crate::misc::elliprd(*self, y, z)
    }
    fn carlson_rf(&self, y: f64, z: f64) -> f64 {
        crate::misc::elliprf(*self, y, z)
    }
    fn carlson_rg(&self, y: f64, z: f64) -> f64 {
        crate::misc::elliprg(*self, y, z)
    }
    fn carlson_rj(&self, y: f64, z: f64, p: f64) -> f64 {
        crate::misc::elliprj(*self, y, z, p)
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