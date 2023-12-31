/*                                                     ellpk.c
 *
 *     Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m1, y, ellpk();
 *
 * y = ellpk( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used internally rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,1        30000       2.5e-16     6.8e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpk domain       x<0, x>1           0.0
 *
 */

/*                                                     ellpk.c */


/*
 * Cephes Math Library, Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#![allow(clippy::excessive_precision)]

use super::consts::MACHEP;
use super::polevl::polevl;

static P: [f64; 11] = [
    1.37982864606273237150E-4,
    2.28025724005875567385E-3,
    7.97404013220415179367E-3,
    9.85821379021226008714E-3,
    6.87489687449949877925E-3,
    6.18901033637687613229E-3,
    8.79078273952743772254E-3,
    1.49380448916805252718E-2,
    3.08851465246711995998E-2,
    9.65735902811690126535E-2,
    1.38629436111989062502E0
];
 
static Q: [f64; 11] = [
    2.94078955048598507511E-5,
    9.14184723865917226571E-4,
    5.94058303753167793257E-3,
    1.54850516649762399335E-2,
    2.39089602715924892727E-2,
    3.01204715227604046988E-2,
    3.73774314173823228969E-2,
    4.88280347570998239232E-2,
    7.03124996963957469739E-2,
    1.24999999999870820058E-1,
    4.99999999999999999821E-1
];

/// Natural log of 4.0
static C1: f64 = 1.3862943611198906188E0;

//$$K(m) = \int_{0}^{\pi / 2} \frac{d\theta}{\sqrt{1 - m\,\sin^2(\theta)}}$$

pub fn ellpk(x: f64) -> f64 {
    //! Complete elliptic integral of the first kind evaluated at `1.0 - x`
    //!
    //! ## Description:
    //!
    //! Approximates the integral
    //!
    #![doc=include_str!("ellpk.svg")]
    //!
    //! where `m = 1 - x`, using the approximation `P(x) - ln(x) * Q(x)`
    //!
    //! The argument `x` is used internally rather than `m` so that the logarithmic singularity at `m = 1` will be shifted to the origin; this preserves maximum accuracy.
    //!
    //! ## Accuracy:
    //!
    //! Relative Error:
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 1</td>
    //!     <td>30000</td>
    //!     <td>2.5e-16</td>
    //!     <td>6.8e-17</td>
    //! </tr>
    //!</table>
    //!
    //! ## Examples
    //! 
    //! ```rust
    //! use spec_math::cephes64::ellpk;
    //!
    //! let x = 0.0;
    //! 
    //! assert_eq!(ellpk(1.0 - x), std::f64::consts::PI * 0.5);
    //! ```
    if x < 0.0 {
        //sf_error("ellpk", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x > 1.0 {
        if x.is_infinite() {
            0.0
        } else {
            ellpk(1.0 / x) / x.sqrt()
        }
    } else if x > MACHEP {
        polevl(x, &P, 10) - x.ln() * polevl(x, &Q, 10)
    } else if x == 0.0 {
        //sf_error("ellpk", SF_ERROR_SINGULAR, NULL);
        f64::INFINITY
    } else {
        C1 - 0.5 * x.ln()
    }
}

#[cfg(test)]
mod ellpk_tests {
    use super::*;

    #[test]
    fn ellpk_trivials() {
        assert_eq!(ellpk(-f64::INFINITY).is_nan(), true);
        assert_eq!(ellpk(-10.0).is_nan(), true);
        assert_eq!(ellpk(-1.0).is_nan(), true);
        assert_eq!(ellpk(-1e-10).is_nan(), true);
        assert_eq!(ellpk(0.0), f64::INFINITY);
        assert_eq!(ellpk(f64::INFINITY), 0.0);
    }

    #[test]
    fn ellpk_large() {
        assert_eq!(ellpk(1.1), 1.5335928197134567);
        assert_eq!(ellpk(1.5), 1.415737208425956);
        assert_eq!(ellpk(2.0), 1.3110287771460598);
        assert_eq!(ellpk(10.0), 0.8152643095897213);
        assert_eq!(ellpk(100.0), 0.3695637362989875);
        assert_eq!(ellpk(1e10), 0.000128992198263876);
    }

    #[test]
    fn ellpk_medium() {
        assert_eq!(ellpk(1.15e-16), 19.73709413388468);
        assert_eq!(ellpk(1e-10), 12.8992198263876);
        assert_eq!(ellpk(0.5), 1.8540746773013719);
        assert_eq!(ellpk(0.9), 1.6124413487202192);
        assert_eq!(ellpk(1.0), 1.5707963267948966);
    }

    #[test]
    fn ellpk_small() {
        assert_eq!(ellpk(1.1e-16), 19.759320015170093); // Should actually be 19.759320015170094
        assert_eq!(ellpk(1.0e-20), 24.412145291060348); // Should actually round to 24.412145291060347
    }
}