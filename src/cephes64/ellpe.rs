/*                                                     ellpe.c
 *
 *     Complete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m, y, ellpe();
 *
 * y = ellpe( m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *            pi/2
 *             -
 *            | |                 2
 * E(m)  =    |    sqrt( 1 - m sin t ) dt
 *          | |
 *           -
 *            0
 *
 * Where m = 1 - m1, using the approximation
 *
 *      P(x)  -  x log x Q(x).
 *
 * Though there are no singularities, the argument m1 is used
 * internally rather than m for compatibility with ellpk().
 *
 * E(1) = 1; E(0) = pi/2.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0, 1       10000       2.1e-16     7.3e-17
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpe domain      x<0, x>1            0.0
 *
 */

/*                                                     ellpe.c         */

/* Elliptic integral of second kind */

/*
 * Cephes Math Library, Release 2.1:  February, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 *
 * Feb, 2002:  altered by Travis Oliphant
 * so that it is called with argument m
 * (which gets immediately converted to m1 = 1-m)
 */

#![allow(clippy::excessive_precision)]

use super::consts::M_PI_2;
use super::polevl::polevl;

static P: [f64; 11] = [
    1.53552577301013293365E-4,
    2.50888492163602060990E-3,
    8.68786816565889628429E-3,
    1.07350949056076193403E-2,
    7.77395492516787092951E-3,
    7.58395289413514708519E-3,
    1.15688436810574127319E-2,
    2.18317996015557253103E-2,
    5.68051945617860553470E-2,
    4.43147180560990850618E-1,
    1.00000000000000000299E0
];

static Q: [f64; 10] = [
    3.27954898576485872656E-5,
    1.00962792679356715133E-3,
    6.50609489976927491433E-3,
    1.68862163993311317300E-2,
    2.61769742454493659583E-2,
    3.34833904888224918614E-2,
    4.27180926518931511717E-2,
    5.85936634471101055642E-2,
    9.37499997197644278445E-2,
    2.49999999999888314361E-1
];

//$$E(m) = \int_{0}^{\pi / 2} \sqrt{1 - m\,\sin^2(\theta)} \,d\theta$$

pub fn ellpe(x: f64) -> f64 {
//! Complete elliptic integral of the second kind
//!
//! ## Description
//!
//! Approximates the integral 
//!
#![doc=include_str!("ellpe.svg")]
//!
//! Where `m = 1 - x`, using the approximation `P(x) - x * ln(x) * Q(x)`
//!
//! Though there are no singularities, the argument `x` is used
//! internally rather than m for compatibility with [`ellpk`](super::ellpk::ellpk).
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
//!     <td>10000</td>
//!     <td>2.1e-16</td>
//!     <td>7.3e-17</td>
//! </tr>
//!</table>
//!
//! ## Examples
//! 
//! ```rust
//! use spec_math::cephes64::ellpe;
//! 
//! assert_eq!(ellpe(1.0), 1.0);
//!
//! assert_eq!(ellpe(0.0), std::f64::consts::PI * 0.5);
//! ```
    let x = 1.0 - x;
    if x <= 0.0 {
	    if x == 0.0 {
            1.0
        } else {
            //sf_error("ellpe", SF_ERROR_DOMAIN, NULL);
            f64::NAN
        }
    } else if x > 1.0 {
        ellpe(1.0 - 1.0 / x) * x.sqrt()
    } else {
        polevl(x, &P, 10) - x.ln() * (x * polevl(x, &Q, 9))
    }
}

#[cfg(test)]
mod ellpe_tests {
    use super::*;

    // #[test]
    // fn timing_tests() {
    //     let mut s = 0.0;
    //     let now = std::time::Instant::now();
    //     for i in 0..100000 {
    //         let x = (i as f64) / 100000.0;
    //         s += ellpe(x);
    //     }
    //     println!("{} {}", now.elapsed().as_micros(), s);
    //     let mut s = 0.0;
    //     let now = std::time::Instant::now();
    //     for i in 0..100000 {
    //         let x = -(i as f64) / 10.0;
    //         s += ellpe(x);
    //     }
    //     println!("{} {}", now.elapsed().as_micros(), s);
    // }

    #[test]
    fn ellpe_trivials() {
        assert_eq!(ellpe(-f64::INFINITY), f64::INFINITY);
        assert_eq!(ellpe(f64::INFINITY).is_nan(), true);
        assert_eq!(ellpe(10.0).is_nan(), true);
        assert_eq!(ellpe(1.01).is_nan(), true);
        assert_eq!(ellpe(1.0), 1.0);
    }

    #[test]
    fn ellpe_negative() {
        assert_eq!(ellpe(-1e-20), 1.5707963267948966);
        assert_eq!(ellpe(-1e-10), 1.5707963268341665);
        assert_eq!(ellpe(-0.1), 1.6093590249375296);
        assert_eq!(ellpe(-1.0), 1.9100988945138562);
        assert_eq!(ellpe(-10.0), 3.639138038417769);
        assert_eq!(ellpe(-1e10), 100000.00006699612);
        assert_eq!(ellpe(-1e20), 10000000000.0);
    }

    #[test]
    fn ellpe_small() {
        assert_eq!(ellpe(0.0), 1.5707963267948966);
        assert_eq!(ellpe(1e-20), 1.5707963267948966);
        assert_eq!(ellpe(1e-10), 1.5707963267556266);
        assert_eq!(ellpe(0.1), 1.5307576368977633);
        assert_eq!(ellpe(0.5), 1.3506438810476755);
        assert_eq!(ellpe(0.9), 1.1047747327040733);
        assert_eq!(ellpe(1.0 - 1e-10), 1.0000000006199612);
        assert_eq!(ellpe(1.0 - 1e-15), 1.000000000000009);
    }
}