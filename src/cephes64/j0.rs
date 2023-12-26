/*                                                     j0.c
*
*     Bessel function of order zero
*
*
*
* SYNOPSIS:
*
* double x, y, j0();
*
* y = j0( x );
*
*
*
* DESCRIPTION:
*
* Returns Bessel function of order zero of the argument.
*
* The domain is divided into the intervals [0, 5] and
* (5, infinity). In the first interval the following rational
* approximation is used:
*
*
*        2         2
* (w - r  ) (w - r  ) P (w) / Q (w)
*       1         2    3       8
*
*            2
* where w = x  and the two r's are zeros of the function.
*
* In the second interval, the Hankel asymptotic expansion
* is employed with two rational functions of degree 6/6
* and 7/7.
*
*
*
* ACCURACY:
*
*                      Absolute error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       60000       4.2e-16     1.1e-16
*
*/
/*							y0.c
*
*	Bessel function of the second kind, order zero
*
*
*
* SYNOPSIS:
*
* double x, y, y0();
*
* y = y0( x );
*
*
*
* DESCRIPTION:
*
* Returns Bessel function of the second kind, of order
* zero, of the argument.
*
* The domain is divided into the intervals [0, 5] and
* (5, infinity). In the first interval a rational approximation
* R(x) is employed to compute
*   y0(x)  = R(x)  +   2 * log(x) * j0(x) / M_PI.
* Thus a call to j0() is required.
*
* In the second interval, the Hankel asymptotic expansion
* is employed with two rational functions of degree 6/6
* and 7/7.
*
*
*
* ACCURACY:
*
*  Absolute error, when y0(x) < 1; else relative error:
*
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       30000       1.3e-15     1.6e-16
*
*/

/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/* Note: all coefficients satisfy the relative error criterion
* except YP, YQ which are designed for absolute error. */

#![allow(clippy::excessive_precision)]

static PP: [f64; 7] = [
    7.96936729297347051624E-4,
    8.28352392107440799803E-2,
    1.23953371646414299388E0,
    5.44725003058768775090E0,
    8.74716500199817011941E0,
    5.30324038235394892183E0,
    9.99999999999999997821E-1,
];

static PQ: [f64; 7] = [
    9.24408810558863637013E-4,
    8.56288474354474431428E-2,
    1.25352743901058953537E0,
    5.47097740330417105182E0,
    8.76190883237069594232E0,
    5.30605288235394617618E0,
    1.00000000000000000218E0,
];

static QP: [f64; 8] = [
    -1.13663838898469149931E-2,
    -1.28252718670509318512E0,
    -1.95539544257735972385E1,
    -9.32060152123768231369E1,
    -1.77681167980488050595E2,
    -1.47077505154951170175E2,
    -5.14105326766599330220E1,
    -6.05014350600728481186E0,
];

static QQ: [f64; 7] = [
    /*  1.00000000000000000000E0, */
    6.43178256118178023184E1,
    8.56430025976980587198E2,
    3.88240183605401609683E3,
    7.24046774195652478189E3,
    5.93072701187316984827E3,
    2.06209331660327847417E3,
    2.42005740240291393179E2,
];

static YP: [f64; 8] = [
    1.55924367855235737965E4,
    -1.46639295903971606143E7,
    5.43526477051876500413E9,
    -9.82136065717911466409E11,
    8.75906394395366999549E13,
    -3.46628303384729719441E15,
    4.42733268572569800351E16,
    -1.84950800436986690637E16,
];

static YQ: [f64; 7] = [
    /* 1.00000000000000000000E0, */
    1.04128353664259848412E3,
    6.26107330137134956842E5,
    2.68919633393814121987E8,
    8.64002487103935000337E10,
    2.02979612750105546709E13,
    3.17157752842975028269E15,
    2.50596256172653059228E17,
];

/*  5.783185962946784521175995758455807035071 */
static DR1: f64 = 5.78318596294678452118E0;

/* 30.47126234366208639907816317502275584842 */
static DR2: f64 = 3.04712623436620863991E1;

static RP: [f64; 4] = [
    -4.79443220978201773821E9,
    1.95617491946556577543E12,
    -2.49248344360967716204E14,
    9.70862251047306323952E15,
];

static RQ: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    4.99563147152651017219E2,
    1.73785401676374683123E5,
    4.84409658339962045305E7,
    1.11855537045356834862E10,
    2.11277520115489217587E12,
    3.10518229857422583814E14,
    3.18121955943204943306E16,
    1.71086294081043136091E18,
];

use super::consts::{M_2_PI, M_PI_4, SQ2OPI};
use super::polevl::{polevl, p1evl};

pub fn j0(x: f64) -> f64 {
    //! Bessel function of the first kind, order zero
    //! 
    //! ## DESCRIPTION:
    //!
    //! Returns Bessel function of the first kind, order zero, of the argument.
    //!
    //! The domain is divided into the intervals `[0, 5]` and `(5, infinity)`. In the first interval the following rational approximation is used:
    //!
    //!
    #![doc=include_str!("j0.svg")]
    //!
    //! where w = x<sup>2</sup>  and the two r's are zeros of the function.
    //!
    //! In the second interval, the Hankel asymptotic expansion is employed with two rational functions of degree 6/6 and 7/7.
    //!
    //! ## ACCURACY:
    //!
    //! Absolute Error:
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
    //!     <td>0, 30</td>
    //!     <td>60000</td>
    //!     <td>4.2e-16</td>
    //!     <td>1.1e-16</td>
    //! </tr>
    //!</table>

    let x = x.abs();

    if x <= 5.0 {
        let z = x * x;
        if x < 1.0e-5 {
            1.0 - z / 4.0
        } else {
            let p = (z - DR1) * (z - DR2);
            p * polevl(z, &RP, 3) / p1evl(z, &RQ, 8)
        }
    } else {
        let w = 5.0 / x;
        let q = 25.0 / (x * x);
        let p = polevl(q, &PP, 6) / polevl(q, &PQ, 6);
        let q = polevl(q, &QP, 7) / p1evl(q, &QQ, 7);
        let xn = x - M_PI_4;
        let p = p * xn.cos() - w * q * xn.sin();
        p * SQ2OPI / x.sqrt()
    }
}

/*                                                     y0() 2  */
/* Bessel function of second kind, order zero  */

/* Rational approximation coefficients YP[], YQ[] are used here.
* The function computed is  y0(x)  -  2 * log(x) * j0(x) / M_PI,
* whose value at x = 0 is  2 * ( log(0.5) + EUL ) / M_PI
* = 0.073804295108687225.
*/

pub fn y0(x: f64) -> f64
{
    //! Bessel function of the second kind, order zero
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns Bessel function of the second kind, of order
    //! zero, of the argument.
    //!
    //! The domain is divided into the intervals [0, 5] and
    //! (5, infinity). In the first interval a rational approximation
    //! R(x) is employed to compute
    //! 
    //! `y0(x)  = R(x)  +   2 * log(x) * j0(x) / M_PI`
    //!
    //! Thus a call to [`cephes64::j0`](crate::cephes64::j0) is required.
    //!
    //! In the second interval, the Hankel asymptotic expansion
    //! is employed with two rational functions of degree 6/6
    //! and 7/7.
    //!
    //! ## ACCURACY:
    //!
    //! Absolute error, when `y0(x) < 1`; else relative error:
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
    //!     <td>0, 30</td>
    //!     <td>30000</td>
    //!     <td>1.3e-15</td>
    //!     <td>1.6e-16</td>
    //! </tr>
    //!</table>

    if x <= 5.0 {
        if x == 0.0 {
            //sf_error("y0", SF_ERROR_SINGULAR, NULL);
            -f64::INFINITY
        } else if x < 0.0 {
            //sf_error("y0", SF_ERROR_DOMAIN, NULL);
            f64::NAN
        } else {
            let z = x * x;
            let w = polevl(z, &YP, 7) / p1evl(z, &YQ, 7);
            w + M_2_PI * x.ln() * j0(x)
        }
    } else {
        let w = 5.0 / x;
        let z = 25.0 / (x * x);
        let p = polevl(z, &PP, 6) / polevl(z, &PQ, 6);
        let q = polevl(z, &QP, 7) / p1evl(z, &QQ, 7);
        let xn = x - M_PI_4;
        let p = p * xn.sin() + w * q * xn.cos();
        p * SQ2OPI / x.sqrt()
    }


}

#[cfg(test)]
mod j0_tests {
    use super::*;

    #[test]
    fn j0_small_x() {
        assert_eq!(j0(0.0), 1.0);

        assert_eq!(j0(1e-20), 1.0);
        assert_eq!(j0(1e-10), 1.0);
        assert_eq!(j0(1e-7), 0.9999999999999974);
        assert_eq!(j0(1e-6), 0.99999999999975);
        assert_eq!(j0(1e-5), 0.9999999999750001);
        assert_eq!(j0(1.01e-5), 0.9999999999744974);
        assert_eq!(j0(0.1), 0.99750156206604);
        assert_eq!(j0(1.0), 0.7651976865579665);
        assert_eq!(j0(2.0), 0.22389077914123562);
        assert_eq!(j0(3.0), -0.2600519549019335);
        assert_eq!(j0(4.0), -0.3971498098638473);
        assert_eq!(j0(5.0), -0.1775967713143383);

        assert_eq!(j0(-1e-20), 1.0);
        assert_eq!(j0(-1e-10), 1.0);
        assert_eq!(j0(-1e-7), 0.9999999999999974);
        assert_eq!(j0(-1e-6), 0.99999999999975);
        assert_eq!(j0(-1e-5), 0.9999999999750001);
        assert_eq!(j0(-1.01e-5), 0.9999999999744974);
        assert_eq!(j0(-0.1), 0.99750156206604);
        assert_eq!(j0(-1.0), 0.7651976865579665);
        assert_eq!(j0(-2.0), 0.22389077914123562);
        assert_eq!(j0(-3.0), -0.2600519549019335);
        assert_eq!(j0(-4.0), -0.3971498098638473);
        assert_eq!(j0(-5.0), -0.1775967713143383);
    }

    #[test]
    fn j0_large_x() {
        assert_eq!(j0(5.0 + 1e-16), -0.1775967713143383);
        assert_eq!(j0(6.0), 0.15064525725099695);
        assert_eq!(j0(10.0), -0.24593576445134832);
        assert_eq!(j0(100.0), 0.01998585030422333);
        assert_eq!(j0(1000.0), 0.02478668615242003);

        assert_eq!(j0(-5.0 - 1e-16), -0.1775967713143383);
        assert_eq!(j0(-6.0), 0.15064525725099695);
        assert_eq!(j0(-10.0), -0.24593576445134832);
        assert_eq!(j0(-100.0), 0.01998585030422333);
        assert_eq!(j0(-1000.0), 0.02478668615242003);
        // TODO: Accuracy reduces for very large x
        //assert_eq!(j0(1e10), 2.175589294792473e-06);
        //assert_eq!(j0(1e20), -5.449273779343996e-11);
    }
}

#[cfg(test)]
mod y0_tests {
    use super::*;

    #[test]
    fn y0_trivials() {
        assert_eq!(y0(0.0), -f64::INFINITY);
        assert_eq!(y0(-1e-20).is_nan(), true);
        assert_eq!(y0(-10.0).is_nan(), true);
        assert_eq!(y0(-f64::INFINITY).is_nan(), true);
    }

    #[test]
    fn y0_small_x() {
        assert_eq!(y0(1e-20), -29.3912282502858);
        assert_eq!(y0(1e-10), -14.732516272697243);
        assert_eq!(y0(0.1), -1.5342386513503667);
        assert_eq!(y0(1.0), 0.08825696421567697);
        assert_eq!(y0(2.0), 0.5103756726497451);
        assert_eq!(y0(3.0), 0.3768500100127906);
        assert_eq!(y0(4.0), -0.016940739325064846);
        assert_eq!(y0(5.0), -0.30851762524903303);
    }

    #[test]
    fn y0_large_x() {
        assert_eq!(y0(5.0 + 1e-16), -0.30851762524903303);
        assert_eq!(y0(6.0), -0.28819468398157916);
        assert_eq!(y0(10.0), 0.05567116728359961);
        assert_eq!(y0(100.0), -0.0772443133650831);
        assert_eq!(y0(1000.0), 0.004715917977623586);
        // TODO: Accuracy reduces for very large x
        //assert_eq!(y0(1e10), -7.676508871690473e-06);
        //assert_eq!(y0(1e20), -5.828155155322492e-11);
    }
}