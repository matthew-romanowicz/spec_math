/*                                                     j1.c
*
*     Bessel function of order one
*
*
*
* SYNOPSIS:
*
* double x, y, j1();
*
* y = j1( x );
*
*
*
* DESCRIPTION:
*
* Returns Bessel function of order one of the argument.
*
* The domain is divided into the intervals [0, 8] and
* (8, infinity). In the first interval a 24 term Chebyshev
* expansion is used. In the second, the asymptotic
* trigonometric representation is employed using two
* rational functions of degree 5/5.
*
*
*
* ACCURACY:
*
*                      Absolute error:
* arithmetic   domain      # trials      peak         rms
*    IEEE      0, 30       30000       2.6e-16     1.1e-16
*
*
*/
/*							y1.c
*
*	Bessel function of second kind of order one
*
*
*
* SYNOPSIS:
*
* double x, y, y1();
*
* y = y1( x );
*
*
*
* DESCRIPTION:
*
* Returns Bessel function of the second kind of order one
* of the argument.
*
* The domain is divided into the intervals [0, 8] and
* (8, infinity). In the first interval a 25 term Chebyshev
* expansion is used, and a call to j1() is required.
* In the second, the asymptotic trigonometric representation
* is employed using two rational functions of degree 5/5.
*
*
*
* ACCURACY:
*
*                      Absolute error:
* arithmetic   domain      # trials      peak         rms
*    IEEE      0, 30       30000       1.0e-15     1.3e-16
*
* (error criterion relative when |y1| > 1).
*
*/


/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

#![allow(clippy::excessive_precision)]
 
static RP: [f64; 4] = [
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
];

static RQ: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    6.20836478118054335476E2,
    2.56987256757748830383E5,
    8.35146791431949253037E7,
    2.21511595479792499675E10,
    4.74914122079991414898E12,
    7.84369607876235854894E14,
    8.95222336184627338078E16,
    5.32278620332680085395E18,
];

static PP: [f64; 7] = [
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
];

static PQ: [f64; 7] = [
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
];

static QP: [f64; 8] = [
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1,
];

static QQ: [f64; 7] = [
    /* 1.00000000000000000000E0, */
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
];

static YP: [f64; 6] = [
    1.26320474790178026440E9,
    -6.47355876379160291031E11,
    1.14509511541823727583E14,
    -8.12770255501325109621E15,
    2.02439475713594898196E17,
    -7.78877196265950026825E17,
];

static YQ: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    5.94301592346128195359E2,
    2.35564092943068577943E5,
    7.34811944459721705660E7,
    1.87601316108706159478E10,
    3.88231277496238566008E12,
    6.20557727146953693363E14,
    6.87141087355300489866E16,
    3.97270608116560655612E18,
];


static Z1: f64 = 1.46819706421238932572E1;
static Z2: f64 = 4.92184563216946036703E1;

use super::consts::{M_2_PI, THPIO4, SQ2OPI};
use super::polevl::{polevl, p1evl};

pub fn j1(x: f64) -> f64 {
    //! Bessel function of the first kind, order one
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns Bessel function of order one of the argument.
    //!
    //! The domain is divided into the intervals `[0, 8]` and
    //! `(8, infinity)`. In the first interval a 24 term Chebyshev
    //! expansion is used. In the second, the asymptotic
    //! trigonometric representation is employed using two
    //! rational functions of degree 5/5.
    //!
    //! ## ACCURACY:
    //! Absolute error:
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
    //!     <td>2.6e-16</td>
    //!     <td>1.1e-16</td>
    //! </tr>
    //!</table>
    let w = x;
    if x < 0.0 {
        -j1(-x)
    } else if w <= 5.0 {
        let z = x * x;
        let w = polevl(z, &RP, 3) / p1evl(z, &RQ, 8);
        w * x * (z - Z1) * (z - Z2)
    } else {
        let w = 5.0 / x;
        let z = w * w;
        let p = polevl(z, &PP, 6) / polevl(z, &PQ, 6);
        let q = polevl(z, &QP, 7) / p1evl(z, &QQ, 7);
        let xn = x - THPIO4;
        let p = p * xn.cos() - w * q * xn.sin();
        p * SQ2OPI / x.sqrt()
    }
}


pub fn y1(x: f64) -> f64 {
    //! Bessel function of second kind, order one
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns Bessel function of the second kind of order one
    //! of the argument.
    //!
    //! The domain is divided into the intervals `[0, 8]` and
    //! `(8, infinity)`. In the first interval a 25 term Chebyshev
    //! expansion is used, and a call to [`cephes64::j1`](crate::cephes64::j1) is required.
    //! In the second, the asymptotic trigonometric representation
    //! is employed using two rational functions of degree 5/5.
    //!
    //! ## ACCURACY:
    //!
    //! Absolute error, when `-1 <= y1(x) <= 1`; else relative error:
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
    //!     <td>1.0e-15</td>
    //!     <td>1.3e-16</td>
    //! </tr>
    //!</table>
    if x <= 5.0 {
        if x == 0.0 {
            //sf_error("y1", SF_ERROR_SINGULAR, NULL);
            -f64::INFINITY
        } else if x <= 0.0 {
            //sf_error("y1", SF_ERROR_DOMAIN, NULL);
            f64::NAN
        } else {
            let z = x * x;
            let w = x * (polevl(z, &YP, 5) / p1evl(z, &YQ, 8));
            w + M_2_PI * (j1(x) * x.ln() - 1.0 / x)
        }
    } else {
        let w = 5.0 / x;
        let z = w * w;
        let p = polevl(z, &PP, 6) / polevl(z, &PQ, 6);
        let q = polevl(z, &QP, 7) / p1evl(z, &QQ, 7);
        let xn = x - THPIO4;
        let p = p * xn.sin() + w * q * xn.cos();
        p * SQ2OPI / x.sqrt()
    }
}

#[cfg(test)]
mod j1_tests {
    use super::*;

    #[test]
    fn j1_small_x() {
        assert_eq!(j1(0.0), 0.0);

        assert_eq!(j1(1e-20), 5.0000000000000005e-21);
        assert_eq!(j1(1e-10), 5.000000000000001e-11);
        assert_eq!(j1(1e-7), 4.999999999999993e-08);
        assert_eq!(j1(1e-6), 4.999999999999375e-07);
        assert_eq!(j1(1e-5), 4.9999999999375e-06);
        assert_eq!(j1(1.01e-5), 5.049999999935606e-06);
        assert_eq!(j1(0.1), 0.049937526036242);
        assert_eq!(j1(1.0), 0.44005058574493355);
        assert_eq!(j1(2.0), 0.5767248077568734);
        assert_eq!(j1(3.0), 0.33905895852593654);
        assert_eq!(j1(4.0), -0.06604332802354912);
        assert_eq!(j1(5.0), -0.3275791375914653);

        assert_eq!(-j1(1e-20), -5.0000000000000005e-21);
        assert_eq!(-j1(1e-10), -5.000000000000001e-11);
        assert_eq!(-j1(1e-7), -4.999999999999993e-08);
        assert_eq!(-j1(1e-6), -4.999999999999375e-07);
        assert_eq!(-j1(1e-5), -4.9999999999375e-06);
        assert_eq!(-j1(1.01e-5), -5.049999999935606e-06);
        assert_eq!(-j1(0.1), -0.049937526036242);
        assert_eq!(-j1(1.0), -0.44005058574493355);
        assert_eq!(-j1(2.0), -0.5767248077568734);
        assert_eq!(-j1(3.0), -0.33905895852593654);
        assert_eq!(-j1(4.0), 0.06604332802354912);
        assert_eq!(-j1(5.0), 0.3275791375914653);
    }

    #[test]
    fn j1_large_x() {
        assert_eq!(j1(5.0 + 1e-16), -0.3275791375914653);
        assert_eq!(j1(6.0), -0.27668385812756563);
        assert_eq!(j1(10.0), 0.04347274616886141);
        assert_eq!(j1(100.0), -0.0771453520141123);
        assert_eq!(j1(1000.0), 0.00472831190708902);

        assert_eq!(j1(-5.0 - 1e-16), 0.3275791375914653);
        assert_eq!(j1(-6.0), 0.27668385812756563);
        assert_eq!(j1(-10.0), -0.04347274616886141);
        assert_eq!(j1(-100.0), 0.0771453520141123);
        assert_eq!(j1(-1000.0), -0.00472831190708902);

        // TODO: Accuracy reduces for very large x
        //assert_eq!(j1(1e10), -7.676506113818561e-06);
        //assert_eq!(j1(1e20), -5.449273779343996e-11);
    }
}

#[cfg(test)]
mod y1_tests {
    use super::*;

    #[test]
    fn y1_trivials() {
        assert_eq!(y1(0.0), -f64::INFINITY);
        assert_eq!(y1(-1e-20).is_nan(), true);
        assert_eq!(y1(-10.0).is_nan(), true);
        assert_eq!(y1(-f64::INFINITY).is_nan(), true);
    }

    #[test]
    fn y1_small_x() {
        assert_eq!(y1(1e-20), -6.3661977236758135e+19);
        assert_eq!(y1(1e-10), -6366197723.675814);
        assert_eq!(y1(0.1), -6.458951094702027);
        assert_eq!(y1(1.0), -0.7812128213002888);
        assert_eq!(y1(2.0), -0.10703243154093756);
        assert_eq!(y1(3.0), 0.3246744247918001);
        assert_eq!(y1(4.0), 0.3979257105571003);
        assert_eq!(y1(5.0), 0.14786314339122691);
    }

    #[test]
    fn y1_large_x() {
        assert_eq!(y1(5.0 + 1e-16), 0.14786314339122691);
        assert_eq!(y1(6.0), -0.17501034430039827);
        assert_eq!(y1(10.0), 0.24901542420695388);
        assert_eq!(y1(100.0), -0.02037231200275932);
        assert_eq!(y1(1000.0), -0.024784331292351868);
        // TODO: Accuracy reduces for very large x
        //assert_eq!(y1(1e10), -2.1755990258465346e-06);
        //assert_eq!(y1(1e20), -5.828155155322492e-11);
    }
}