/*                                                     k1.c
*
*     Modified Bessel function, third kind, order one
*
*
*
* SYNOPSIS:
*
* double x, y, k1();
*
* y = k1( x );
*
*
*
* DESCRIPTION:
*
* Computes the modified Bessel function of the third kind
* of order one of the argument.
*
* The range is partitioned into the two intervals [0,2] and
* (2, infinity).  Chebyshev polynomial expansions are employed
* in each interval.
*
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       30000       1.2e-15     1.6e-16
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* k1 domain          x <= 0          INFINITY
*
*/
/*							k1e.c
*
*	Modified Bessel function, third kind, order one,
*	exponentially scaled
*
*
*
* SYNOPSIS:
*
* double x, y, k1e();
*
* y = k1e( x );
*
*
*
* DESCRIPTION:
*
* Returns exponentially scaled modified Bessel function
* of the third kind of order one of the argument:
*
*      k1e(x) = exp(x) * k1(x).
*
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       30000       7.8e-16     1.2e-16
* See k1().
*
*/

/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#![allow(clippy::excessive_precision)]

/* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
* in the interval [0,2].
*
* lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
*/

const A: [f64; 11] = [
    -7.02386347938628759343E-18,
    -2.42744985051936593393E-15,
    -6.66690169419932900609E-13,
    -1.41148839263352776110E-10,
    -2.21338763073472585583E-8,
    -2.43340614156596823496E-6,
    -1.73028895751305206302E-4,
    -6.97572385963986435018E-3,
    -1.22611180822657148235E-1,
    -3.53155960776544875667E-1,
    1.52530022733894777053E0
];

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
* in the interval [2,infinity].
*
* lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
*/
const B: [f64; 25] = [
    -5.75674448366501715755E-18,
    1.79405087314755922667E-17,
    -5.68946255844285935196E-17,
    1.83809354436663880070E-16,
    -6.05704724837331885336E-16,
    2.03870316562433424052E-15,
    -7.01983709041831346144E-15,
    2.47715442448130437068E-14,
    -8.97670518232499435011E-14,
    3.34841966607842919884E-13,
    -1.28917396095102890680E-12,
    5.13963967348173025100E-12,
    -2.12996783842756842877E-11,
    9.21831518760500529508E-11,
    -4.19035475934189648750E-10,
    2.01504975519703286596E-9,
    -1.03457624656780970260E-8,
    5.74108412545004946722E-8,
    -3.50196060308781257119E-7,
    2.40648494783721712015E-6,
    -1.93619797416608296024E-5,
    1.95215518471351631108E-4,
    -2.85781685962277938680E-3,
    1.03923736576817238437E-1,
    2.72062619048444266945E0
];

use super::chbevl::chbevl;
use super::i1::i1;

pub fn k1(x: f64) -> f64
{
    //! Modified Bessel function, third kind, order one
    //!
    //! ## DESCRIPTION:
    //!
    //! Computes the modified Bessel function of the third kind
    //! of order one of the argument.
    //!
    //! The range is partitioned into the two intervals [0,2] and
    //! (2, infinity).  Chebyshev polynomial expansions are employed
    //! in each interval.
    //!
    //! ## ACCURACY:
    //!
    //! Relative error:
    //! 
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
    //!     <td>1.2e-15</td>
    //!     <td>1.6e-16</td>
    //! </tr>
    //!</table>

    if x == 0.0 {
        //sf_error("k1", SF_ERROR_SINGULAR, NULL);
        f64::INFINITY
    } else if x < 0.0 {
        //sf_error("k1", SF_ERROR_DOMAIN, NULL);
        f64:: NAN
    } else if x <= 2.0 {
        (0.5 * x).ln() * i1(x) + chbevl(x * x - 2.0, &A, 11) / x
    } else {
        (-x).exp() * chbevl(8.0 / x - 2.0, &B, 25) / x.sqrt()
    }
}




pub fn k1e(x: f64) -> f64
{
    //! Modified Bessel function, third kind, order one,
    //! exponentially scaled
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns exponentially scaled modified Bessel function
    //! of the third kind of order one of the argument:
    //!
    //! k1e(x) = exp(x) * k1(x).
    //!
    //! ## ACCURACY:
    //!
    //! Relative error:
    //!
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
    //!     <td>7.8e-16</td>
    //!     <td>1.2e-16</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::k1`](crate::cephes64::k1)

    if x == 0.0 {
        //sf_error("k1e", SF_ERROR_SINGULAR, NULL);
        f64::INFINITY
    } else if x < 0.0 {
        //sf_error("k1e", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x <= 2.0 {
        let y = (0.5 * x).ln() * i1(x) + chbevl(x * x - 2.0, &A, 11) / x;
        y * x.exp()
    } else {
        chbevl(8.0 / x - 2.0, &B, 25) / x.sqrt()
    }
}

#[cfg(test)]
mod k1_tests {
    use super::*;

    #[test]
    fn k1_trivials() {
        assert_eq!(k1(f64::INFINITY), 0.0);
        assert_eq!(k1(-f64::INFINITY).is_nan(), true);
        assert_eq!(k1(f64::NAN).is_nan(), true);
        assert_eq!(k1(-1e-20).is_nan(), true);
        assert_eq!(k1(-1.0).is_nan(), true);
        assert_eq!(k1(0.0), f64::INFINITY);
    }

    #[test]
    fn k1_small_x() {
        assert_eq!(k1(1e-20), 1e+20);

        assert_eq!(k1(1e-10), 10000000000.0);
        assert_eq!(k1(0.1), 9.853844780870606);
        assert_eq!(k1(0.5), 1.6564411200033007);
        assert_eq!(k1(1.0), 0.6019072301972346);
        assert_eq!(k1(1.5), 0.2773878004568438);
        assert_eq!(k1(2.0), 0.13986588181652246);
    }

    #[test]
    fn k1_large_x() {
        assert_eq!(k1(2.0 + 1e-16), 0.13986588181652246);
        assert_eq!(k1(3.0), 0.04015643112819419);
        assert_eq!(k1(10.0), 1.8648773453825585e-05);
        assert_eq!(k1(20.0), 5.883057969557038e-10);
        assert_eq!(k1(100.0), 4.67985373563691e-45);
        assert_eq!(k1(1000.0), 0.0);
    }
}

#[cfg(test)]
mod k1e_tests {
    use super::*;

    #[test]
    fn k0e_trivials() {
        assert_eq!(k1e(f64::INFINITY), 0.0);
        assert_eq!(k1e(-f64::INFINITY).is_nan(), true);
        assert_eq!(k1e(f64::NAN).is_nan(), true);
        assert_eq!(k1e(-1e-20).is_nan(), true);
        assert_eq!(k1e(-1.0).is_nan(), true);
        assert_eq!(k1e(0.0), f64::INFINITY);
    }

    #[test]
    fn k0e_small_x() {
        assert_eq!(k1e(1e-20), 1e+20);

        assert_eq!(k1e(1e-10), 10000000001.0);
        assert_eq!(k1e(0.1), 10.890182683049698);
        assert_eq!(k1e(0.5), 2.7310097082117855);
        assert_eq!(k1e(1.0), 1.636153486263258);
        assert_eq!(k1e(1.5), 1.2431658735525528);
        assert_eq!(k1e(2.0), 1.0334768470686888);
    }

    #[test]
    fn k0e_large_x() {
        assert_eq!(k1e(2.0 + 1e-16), 1.0334768470686888);
        assert_eq!(k1e(3.0), 0.806563480128787);
        assert_eq!(k1e(10.0), 0.4107665705957887);
        assert_eq!(k1e(20.0), 0.28542549694072644);
        assert_eq!(k1e(100.0), 0.12579995047957854);
        assert_eq!(k1e(1000.0), 0.03964813081296021);
        assert_eq!(k1e(1e10), 1.2533141373624994e-05);
        assert_eq!(k1e(1e20), 1.2533141373155e-10);
    }
}