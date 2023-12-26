/*                                                     k0.c
*
*     Modified Bessel function, third kind, order zero
*
*
*
* SYNOPSIS:
*
* double x, y, k0();
*
* y = k0( x );
*
*
*
* DESCRIPTION:
*
* Returns modified Bessel function of the third kind
* of order zero of the argument.
*
* The range is partitioned into the two intervals [0,8] and
* (8, infinity).  Chebyshev polynomial expansions are employed
* in each interval.
*
*
*
* ACCURACY:
*
* Tested at 2000 random points between 0 and 8.  Peak absolute
* error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       30000       1.2e-15     1.6e-16
*
* ERROR MESSAGES:
*
*   message         condition      value returned
*  K0 domain          x <= 0          INFINITY
*
*/
/*							k0e()
*
*	Modified Bessel function, third kind, order zero,
*	exponentially scaled
*
*
*
* SYNOPSIS:
*
* double x, y, k0e();
*
* y = k0e( x );
*
*
*
* DESCRIPTION:
*
* Returns exponentially scaled modified Bessel function
* of the third kind of order zero of the argument.
*
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       30000       1.4e-15     1.4e-16
* See k0().
*
*/

/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#![allow(clippy::excessive_precision)]

use crate::cephes64::chbevl;
use crate::cephes64::i0;
 
/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
* in the interval [0,2].  The odd order coefficients are all
* zero; only the even order coefficients are listed.
*
* lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
*/
 
const A: [f64; 10] = [
    1.37446543561352307156E-16,
    4.25981614279661018399E-14,
    1.03496952576338420167E-11,
    1.90451637722020886025E-9,
    2.53479107902614945675E-7,
    2.28621210311945178607E-5,
    1.26461541144692592338E-3,
    3.59799365153615016266E-2,
    3.44289899924628486886E-1,
    -5.35327393233902768720E-1
];

/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
* in the inverted interval [2,infinity].
*
* lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
*/
const B: [f64; 25] = [
    5.30043377268626276149E-18,
    -1.64758043015242134646E-17,
    5.21039150503902756861E-17,
    -1.67823109680541210385E-16,
    5.51205597852431940784E-16,
    -1.84859337734377901440E-15,
    6.34007647740507060557E-15,
    -2.22751332699166985548E-14,
    8.03289077536357521100E-14,
    -2.98009692317273043925E-13,
    1.14034058820847496303E-12,
    -4.51459788337394416547E-12,
    1.85594911495471785253E-11,
    -7.95748924447710747776E-11,
    3.57739728140030116597E-10,
    -1.69753450938905987466E-9,
    8.57403401741422608519E-9,
    -4.66048989768794782956E-8,
    2.76681363944501510342E-7,
    -1.83175552271911948767E-6,
    1.39498137188764993662E-5,
    -1.28495495816278026384E-4,
    1.56988388573005337491E-3,
    -3.14481013119645005427E-2,
    2.44030308206595545468E0
];

pub fn k0(x: f64) -> f64 {
    //! Modified Bessel function of the third kind, order zero
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns modified Bessel function of the third kind
    //! of order zero of the argument.
    //!
    //! The range is partitioned into the two intervals [0,8] and
    //! (8, infinity).  Chebyshev polynomial expansions are employed
    //! in each interval.
    //!
    //! ## ACCURACY:
    //!
    //! Tested at 2000 random points between 0 and 8.  Peak absolute
    //! error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
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
        //sf_error("k0", SF_ERROR_SINGULAR, NULL);
        f64::INFINITY
    } else if x < 0.0 {
        //sf_error("k0", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x <= 2.0 {
        chbevl(x * x - 2.0, &A, 10) - (0.5 * x).ln() * i0(x)
    } else {
        (-x).exp() * chbevl(8.0 / x - 2.0, &B, 25) / x.sqrt()
    }
}

pub fn k0e(x: f64) -> f64 {
    //! Modified Bessel function of the third kind, order zero, exponentially scaled
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns exponentially scaled modified Bessel function
    //! of the third kind of order zero of the argument.
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
    //!     <td>1.4e-15</td>
    //!     <td>1.4e-16</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::k0`](crate::cephes64::k0)

    if x == 0.0 {
        //sf_error("k0e", SF_ERROR_SINGULAR, NULL);
        f64::INFINITY
    } else if x < 0.0 {
        //sf_error("k0e", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x <= 2.0 {
        let y = chbevl(x * x - 2.0, &A, 10) - (0.5 * x).ln() * i0(x);
        y * x.exp()
    } else {
        chbevl(8.0 / x - 2.0, &B, 25) / x.sqrt()
    }
}

#[cfg(test)]
mod k0_tests {
    use super::*;

    #[test]
    fn k0_trivials() {
        assert_eq!(k0(f64::INFINITY), 0.0);
        assert_eq!(k0(-f64::INFINITY).is_nan(), true);
        assert_eq!(k0(f64::NAN).is_nan(), true);
        assert_eq!(k0(-1e-20).is_nan(), true);
        assert_eq!(k0(-1.0).is_nan(), true);
        assert_eq!(k0(0.0), f64::INFINITY);
    }

    #[test]
    fn k0_small_x() {
        assert_eq!(k0(1e-20), 46.16763337553933);

        assert_eq!(k0(1e-10), 23.141782445598864);
        assert_eq!(k0(0.1), 2.4270690247020164);
        assert_eq!(k0(0.5), 0.9244190712276656);
        assert_eq!(k0(1.0), 0.42102443824070823);
        assert_eq!(k0(1.5), 0.21380556264752565);
        assert_eq!(k0(2.0), 0.1138938727495334);
    }

    #[test]
    fn k0_large_x() {
        assert_eq!(k0(2.0 + 1e-16), 0.1138938727495334);
        assert_eq!(k0(3.0), 0.03473950438627925);
        assert_eq!(k0(10.0), 1.778006231616765e-05);
        assert_eq!(k0(20.0), 5.741237815336524e-10);
        assert_eq!(k0(100.0), 4.6566282291759025e-45);
        assert_eq!(k0(1000.0), 0.0);
    }
}

#[cfg(test)]
mod k0e_tests {
    use super::*;

    #[test]
    fn k0e_trivials() {
        assert_eq!(k0e(f64::INFINITY), 0.0);
        assert_eq!(k0e(-f64::INFINITY).is_nan(), true);
        assert_eq!(k0e(f64::NAN).is_nan(), true);
        assert_eq!(k0e(-1e-20).is_nan(), true);
        assert_eq!(k0e(-1.0).is_nan(), true);
        assert_eq!(k0e(0.0), f64::INFINITY);
    }

    #[test]
    fn k0e_small_x() {
        assert_eq!(k0e(1e-20), 46.16763337553933);

        assert_eq!(k0e(1e-10), 23.14178244791304);
        assert_eq!(k0e(0.1), 2.6823261022628944);
        assert_eq!(k0e(0.5), 1.5241093857739092);
        assert_eq!(k0e(1.0), 1.1444630798068947);
        assert_eq!(k0e(1.5), 0.9582100532948961);
        assert_eq!(k0e(2.0), 0.8415682150707712);
    }

    #[test]
    fn k0e_large_x() {
        assert_eq!(k0e(2.0 + 1e-16), 0.8415682150707712);
        assert_eq!(k0e(3.0), 0.6977615980438517);
        assert_eq!(k0e(10.0), 0.39163193443659866);
        assert_eq!(k0e(20.0), 0.2785448766571822);
        assert_eq!(k0e(100.0), 0.1251756216591266);
        assert_eq!(k0e(1000.0), 0.03962832160075422);
        assert_eq!(k0e(1e10), 1.2533141372998337e-05);
        assert_eq!(k0e(1e20), 1.2533141373155003e-10);
    }
}