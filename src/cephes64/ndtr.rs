/*                                                     ndtr.c
*
*     Normal distribution function
*
*
*
* SYNOPSIS:
*
* double x, y, ndtr();
*
* y = ndtr( x );
*
*
*
* DESCRIPTION:
*
* Returns the area under the Gaussian probability density
* function, integrated from minus infinity to x:
*
*                            x
*                             -
*                   1        | |          2
*    ndtr(x)  = ---------    |    exp( - t /2 ) dt
*               sqrt(2pi)  | |
*                           -
*                          -inf.
*
*             =  ( 1 + erf(z) ) / 2
*             =  erfc(z) / 2
*
* where z = x/sqrt(2). Computation is via the functions
* erf and erfc.
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE     -13,0        30000       3.4e-14     6.7e-15
*
*
* ERROR MESSAGES:
*
*   message         condition         value returned
* erfc underflow    x > 37.519379347       0.0
*
*/


/*
* Cephes Math Library Release 2.2:  June, 1992
* Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

//#include <float.h>		/* DBL_EPSILON */
//#include "mconf.h"

use crate::cephes64::consts::{MAXLOG, M_SQRT1_2};
use crate::cephes64::polevl::{polevl, p1evl};

const P: [f64; 9] = [
    2.46196981473530512524E-10,
    5.64189564831068821977E-1,
    7.46321056442269912687E0,
    4.86371970985681366614E1,
    1.96520832956077098242E2,
    5.26445194995477358631E2,
    9.34528527171957607540E2,
    1.02755188689515710272E3,
    5.57535335369399327526E2
];

const Q: [f64; 8] = [
    /* 1.00000000000000000000E0, */
    1.32281951154744992508E1,
    8.67072140885989742329E1,
    3.54937778887819891062E2,
    9.75708501743205489753E2,
    1.82390916687909736289E3,
    2.24633760818710981792E3,
    1.65666309194161350182E3,
    5.57535340817727675546E2
];

const R: [f64; 6] = [
    5.64189583547755073984E-1,
    1.27536670759978104416E0,
    5.01905042251180477414E0,
    6.16021097993053585195E0,
    7.40974269950448939160E0,
    2.97886665372100240670E0
];

const S: [f64; 6] = [
    /* 1.00000000000000000000E0, */
    2.26052863220117276590E0,
    9.39603524938001434673E0,
    1.20489539808096656605E1,
    1.70814450747565897222E1,
    9.60896809063285878198E0,
    3.36907645100081516050E0
];

const T: [f64; 5] = [
    9.60497373987051638749E0,
    9.00260197203842689217E1,
    2.23200534594684319226E3,
    7.00332514112805075473E3,
    5.55923013010394962768E4
];

const U: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    3.35617141647503099647E1,
    5.21357949780152679795E2,
    4.59432382970980127987E3,
    2.26290000613890934246E4,
    4.92673942608635921086E4
];

const UTHRESH: f64 = 37.519379347;


pub fn ndtr(a: f64) -> f64
{
    //double x, y, z;

    if a.is_nan() {
        //sf_error("ndtr", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    }

    let x = a * M_SQRT1_2;
    let z = x.abs();

    if z < M_SQRT1_2 {
        0.5 + 0.5 * erf(x)
    } else {
        if x > 0.0 {
            1.0 - 0.5 * erfc(z)
        } else {
            0.5 * erfc(z)
        }
    }
}

// $$\mathrm{erfc}(x)= 1 - \mathrm{erf}(x) = \frac{2}{\sqrt{\pi}}\int_{x}^{\infty}{\exp(-t^2)\,dt}$$

pub fn erfc(a: f64) -> f64 {
    //! Complementary error function
    //!
    //! ## DESCRIPTION:
    //!
    #![doc=include_str!("erfc.svg")]
    //!
    //! For small x, erfc(x) = 1 - erf(x); otherwise rational
    //! approximations are computed.
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
    //!     <td>0, 26.6417</td>
    //!     <td>30000</td>
    //!     <td>5.7e-14</td>
    //!     <td>1.5e-14</td>
    //! </tr>
    //!</table>

    if a.is_nan() {
        //sf_error("erfc", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    }

    let x = a.abs();

    if x < 1.0 {
        return 1.0 - erf(a);
    }

    let z = -a * a;

    if z < -MAXLOG {
        //sf_error("erfc", SF_ERROR_UNDERFLOW, NULL);
        if a < 0.0 {
            return 2.0;
        } else {
            return 0.0;
        }
    }

    let z = z.exp();

    let (p, q) = if x < 8.0 {
        (polevl(x, &P, 8), p1evl(x, &Q, 8))
    } else {
        (polevl(x, &R, 5), p1evl(x, &S, 6))
    };
    
    let mut y = (z * p) / q;

    if a < 0.0 {
        y = 2.0 - y;
    }

    if y != 0.0 {
        y
    } else if a < 0.0 { // Underflow conditions
        2.0
    } else {
        0.0
    }
}



pub fn erf(x: f64) -> f64 {
    //! Error function
    //!
    //! ## DESCRIPTION:
    //!
    //! The integral is
    //!
    #![doc=include_str!("erf.svg")]
    //!
    //! For 0 <= |x| < 1, erf(x) = x * P4(x^2)/Q5(x^2); otherwise
    //! erf(x) = 1 - erfc(x).
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
    //!     <td>0, 1</td>
    //!     <td>30000</td>
    //!     <td>3.7e-16</td>
    //!     <td>1.0e-16</td>
    //! </tr>
    //!</table>

    if x.is_nan() {
        //sf_error("erf", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x < 0.0 {
        -erf(-x)
    } else if x > 1.0 { 
        // Original had x.abs() > 1.0, but since negative x is 
        // captured in previous case, abs function was removed
        1.0 - erfc(x)
    } else {
        let z = x * x;
        x * polevl(z, &T, 4) / p1evl(z, &U, 5)
    }
    

    
}

#[cfg(test)]
mod erfc_tests {
    use super::*;

    #[test]
    fn erfc_trivials() {
        assert_eq!(erfc(f64::NAN).is_nan(), true);
        assert_eq!(erfc(f64::INFINITY), 0.0);
        assert_eq!(erfc(-f64::INFINITY), 2.0);
    }

    #[test]
    fn erfc_small() {
        assert_eq!(erfc(0.0), 1.0);

        assert_eq!(erfc(0.1), 0.8875370839817152);
        assert_eq!(erfc(0.5), 0.4795001221869535);
        assert_eq!(erfc(1.0), 0.15729920705028516);

        assert_eq!(erfc(-0.1), 1.1124629160182848);
        assert_eq!(erfc(-0.5), 1.5204998778130465);
        assert_eq!(erfc(-1.0), 1.8427007929497148);
    }

    #[test]
    fn erfc_large() {
        assert_eq!(erfc(1.1), 0.11979493042591832);
        assert_eq!(erfc(5.0), 1.5374597944280347e-12);
        assert_eq!(erfc(10.0), 2.0884875837625446e-45);
        assert_eq!(erfc(15.0), 7.212994172451209e-100);
        assert_eq!(erfc(25.0), 8.300172571196523e-274);
        assert_eq!(erfc(1e20), 0.0);

        assert_eq!(erfc(-1.1), 1.8802050695740817);
        assert_eq!(erfc(-5.0), 1.9999999999984626);
        assert_eq!(erfc(-1e20), 2.0);
    }
}

#[cfg(test)]
mod erf_tests {
    use super::*;

    #[test]
    fn erf_trivials() {
        assert_eq!(erf(f64::NAN).is_nan(), true);
        assert_eq!(erf(f64::INFINITY), 1.0);
        assert_eq!(erf(-f64::INFINITY), -1.0);
    }

    #[test]
    fn erf_small() {
        assert_eq!(erf(0.0), 0.0);

        assert_eq!(erf(0.1), 0.1124629160182849);
        assert_eq!(erf(0.5), 0.5204998778130465);
        assert_eq!(erf(1.0), 0.8427007929497148);

        assert_eq!(erf(-0.1), -0.1124629160182849);
        assert_eq!(erf(-0.5), -0.5204998778130465);
        assert_eq!(erf(-1.0), -0.8427007929497148);
    }

    #[test]
    fn erf_large() {
        assert_eq!(erf(1.1), 0.8802050695740817);
        assert_eq!(erf(5.0), 0.9999999999984626);
        assert_eq!(erf(1e20), 1.0);

        assert_eq!(erf(-1.1), -0.8802050695740817);
        assert_eq!(erf(-5.0), -0.9999999999984626);
        assert_eq!(erf(-1e20), -1.0);
    }
}