/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1985, 1987 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#![allow(clippy::excessive_precision)]

/* Chebyshev coefficients for reciprocal Gamma function
* in interval 0 to 1.  Function is 1/(x Gamma(x)) - 1
*/

const R: [f64; 16] = [
    3.13173458231230000000E-17,
    -6.70718606477908000000E-16,
    2.20039078172259550000E-15,
    2.47691630348254132600E-13,
    -6.60074100411295197440E-12,
    5.13850186324226978840E-11,
    1.08965386454418662084E-9,
    -3.33964630686836942556E-8,
    2.68975996440595483619E-7,
    2.96001177518801696639E-6,
    -8.04814124978471142852E-5,
    4.16609138709688864714E-4,
    5.06579864028608725080E-3,
    -6.41925436109158228810E-2,
    -4.98558728684003594785E-3,
    1.27546015610523951063E-1
];

use crate::cephes64::consts::{M_PI, MAXLOG};
use crate::cephes64::sinpi::sinpi;
use crate::cephes64::chbevl::chbevl;
use crate::cephes64::gamma::lgam;


pub fn rgamma(x: f64) -> f64 {
    //! Reciprocal Gamma function
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns one divided by the Gamma function of the argument.
    //!
    //! The function is approximated by a Chebyshev expansion in
    //! the interval [0,1].  Range reduction is by recurrence
    //! for arguments between -34.034 and +34.84425627277176174.
    //! 0 is returned for positive arguments outside this
    //! range.  For arguments less than -34.034 the cosecant
    //! reflection formula is applied; lograrithms are employed
    //! to avoid unnecessary overflow.
    //!
    //! The reciprocal Gamma function has no singularities,
    //! but overflow and underflow may occur for large arguments.
    //! These conditions return either INFINITY or 0 with
    //! appropriate sign.
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
    //!     <td>-30, 30</td>
    //!     <td>30000</td>
    //!     <td>1.1e-15</td>
    //!     <td>2.0e-16</td>
    //! </tr>
    //!</table>
    //!
    //! For arguments less than -34.034 the peak error is on the
    //! order of 5e-15 (DEC), excepting overflow or underflow.

    if x > 34.84425627277176174 {
        (-lgam(x)).exp()
    } else if x < -34.034 {
        let w = -x;
        let mut z = sinpi(w);
        if z == 0.0 {
            return 0.0;
        } 

        let sign: f64 = if z < 0.0 {
            z = -z;
            1.0
        } else {
            -1.0
        };

        let y = (w * z).ln() - M_PI.ln() + lgam(w);
        if y < -MAXLOG {
            //sf_error(name, SF_ERROR_UNDERFLOW, NULL);
            sign * 0.0
        } else if y > MAXLOG {
            //sf_error(name, SF_ERROR_OVERFLOW, NULL);
            sign * f64::INFINITY
        } else {
            sign * y.exp()
        }
    } else {
        let mut z = 1.0;
        let mut w = x;
    
        while w > 1.0 {		/* Downward recurrence */
            w -= 1.0;
            z *= w;
        }

        while w < 0.0 {		/* Upward recurrence */
            z /= w;
            w += 1.0;
        }

        if w == 0.0 {		/* Nonpositive integer */
            0.0
        } else if w == 1.0 {		/* Other integer */
            1.0 / z
        } else {
            w * (1.0 + chbevl(4.0 * w - 2.0, &R, 16)) / z
        }
    }

}

#[cfg(test)]
mod rgamma_tests {
    use super::*;

    #[test]
    fn rgamma_trivials() {
        assert_eq!(rgamma(f64::INFINITY), 0.0);
        assert_eq!(rgamma(-f64::INFINITY).is_nan(), true);
    }

    #[test]
    fn rgamma_large() {
        assert_eq!(rgamma(34.9), 4.825650090664939e-39);
        assert_eq!(rgamma(35.0), 3.3871575355212e-39);
        assert_eq!(rgamma(50.0), 1.6439747083166107e-63);
        assert_eq!(rgamma(100.0), 1.0715102881254683e-156);
        assert_eq!(rgamma(150.0), 2.6254143103890303e-261);
        assert_eq!(rgamma(200.0), 0.0);
    }

    #[test]
    fn rgamma_large_neg() {
        assert_eq!(rgamma(-34.1), -4.138504207854243e+37);
        assert_eq!(rgamma(-34.5), -5.5398467958342174e+38);
        assert_eq!(rgamma(-34.9), -7.113807196676161e+38);
        assert_eq!(rgamma(-35.0), 0.0);
        //assert_eq!(rgamma(-35.0 - 1e-13), 6.769955188053906e+25);
        assert_eq!(rgamma(-100.1), -1.4557088784217443e+157);
        assert_eq!(rgamma(-100.5), -2.981789478307643e+158);
        assert_eq!(rgamma(-100.9), -5.841604372697383e+158);
        assert_eq!(rgamma(-100.0), 0.0);
    }

    #[test]
    fn rgamma_nominal() {
        assert_eq!(rgamma(-34.034), -1.130084754158773e+37);
        assert_eq!(rgamma(-34.0), 0.0);
        assert_eq!(rgamma(-29.9), 1.854094667574915e+31);
        assert_eq!(rgamma(-29.5), 1.5351121119984077e+31);
        assert_eq!(rgamma(-29.1), 1.2201946354944048e+30);
        assert_eq!(rgamma(-19.9), 1.7696341682144826e+17);
        assert_eq!(rgamma(-19.5), 1.7208605883889952e+17);
        assert_eq!(rgamma(-19.1), 1.6108179996116448e+16);
        assert_eq!(rgamma(-19.0), 0.0);
        assert_eq!(rgamma(-5.9), 58.77141310120965);
        assert_eq!(rgamma(-5.5), 91.63673001529571);
        assert_eq!(rgamma(-5.1), 14.0120239801757);
        assert_eq!(rgamma(-5.0), 0.0);
        assert_eq!(rgamma(-0.9), -0.09460233055005997);
        assert_eq!(rgamma(-0.5), -0.28209479177387814);
        assert_eq!(rgamma(-0.1), -0.09357787209128729);
        assert_eq!(rgamma(0.0), 0.0);
        assert_eq!(rgamma(0.5), 0.5641895835477563);
        assert_eq!(rgamma(1.0), 1.0);
        assert_eq!(rgamma(5.0), 0.041666666666666664);
        assert_eq!(rgamma(10.5), 8.823957200203801e-07);
        assert_eq!(rgamma(15.0), 1.1470745597729725e-11);
        assert_eq!(rgamma(20.1), 6.106410800727484e-18);
        assert_eq!(rgamma(25.9), 8.91098812597932e-26);
        assert_eq!(rgamma(30.0), 1.1309962886447718e-31);
        assert_eq!(rgamma(34.8442562727717617), 5.877471754111501e-39);
    }
}