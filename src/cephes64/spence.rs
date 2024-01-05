/*
* Cephes Math Library Release 2.1:  January, 1989
* Copyright 1985, 1987, 1989 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#![allow(clippy::excessive_precision)]


const A: [f64; 8] = [
    4.65128586073990045278E-5,
    7.31589045238094711071E-3,
    1.33847639578309018650E-1,
    8.79691311754530315341E-1,
    2.71149851196553469920E0,
    4.25697156008121755724E0,
    3.29771340985225106936E0,
    1.00000000000000000126E0,
];

const B: [f64; 8] = [
    6.90990488912553276999E-4,
    2.54043763932544379113E-2,
    2.82974860602568089943E-1,
    1.41172597751831069617E0,
    3.63800533345137075418E0,
    5.03278880143316990390E0,
    3.54771340985225096217E0,
    9.99999999999999998740E-1,
];

use crate::cephes64::consts::M_PI;
use crate::cephes64::polevl::polevl;

pub fn spence(x: f64) -> f64 {
    //! Dilogarithm
    //!
    //! ## DESCRIPTION:
    //!
    //! Computes the integral
    //!
    #![doc=include_str!("spence.svg")]
    //!
    //! for x >= 0.  A rational approximation gives the integral in
    //! the interval (0.5, 1.5).  Transformation formulas for 1/x
    //! and 1-x are employed outside the basic expansion range.
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
    //!     <td>0, 4</td>
    //!     <td>30000</td>
    //!     <td>3.9e-15</td>
    //!     <td>5.4e-16</td>
    //! </tr>
    //!</table>

    if x < 0.0 {
        //sf_error("spence", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x == 1.0 {
        0.0
    } else if x == 0.0 {
        M_PI * M_PI / 6.0
    } else {

        let mut flag = 0;
        let mut x = x;

        if x > 2.0 {
            x = 1.0 / x;
            flag |= 2;
        }

        let w = if x > 1.5 {
            flag |= 2;
            1.0 / x - 1.0
        } else if x < 0.5 {
            flag |= 1;
            -x
        } else {
            x - 1.0
        };


        let mut y = -w * polevl(w, &A, 7) / polevl(w, &B, 7);

        if flag & 1 != 0 {
            y = (M_PI * M_PI) / 6.0 - x.ln() * (1.0 - x).ln() - y;
        }

        if flag & 2 != 0 {
            let z = x.ln();
            y = -0.5 * z * z - y;
        }

        y
    }
}

#[cfg(test)]
mod spence_tests {
    use super::*;

    #[test]
    fn spence_trivials() {
        assert_eq!(spence(-1e-20).is_nan(), true);
        assert_eq!(spence(0.0), M_PI * M_PI / 6.0);
        assert_eq!(spence(1.0), 0.0);
    }

    #[test]
    fn spence_asy() { // x > 2.0
        assert_eq!(spence(2.0 + 1e-15), -0.822467033424113);
        assert_eq!(spence(2.5), -1.1473806603755703);
        assert_eq!(spence(3.0), -1.436746366883681);
        assert_eq!(spence(10.0), -3.950663778244157);
        assert_eq!(spence(100.0), -12.192421669033168);
        assert_eq!(spence(1e5), -67.91853531797291);
        assert_eq!(spence(1e10), -266.73983958836556);
        assert_eq!(spence(1e50), -6629.017572164846);
        assert_eq!(spence(1e150), -59647.998676948824);
        assert_eq!(spence(1e300), -238587.05990559477);
    }

    #[test]
    fn spence_large() { // 1.5 < x <= 2.0
        assert_eq!(spence(1.5 + 1e-15), -0.44841420692364703);
        assert_eq!(spence(1.6), -0.5281071740446668);
        assert_eq!(spence(1.7), -0.6051584023377058);
        assert_eq!(spence(1.8), -0.6797815878346817);
        assert_eq!(spence(1.9), -0.7521631792172623);
        assert_eq!(spence(2.0), -0.8224670334241142);
    }

    #[test]
    fn spence_medium() { // 0.5 <= x <= 1.5
        assert_eq!(spence(0.5), 0.5822405264650135);
        assert_eq!(spence(0.6), 0.44928297447128196);
        assert_eq!(spence(0.7), 0.3261295100754763);
        assert_eq!(spence(0.8), 0.21100377543970478);
        assert_eq!(spence(0.9), 0.1026177910993911);
        assert_eq!(spence(1.0), 0.0);
    }

    #[test]
    fn spence_small() { // 0.0 < x < 0.5
        assert_eq!(spence(0.5 - 1e-15), 0.5822405264650127);
        assert_eq!(spence(0.4), 0.727586307716333);
        assert_eq!(spence(0.3), 0.8893776242860385);
        assert_eq!(spence(0.2), 1.0747946000082482);
        assert_eq!(spence(0.1), 1.2997147230049588);
        assert_eq!(spence(1e-15), 1.6449340668481909);
    }
}