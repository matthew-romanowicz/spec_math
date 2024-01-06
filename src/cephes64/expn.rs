/* Cephes Math Library Release 1.1:  March, 1985
* Copyright 1985 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

/* Sources
* [1] NIST, "The Digital Library of Mathematical Functions", dlmf.nist.gov
*/

/* Scipy changes:
* - 09-10-2016: improved asymptotic expansion for large n
*/

#![allow(clippy::excessive_precision)]

use crate::cephes64::consts::{MACHEP, MAXLOG};
use crate::cephes64::polevl::polevl;
use crate::cephes64::gamma::gamma;

const N_A: usize = 13;

const A0: [f64; 1] = [
    1.00000000000000000
];

const A1: [f64; 1] = [
    1.00000000000000000
];

const A2: [f64; 2] = [
    -2.00000000000000000, 
    1.00000000000000000
];

const A3: [f64; 3] = [
    6.00000000000000000, 
    -8.00000000000000000, 
    1.00000000000000000
];

const A4: [f64; 4] = [
    -24.0000000000000000, 
    58.0000000000000000, 
    -22.0000000000000000, 
    1.00000000000000000
];

const A5: [f64; 5] = [
    120.000000000000000, 
    -444.000000000000000, 
    328.000000000000000, 
    -52.0000000000000000, 
    1.00000000000000000
];

const A6: [f64; 6] = [
    -720.000000000000000, 
    3708.00000000000000, 
    -4400.00000000000000, 
    1452.00000000000000, 
    -114.000000000000000, 
    1.00000000000000000
];

const A7: [f64; 7] = [
    5040.00000000000000, 
    -33984.0000000000000, 
    58140.0000000000000, 
    -32120.0000000000000, 
    5610.00000000000000, 
    -240.000000000000000, 
    1.00000000000000000
];

const A8: [f64; 8] = [
    -40320.0000000000000, 
    341136.000000000000, 
    -785304.000000000000, 
    644020.000000000000, 
    -195800.000000000000, 
    19950.0000000000000, 
    -494.000000000000000, 
    1.00000000000000000
];

const A9: [f64; 9] = [
    362880.000000000000, 
    -3733920.00000000000, 
    11026296.0000000000, 
    -12440064.0000000000, 
    5765500.00000000000, 
    -1062500.00000000000, 
    67260.0000000000000, 
    -1004.00000000000000, 
    1.00000000000000000
];

const A10: [f64; 10] = [
    -3628800.00000000000, 
    44339040.0000000000, 
    -162186912.000000000, 
    238904904.000000000, 
    -155357384.000000000, 
    44765000.0000000000, 
    -5326160.00000000000, 
    218848.000000000000, 
    -2026.00000000000000, 
    1.00000000000000000
];

const A11: [f64; 11] = [
    39916800.0000000000, 
    -568356480.000000000, 
    2507481216.00000000, 
    -4642163952.00000000, 
    4002695088.00000000, 
    -1648384304.00000000, 
    314369720.000000000, 
    -25243904.0000000000, 
    695038.000000000000, 
    -4072.00000000000000, 
    1.00000000000000000
];

const A12: [f64; 12] = [
    -479001600.000000000, 
    7827719040.00000000, 
    -40788301824.0000000, 
    92199790224.0000000, 
    -101180433024.000000, 
    56041398784.0000000, 
    -15548960784.0000000, 
    2051482776.00000000, 
    -114876376.000000000, 
    2170626.00000000000, 
    -8166.00000000000000, 
    1.00000000000000000
];

const EUL: f64 = 0.57721566490153286060;
const BIG: f64 = 1.44115188075855872E+17;

pub fn expn(n: i32, x: f64) -> f64 {
    //! Exponential integral En
    //!
    //! ## DESCRIPTION:
    //!
    //! Evaluates the exponential integral
    //!
    #![doc=include_str!("expn.svg")]
    //!
    //! Both `n` and `x` must be nonnegative.
    //!
    //! The routine employs either a power series, a continued
    //! fraction, or an asymptotic formula depending on the
    //! relative values of `n` and `x`.
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
    //!     <td>10000</td>
    //!     <td>1.7e-15</td>
    //!     <td>3.6e-16</td>
    //! </tr>
    //!</table>

    if x.is_nan() || n < 0 || x < 0.0 {
        //sf_error("expn", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if x > MAXLOG {
        0.0
    } else if x == 0.0 {
        if n < 2 {
            //sf_error("expn", SF_ERROR_SINGULAR, NULL);
            f64::INFINITY
        } else {
            1.0 / (n as f64 - 1.0)
        }
    } else if n == 0 {
        (-x).exp() / x
    } else if n > 50 {
        /* Asymptotic expansion for large n, DLMF 8.20(ii) */
        expn_large_n(n, x)
    } else if x > 1.0 {
        /* Continued fraction, DLMF 8.19.17 */
        let mut k: i32 = 1;
        let mut pkm2: f64 = 1.0;
        let mut qkm2 = x;
        let mut pkm1: f64 = 1.0;
        let mut qkm1 = x + n as f64;
        let mut ans = pkm1 / qkm1;
        let mut t: f64;
        let mut yk: f64;
        let mut xk: f64;

        loop {
            k += 1;
            if k & 1 != 0 {
                yk = 1.0;
                xk = (n + (k - 1) / 2) as f64;
            } else {
                yk = x;
                xk = (k / 2) as f64;
            }
            let pk = pkm1 * yk + pkm2 * xk;
            let qk = qkm1 * yk + qkm2 * xk;
            if qk != 0.0 {
                let r = pk / qk;
                t = ((ans - r) / r).abs();
                ans = r;
            } else {
                t = 1.0;
            }
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            if pk.abs() > BIG {
                pkm2 /= BIG;
                pkm1 /= BIG;
                qkm2 /= BIG;
                qkm1 /= BIG;
            }

            if !(t > MACHEP) {
                break;
            }
        }

        ans * (-x).exp()

    } else {

        /* Power series expansion, DLMF 8.19.8 */
        let mut psi = -EUL - x.ln();
        for i in 1..n {
            psi += 1.0 / i as f64;
        }

        let z = -x;
        let mut xk: f64 = 0.0;
        let mut yk: f64 = 1.0;
        let mut pk = 1.0 - n as f64;
        let mut ans = if n == 1 {
            0.0
        } else {
            1.0 / pk
        };

        loop {
            xk += 1.0;
            yk *= z / xk;
            pk += 1.0;
            if pk != 0.0 {
                ans += yk / pk;
            }
            let t = if ans != 0.0 {
                (yk / ans).abs()
            } else {
                1.0
            };
            if !(t > MACHEP) {
                break;
            }
        }

        //let k = xk;
        let t = n;
        let r = n - 1;
        z.powi(r) * psi / gamma(t as f64) - ans
    }
}


/* Asymptotic expansion for large n, DLMF 8.20(ii) */
fn expn_large_n(n: i32, x: f64) -> f64
{
    let p = n as f64;
    let lambda = x / p;
    let multiplier = 1.0 / p / (lambda + 1.0) / (lambda + 1.0);
    let mut fac: f64 = 1.0;
    let mut res: f64 = 1.0; /* A[0] = 1 */

    let expfac = (-lambda * p).exp() / (lambda + 1.0) / p;
    if expfac == 0.0 {
        //sf_error("expn", SF_ERROR_UNDERFLOW, NULL);
        return 0.0;
    }

    /* Do the k = 1 term outside the loop since A[1] = 1 */
    fac *= multiplier;
    res += fac;

    for k in 2..N_A {
        fac *= multiplier;
        let term = match k {
            0 => fac * polevl(lambda, &A0, 0),
            1 => fac * polevl(lambda, &A1, 0),
            2 => fac * polevl(lambda, &A2, 1),
            3 => fac * polevl(lambda, &A3, 2),
            4 => fac * polevl(lambda, &A4, 3),
            5 => fac * polevl(lambda, &A5, 4),
            6 => fac * polevl(lambda, &A6, 5),
            7 => fac * polevl(lambda, &A7, 6),
            8 => fac * polevl(lambda, &A8, 7),
            9 => fac * polevl(lambda, &A9, 8),
            10 => fac * polevl(lambda, &A10, 9),
            11 => fac * polevl(lambda, &A11, 10),
            12 => fac * polevl(lambda, &A12, 11),
            _ => panic!("Impossible State")
        };
        //fac * polevl(lambda, &A[k], Adegs[k]);

        res += term;
        if term.abs() < MACHEP * res.abs() {
            break;
        }
    }

    expfac * res
}

#[cfg(test)]
mod expn_tests {
    use super::*;

    #[test]
    fn expn_trivials() {
        assert_eq!(expn(10, f64::NAN).is_nan(), true);
        assert_eq!(expn(10, -1e-10).is_nan(), true);//TODO: Should this be infinity?
        assert_eq!(expn(-1, 1.0).is_nan(), true);//TODO: Should this be infinity?
        assert_eq!(expn(0, 0.0), f64::INFINITY);
        assert_eq!(expn(1, 0.0), f64::INFINITY);
        assert_eq!(expn(1, MAXLOG + 1.0), 0.0);
    }

    #[test]
    fn expn_x_0() { // x = 0.0, n >= 2
        assert_eq!(expn(2, 0.0), 1.0);
        assert_eq!(expn(10, 0.0), 0.1111111111111111);
        assert_eq!(expn(1000, 0.0), 0.001001001001001001);
        assert_eq!(expn(2147483647, 0.0), 4.656612877414201e-10);
    }

    #[test]
    fn expn_n_0() { // n = 0, x != 0.0
        assert_eq!(expn(0, 1e-300), 9.999999999999999e+299);
        assert_eq!(expn(0, 1e-20), 1e20);
        assert_eq!(expn(0, 0.1), 9.048374180359595);
        assert_eq!(expn(0, 1.0), 0.36787944117144233);
        assert_eq!(expn(0, 100.0), 3.720075976020836e-46);
        assert_eq!(expn(0, 300.0), 1.7160667408040045e-133);
        assert_eq!(expn(0, 700.0), 1.4085252205371101e-307);
    }

    #[test]
    fn expn_x_large() { // 1.0 < x, 0 < n <= 50
        assert_eq!(expn(1, 1.0 + 1e-15), 0.21938393439551981);
        assert_eq!(expn(1, 100.0), 3.683597761682032e-46);
        assert_eq!(expn(1, 300.0), 1.71038427680451e-133);
        assert_eq!(expn(1, 700.0), 1.406518766234033e-307);

        assert_eq!(expn(20, 1.0 + 1e-15), 0.01834597120675585);
        assert_eq!(expn(20, 130.0), 2.322779091352467e-59);
        assert_eq!(expn(20, 350.0), 2.6840632099904437e-155);
        assert_eq!(expn(20, 700.0), 1.3694522116512558e-307);

        assert_eq!(expn(30, 1.0 + 1e-15), 0.01224860364044181);
        assert_eq!(expn(30, 110.0), 1.2081944304336529e-50);
        assert_eq!(expn(30, 310.0), 6.876125800674272e-138);
        assert_eq!(expn(30, 700.0), 1.350716463010757e-307);

        assert_eq!(expn(50, 1.0 + 1e-15), 0.007354589497231292);
        assert_eq!(expn(50, 90.0), 5.867698499297802e-42);
        assert_eq!(expn(50, 250.0), 8.902219116538637e-112);
        assert_eq!(expn(50, 700.0), 1.3147401151194962e-307);
    }

    #[test]
    fn expn_x_small() { // 0.0 < x <= 1.0, 0 < n <= 50
        assert_eq!(expn(1, 1e-300), 690.1983122333121);
        assert_eq!(expn(1, 1e-20), 45.474486194979384);
        assert_eq!(expn(1, 0.1), 1.8229239584193906);
        assert_eq!(expn(1, 1.0), 0.2193839343955205);

        assert_eq!(expn(20, 1e-300), 0.05263157894736842);
        assert_eq!(expn(20, 1e-20), 0.05263157894736842);
        assert_eq!(expn(20, 0.15), 0.04492621630828857);
        assert_eq!(expn(20, 1.0), 0.018345971206755872);

        assert_eq!(expn(40, 1e-300), 0.02564102564102564);
        assert_eq!(expn(40, 3e-20), 0.02564102564102564);
        assert_eq!(expn(40, 0.5), 0.015350163158933868);
        assert_eq!(expn(40, 1.0), 0.009191102220599899);

        assert_eq!(expn(50, 1e-300), 0.02040816326530612);
        assert_eq!(expn(50, 5e-20), 0.02040816326530612);
        assert_eq!(expn(50, 0.6), 0.011061998769259205);
        assert_eq!(expn(50, 1.0), 0.0073545894972313);
    }

    #[test]
    fn expn_asy() { // n > 50, x != 0.0
        assert_eq!(expn(51, 1e-20), 0.019999999999999997);
        assert_eq!(expn(51, 0.5), 0.012008107250881738);
        assert_eq!(expn(51, 1.0), 0.0072104970334842195);
        assert_eq!(expn(51, 5.0), 0.00012230294311891697);
        assert_eq!(expn(51, 50.0), 1.919154543178065e-24);
        assert_eq!(expn(51, 300.0), 1.4673281827194488e-133);
        assert_eq!(expn(51, 700.0), 1.3129914758465938e-307);

        assert_eq!(expn(100, 1e-20), 0.0101010101010101);
        assert_eq!(expn(100, 0.4), 0.006743386681914418);
        assert_eq!(expn(100, 1.5), 0.0022198677118589294);
        assert_eq!(expn(100, 8.0), 3.132958102010227e-06);
        assert_eq!(expn(100, 80.0), 1.0057841775490805e-37);
        assert_eq!(expn(100, 350.0), 2.2076620942609717e-155);
        assert_eq!(expn(100, 700.0), 1.2326517497355776e-307);

        assert_eq!(expn(2147483647, 1e-20), 4.656612877414201e-10);
        assert_eq!(expn(2147483647, 0.6), 2.5556033311961735e-10);
        assert_eq!(expn(2147483647, 1.1), 1.5500517740733187e-10);
        assert_eq!(expn(2147483647, 7.5), 2.575499791332095e-13);
        assert_eq!(expn(2147483647, 150.0), 3.3411641971543943e-75);
        assert_eq!(expn(2147483647, 320.0), 4.941239005955638e-149);
        assert_eq!(expn(2147483647, 700.0), 4.591268179e-314);
    }
}