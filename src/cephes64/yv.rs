/*
* Cephes Math Library Release 2.8: June, 2000
* Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

use crate::cephes64::consts::M_PI;
use crate::cephes64::yn::yn;
use crate::cephes64::jv::jv;

pub fn yv(v: f64, x: f64) -> f64 {
    //! Bessel function of noninteger order

    let n = v as isize;
    if n as f64 == v {
        yn(n, x)
    } else if v == v.floor() {
        /* Zero in denominator. */
        //sf_error("yv", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let t = M_PI * v;
        let y = (t.cos() * jv(v, x) - jv(-v, x)) / t.sin();

        if y.is_infinite() {
            if v > 0.0 {
                //sf_error("yv", SF_ERROR_OVERFLOW, NULL);
                -f64::INFINITY
            } else if v < -1e10 {
                /* Whether it's +inf or -inf is numerically ill-defined. */
                //sf_error("yv", SF_ERROR_DOMAIN, NULL);
                f64::NAN
            } else {
                y
            }
        } else {
            y
        }
    }
}

#[cfg(test)]
mod yv_tests {
    use super::*;

    #[test]
    fn yv_trivials() {
        assert_eq!(yv(f64::NAN, 10.5).is_nan(), true);
        assert_eq!(yv(1.0, f64::NAN).is_nan(), true);
        assert_eq!(yv(1.1, f64::INFINITY).is_nan(), true);
        assert_eq!(yv(1.1, -f64::INFINITY).is_nan(), true);
        assert_eq!(yv(f64::INFINITY, 10.3).is_nan(), true);
        assert_eq!(yv(-f64::INFINITY, -12.0).is_nan(), true);
    }

    #[test]
    fn yv_tests() { // Very few tests needed due to limited cases
        assert_eq!(yv(1.0, 10.5), 0.2337042283572685);
        assert_eq!(yv(15.0, 0.01), -9.093041849196235e44);
        assert_eq!(yv(-100.0, 100.0), -0.16692141141757652);

        assert_eq!(yv(-0.5, 23.1), -0.14860999401167396);
        assert_eq!(yv(12.9, 1e-4), -3.605334427668831e63);
        assert_eq!(yv(-156.3, 1.234e6), 0.0005846478307984981);
    }
}