
// Adapted from https://github.com/scipy/scipy/blob/d4c0eaf0070014bf3ccf1acc91f8a1c8768615f1/scipy/special/cephes/besselpoly.c

use crate::cephes64::gamma;

const EPS: f64 = 1.0e-17;

/// Weighted integral of the Bessel function of the first kind
pub fn besselpoly(a: f64, lambda: f64, nu: f64) -> f64 {

    let mut factor: isize = 0;
    let mut sum: f64 = 0.0;
    let mut nu = nu;

    /* Special handling for a = 0.0 */
    if a == 0.0 {
        if nu == 0.0 {
            1.0 / (lambda + 1.0)
        } else {
            0.0
        }
    } else {
        /* Special handling for negative and integer nu */
        if nu < 0.0 && nu.floor() == nu {
            nu = -nu;
            factor = (nu as isize) % 2;
        }
        let mut sm = (nu * a.ln()).exp()/(gamma(nu + 1.0) * (lambda + nu + 1.0));
        for m in 0..1000 {
            sum += sm;
            let sol = sm;
            let m_f = m as f64 + 1.0;
            let s = lambda + nu + 2.0 * m_f;
            sm *= -a * a * (s - 1.0) / ((nu + m_f) * m_f * (s + 1.0));
            if !(((sm - sol) / sm).abs() > EPS) { // Needs to be negative to account for NAN
                break;
            }
        }
        if factor == 0 {
            sum
        } else {
            -sum
        }
    }
}

#[cfg(test)]
mod besselpoly_tests {
    use super::*;

    #[test]
    fn besselpoly_trivials() {
        assert_eq!(besselpoly(0.0, 10.0, 1e-20), 0.0);
        assert_eq!(besselpoly(0.0, 10.0, -1e-20), 0.0);
        assert_eq!(besselpoly(0.0, -1.0, 0.0), f64::INFINITY);
    }

    #[test]
    fn besselpoly_a_0() {
        assert_eq!(besselpoly(0.0, 1e20, 0.0), 1e-20);
        assert_eq!(besselpoly(0.0, 10.0, 0.0), 0.09090909090909091);
        assert_eq!(besselpoly(0.0, 1.0, 0.0), 0.5);
        assert_eq!(besselpoly(0.0, 0.0, 0.0), 1.0);
        assert_eq!(besselpoly(0.0, -10.0, 0.0), -0.1111111111111111);
        assert_eq!(besselpoly(0.0, -1e20, 0.0), -1e-20);
    }

    #[test]
    fn besselpoly_neg_int_nu() {
        assert_eq!(besselpoly(0.0, 0.0, -1.0), 0.0);
        assert_eq!(besselpoly(10.0, 10.0, -1.0), 0.005443028118097967);
        assert_eq!(besselpoly(10.0, -10.0, -2.0), 469089.99298284587);
        assert_eq!(besselpoly(1e5, 10.0, -3.0), 3.792702761141381e+305);
        assert_eq!(besselpoly(1e5, -10.0, -5.0), -6.944444442361129e+31);
        assert_eq!(besselpoly(10.0, 1e10, -15.0), 8.120687623664639e-14);
        assert_eq!(besselpoly(10.0, -1e10, -25.0), 9.781165808953606e-13);
    }

    #[test]
    fn besselpoly_nominal() {
        assert_eq!(besselpoly(0.0, 0.0, 1.0), 0.0);
        assert_eq!(besselpoly(10.0, 10.0, 1.0), -0.005443028118097967);
        assert_eq!(besselpoly(10.0, -10.0, 2.0), 469089.99298284587);
        assert_eq!(besselpoly(1e5, 10.0, 3.0), -3.792702761141381e+305);
        assert_eq!(besselpoly(1e5, -10.0, 5.0), 6.944444442361129e+31);
        assert_eq!(besselpoly(10.0, 1e10, 15.0), -8.120687623664639e-14);
        assert_eq!(besselpoly(10.0, -1e10, 25.0), -9.781165808953606e-13);

        assert_eq!(besselpoly(0.0, 0.0, 0.5), 0.0);
        assert_eq!(besselpoly(1.0, -10.0, 0.5), -0.04633245834783164);
        assert_eq!(besselpoly(1.0, -10.0, -0.5), 0.034299575109575424);
        assert_eq!(besselpoly(10.0, 10.0, 1.5), -0.008167704119219283);
        assert_eq!(besselpoly(10.0, 10.0, -1.5), -0.0007316216802383398);
        assert_eq!(besselpoly(1e5, 10.0, 13.1), 3.0562740894880096e+305);
        assert_eq!(besselpoly(1e5, 10.0, -13.1), -6.420556839742189e+304);
        assert_eq!(besselpoly(10.0, 1e5, 25.8), 5.301091679335061e-08);
        assert_eq!(besselpoly(10.0, 1e5, -25.8), -2.183175315099723e-05);
        assert_eq!(besselpoly(50.0, 1e10, 50.8), 495433739.236827);
        assert_eq!(besselpoly(50.0, 1e10, -50.8), 4864756841.840624);
        assert_eq!(besselpoly(500.0, 1e20, 100.8), 5.519508017431162e+307);
        assert_eq!(besselpoly(500.0, 1e20, -100.8), -7.981400076198974e+307);
    }
}