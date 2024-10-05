/*
* Pochhammer symbol (a)_m = gamma(a + m) / gamma(a)
*/
use crate::cephes64::lgam;
use crate::misc::gammasgn;

fn is_nonpos_int(x: f64) -> bool
{
    x <= 0.0 && x == x.ceil() && x.abs() < 1e13
}

pub fn poch(a: f64, mut m: f64) -> f64 {
    //! Pochhammer Symbol (a)_m
    //!
    //! # Domain:
    //!
    //! Returns `NAN` if `a` and `m` sum to 0.0 or a negative integer.

    let mut r: f64 = 1.0;

    /*
    * 1. Reduce magnitude of `m` to |m| < 1 by using recurrence relations.
    *
    * This may end up in over/underflow, but then the function itself either
    * diverges or goes to zero. In case the remainder goes to the opposite
    * direction, we end up returning 0*INF = NAN, which is OK.
    */

    /* Recurse down */
    while m >= 1.0 {
        if a + m == 1.0 {
            break;
        }
        m -= 1.0;
        r *= a + m;
        if !r.is_finite() || r == 0.0 {
            break;
        }
    }

    /* Recurse up */
    while m <= -1.0 {
        if a + m == 0.0 {
            break;
        }
        r /= a + m;
        m += 1.0;
        if !r.is_finite() || r == 0.0 {
            break;
        }
    }

    /*
    * 2. Evaluate function with reduced `m`
    *
    * Now either `m` is not big, or the `r` product has over/underflown.
    * If so, the function itself does similarly.
    */

    if m == 0.0 {
        /* Easy case */
        r
    } else if a > 1e4 && m.abs() <= 1.0 {
        /* Avoid loss of precision */
        r * a.powf(m) * (
            1.0
            + m*(m-1.0)/(2.0*a)
            + m*(m-1.0)*(m-2.0)*(3.0*m-1.0)/(24.0*a*a)
            + m*m*(m-1.0)*(m-1.0)*(m-2.0)*(m-3.0)/(48.0*a*a*a)
            )
    } else if is_nonpos_int(a + m) && !is_nonpos_int(a) && a + m != m {
        /* Check for infinity */
        // M Romanowicz: this should return NAN, not infinity
        //return f64::INFINITY;
        f64::NAN
    } else if !is_nonpos_int(a + m) && is_nonpos_int(a) {
        /* Check for zero */
        0.0
    } else {
        r * (lgam(a + m) - lgam(a)).exp() * gammasgn(a + m) * gammasgn(a)
    }
}

#[cfg(test)]
mod poch_tests {
    use super::*;

    #[test]
    fn poch_singularities() {
        assert!(poch(-10.5, 10.5).is_nan());
        assert!(poch(100.0, -110.0).is_nan());
        assert!(poch(-10.5, -0.5).is_nan());
        assert!(poch(-10.5, -1.5).is_nan());
        assert_eq!(poch(-10.5, -1.5 - 1e-10), 79074969.54406857);
        assert_eq!(poch(-10.5, -1.5 + 1e-10), -79074969.58401728);
    }

    #[test]
    fn poch_pos_m() { // m <= -1
        assert_eq!(poch(-1e100, 1.0), -1e100);
        assert_eq!(poch(-1e10, 1.0), -1e10);
        assert_eq!(poch(-30.0, 1.0), -30.0);
        assert_eq!(poch(0.0, 1.0), 0.0);
        assert_eq!(poch(30.0, 1.0), 30.0);
        assert_eq!(poch(1e10, 1.0), 1e10);
        assert_eq!(poch(1e100, 1.0), 1e100);

        assert!(poch(-1e100, 1.5).is_nan());
        assert_eq!(poch(-1e10, 1.5), 0.0);
        assert_eq!(poch(-30.0, 1.5), 0.0);
        assert_eq!(poch(0.0, 1.5), 0.0);
        assert_eq!(poch(30.0, 1.5), 166.3607961576209);
        assert_eq!(poch(1e10, 1.5), 1000000000037500.0);
        assert_eq!(poch(1e100, 1.5), 1.0000000000000002e+150);

        assert!(poch(-1e100, 100.1).is_nan());
        assert_eq!(poch(-1e10, 100.2), 0.0);
        assert_eq!(poch(-30.0, 100.3), 0.0);
        assert_eq!(poch(0.0, 100.4), 0.0);
        assert_eq!(poch(30.0, 100.5), 6.4086281038236536e+187);
        assert_eq!(poch(1e10, 100.6), f64::INFINITY);
        assert_eq!(poch(1e100, 100.7), f64::INFINITY);
    }

    #[test]
    fn poch_neg_m() { // m <= -1
        assert_eq!(poch(-1e100, -1.0), -1e-100);
        assert_eq!(poch(-1e10, -1.0), -9.999999999e-11);
        assert_eq!(poch(-30.0, -1.0), -0.03225806451612903);
        assert_eq!(poch(0.0, -1.0), -1.0);
        assert_eq!(poch(30.0, -1.0), 0.034482758620689655);
        assert_eq!(poch(1e10, -1.0), 1.0000000001e-10);
        assert_eq!(poch(1e100, -1.0), 1e-100);

        assert!(poch(-1e100, -1.5).is_nan());
        assert_eq!(poch(-1e10, -1.5), 0.0);
        assert_eq!(poch(-30.0, -1.5), 0.0);
        assert_eq!(poch(0.0, -1.5), 0.0);
        assert_eq!(poch(30.0, -1.5), 0.006487603131771398);
        assert_eq!(poch(1e10, -1.5), 1.0000000001875003e-15);
        assert_eq!(poch(1e100, -1.5), 1e-150);

        assert!(poch(-1e100, -100.1).is_nan());
        assert_eq!(poch(-1e10, -100.2), 0.0);
        assert_eq!(poch(-30.0, -100.3), 0.0);
        assert_eq!(poch(0.0, -100.4), 0.0);
        assert_eq!(poch(30.0, -100.5), -3.526480113279172e-132);
        assert_eq!(poch(1e10, -100.6), 0.0);
        assert_eq!(poch(1e100, -100.7), 0.0);
    }
}