use crate::cephes64::{M_2_SQRTPI, M_SQRT1_2, ndtri};

///Inverse of the error function.
///
/// Computes the inverse of the error function on the restricted domain
/// -1 < y < 1. This restriction ensures the existence of a unique result
/// such that erf(erfinv(y)) = y.
pub fn erfinv(y: f64) -> f64 {
    /* 
    * For small arguments, use the Taylor expansion
    * erf(y) = 2/\sqrt{\pi} (y - y^3 / 3 + O(y^5)),    y\to 0 
    * where we only retain the linear term.
    * Otherwise, y + 1 loses precision for |y| << 1.
    */
    if -1e-7 < y && y < 1e-7 {
        y / M_2_SQRTPI
    } else if -1.0 < y && y < 1.0 {
        ndtri(0.5 * (y + 1.0)) * M_SQRT1_2
    } else if y == -1.0 {
        -f64::INFINITY
    } else if y == 1.0 {
        f64::INFINITY
    }
    else if y.is_nan() {
        //sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        y
    } else {
        //sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    }
}


/// Inverse of the complementary error function.
///
/// Computes the inverse of the complimentary error function on the restricted
/// domain 0 < y < 2. This restriction ensures the existence of a unique result
/// such that erfc(erfcinv(y)) = y.
pub fn erfcinv(y: f64) -> f64 {
    
    if 0.0 < y && y < 2.0 {
        -ndtri(0.5 * y) * M_SQRT1_2
    } else if y == 0.0 {
        f64::INFINITY
    } else if y == 2.0 {
        -f64::INFINITY
    }
    else if y.is_nan() {
        //sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        y
    }
    else {
        //sf_error("erfcinv", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    }
}

#[cfg(test)]
mod erfinv_tests {
    use super::*;

    #[test]
    fn erfinv_trivials() {
        assert_eq!(erfinv(1.0), f64::INFINITY);
        assert_eq!(erfinv(-1.0), -f64::INFINITY);
        assert_eq!(erfinv(f64::NAN).is_nan(), true);
        assert_eq!(erfinv(1.0 + 1e-10).is_nan(), true);
        assert_eq!(erfinv(-1.0 - 1e-10).is_nan(), true);
    }

    #[test]
    fn erfinv_small() {
        assert_eq!(erfinv(0.0), 0.0);

        assert_eq!(erfinv(5e-8), 4.43113462726379e-8);
        assert_eq!(erfinv(1e-7 - 1e-20), 8.862269254526693e-8);

        assert_eq!(erfinv(-5e-8), -4.43113462726379e-8);
        assert_eq!(erfinv(-1e-7 + 1e-20), -8.862269254526693e-8);
    }

    #[test]
    fn erfinv_large() {
        // TODO: This is inaccurate for small x
        assert_eq!(erfinv(1e-7), 8.862269259701993e-8);
        assert_eq!(erfinv(1e-5), 8.862269254817653e-06);
        assert_eq!(erfinv(1e-3), 0.0008862271574664545);
        assert_eq!(erfinv(0.1), 0.08885599049425778);
        assert_eq!(erfinv(0.5), 0.4769362762044699);
        assert_eq!(erfinv(0.9), 1.1630871536766738);
        assert_eq!(erfinv(1.0 - 1e-3), 2.3267537655135464);
        assert_eq!(erfinv(1.0 - 1e-5), 3.1234132743398733);
        assert_eq!(erfinv(1.0 - 1e-10), 4.572824958544925);
        assert_eq!(erfinv(1.0 - 1e-15), 5.68612844131039);

        // TODO: This is inaccurate, especially as x approaches 1
        assert_eq!(erfinv(-1e-7), -8.862269249862897e-8);
        assert_eq!(erfinv(-1e-5), -8.862269254719263e-06);
        assert_eq!(erfinv(-1e-3), -0.000886227157466553);
        assert_eq!(erfinv(-0.1), -0.08885599049425769);
        assert_eq!(erfinv(-0.5), -0.4769362762044699);
        assert_eq!(erfinv(-0.9), -1.1630871536766743);
        assert_eq!(erfinv(-1.0 + 1e-3), -2.3267537655135246);
        assert_eq!(erfinv(-1.0 + 1e-5), -3.123413274341571);
        assert_eq!(erfinv(-1.0 + 1e-10), -4.572824958544925);
        assert_eq!(erfinv(-1.0 + 1e-15), -5.675915739744713);

    }
}

#[cfg(test)]
mod erfcinv_tests {
    use super::*;

    #[test]
    fn erfcinv_trivials() {
        assert_eq!(erfcinv(0.0), f64::INFINITY);
        assert_eq!(erfcinv(2.0), -f64::INFINITY);
        assert_eq!(erfcinv(f64::NAN).is_nan(), true);
        assert_eq!(erfcinv(2.0 + 1e-10).is_nan(), true);
        assert_eq!(erfcinv(-1e-10).is_nan(), true);
    }

    #[test]
    fn erfcinv_values() {
        assert_eq!(erfcinv(1e-300), 26.209469960516124);
        assert_eq!(erfcinv(1e-200), 21.374783049026263);
        assert_eq!(erfcinv(1e-100), 15.065574702592647);
        assert_eq!(erfcinv(1e-50), 10.592090169527365);
        assert_eq!(erfcinv(1e-20), 6.601580622355143);
        assert_eq!(erfcinv(1e-10), 4.572824967389486);
        assert_eq!(erfcinv(1e-5), 3.1234132743408756);
        assert_eq!(erfcinv(0.1), 1.1630871536766743);
        assert_eq!(erfcinv(0.5), 0.4769362762044699);
        assert_eq!(erfcinv(1.0), 0.0);
        assert_eq!(erfcinv(1.5), -0.4769362762044699);
        assert_eq!(erfcinv(1.9), -1.1630871536766738);
        assert_eq!(erfcinv(2.0 - 1e-5), -3.1234132743398733);
        assert_eq!(erfcinv(2.0 - 1e-10), -4.572824958544925);
        assert_eq!(erfcinv(2.0 - 1e-15), -5.666765021916333);
    }
}