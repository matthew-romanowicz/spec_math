/*  pdtri()
*
*  Inverse Poisson distribution
*
*
*
* SYNOPSIS:
*
* int k;
* double m, y, pdtr();
*
* m = pdtri( k, y );
*
*
*
*
* DESCRIPTION:
*
* Finds the Poisson variable x such that the integral
* from 0 to x of the Poisson density is equal to the
* given probability y.
*
* This is accomplished using the inverse Gamma integral
* function and the relation
*
*    m = igamci( k+1, y ).
*
*
*
*
* ACCURACY:
*
* See igami.c.
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* pdtri domain    y < 0 or y >= 1       0.0
*                     k < 0
*
*/

/*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

use crate::cephes64::igam::{igam, igamc};
use crate::cephes64::igami::igamci;

pub fn pdtrc(k: f64, m: f64) -> f64 {
    //! Complemented poisson distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the sum of the terms k+1 to infinity of the Poisson
    //! distribution:
    //!
    #![doc=include_str!("pdtrc.svg")]
    //!
    //! The terms are not summed directly; instead the incomplete
    //! Gamma integral is employed, according to the formula
    //!
    //! `y = pdtrc( k, m ) = igam( k+1, m )`
    //!
    //! The arguments must both be nonnegative.
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igam`](crate::cephes64::igam).

    if k < 0.0 || m < 0.0 {
        //sf_error("pdtrc", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if m == 0.0 {
        0.0
    } else {
        let v = k.floor() + 1.0;
        igam(v, m)
    }
}

// $$\sum_{j=0}^k{e^{-m}\,\frac{m^j}{j!}}$$

pub fn pdtr(k: f64, m: f64) -> f64 {
    //! Poisson distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the sum of the first k terms of the Poisson
    //! distribution:
    //!
    #![doc=include_str!("pdtr.svg")]
    //!
    //! The terms are not summed directly; instead the incomplete
    //! Gamma integral is employed, according to the relation
    //!
    //! `y = pdtr( k, m ) = igamc( k+1, m )`
    //!
    //! The arguments must both be nonnegative.
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igamc`](crate::cephes64::igamc).

    if k < 0.0 || m < 0.0 {
        //sf_error("pdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if m == 0.0 {
        1.0
    } else {
        let v = k.floor() + 1.0;
        igamc(v, m)
    }
}


pub fn pdtri(k: isize, y: f64) -> f64 {
    //! Inverse Poisson distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Finds the Poisson variable `x` such that the integral
    //! from `0` to `x` of the Poisson density is equal to the
    //! given probability `y`.
    //!
    //! This is accomplished using the inverse Gamma integral
    //! function and the relation
    //!
    //! `m = igamci( k+1, y )`
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igami`](crate::cephes64::igami).
    
    if k < 0 || !(0.0..1.0).contains(&y) {
        //sf_error("pdtri", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let v = (k + 1) as f64;
        igamci(v, y)
    }
}

// Very few tests due to simplicity of implementation
#[cfg(test)]
mod pdtr_tests {
    use super::*;

    #[test]
    fn pdtr_trivials() {
        assert!(pdtrc(0.0, -1e-20).is_nan());
        assert!(pdtrc(-1e-20, 0.0).is_nan());
        assert_eq!(pdtrc(0.0, 0.0), 0.0);
        assert_eq!(pdtrc(0.0, f64::INFINITY), 1.0);

        assert!(pdtr(0.0, -1e-20).is_nan());
        assert!(pdtr(-1e-20, 0.0).is_nan());
        assert_eq!(pdtr(0.0, 0.0), 1.0);
        assert_eq!(pdtr(0.0, f64::INFINITY), 0.0);

        assert!(pdtri(0, -1e-20).is_nan());
        assert!(pdtri(-1, 0.0).is_nan());
        assert!(pdtri(0, 1.0).is_nan()); // TODO: Should this be infinity?
        assert_eq!(pdtri(0, 0.0), f64::INFINITY); // TODO: Should this be 0?
    }

    #[test]
    fn pdtr_tests() {
        assert_eq!(pdtrc(1.5, 3.0), 0.8008517265285442);
        assert_eq!(pdtrc(3.0, 1.5), 0.06564245437845008);

        assert_eq!(pdtr(1.5, 3.0), 0.1991482734714558);
        assert_eq!(pdtr(3.0, 1.5), 0.9343575456215499);

        assert_eq!(pdtri(1, 0.5), 1.6783469900166612);
        assert_eq!(pdtri(3, 0.5), 3.672060748850897);
    }
}