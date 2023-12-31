/*                                                     chdtr.c
*
*     Chi-square distribution
*
*
*
* SYNOPSIS:
*
* double df, x, y, chdtr();
*
* y = chdtr( df, x );
*
*
*
* DESCRIPTION:
*
* Returns the area under the left hand tail (from 0 to x)
* of the Chi square probability density function with
* v degrees of freedom.
*
*
*                                  inf.
*                                    -
*                        1          | |  v/2-1  -t/2
*  P( x | v )   =   -----------     |   t      e     dt
*                    v/2  -       | |
*                   2    | (v/2)   -
*                                   x
*
* where x is the Chi-square variable.
*
* The incomplete Gamma integral is used, according to the
* formula
*
*     y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
*
*
* The arguments must both be positive.
*
*
*
* ACCURACY:
*
* See igam().
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* chdtr domain   x < 0 or v < 1        0.0
*/
/*							chdtrc()
*
*	Complemented Chi-square distribution
*
*
*
* SYNOPSIS:
*
* double v, x, y, chdtrc();
*
* y = chdtrc( v, x );
*
*
*
* DESCRIPTION:
*
* Returns the area under the right hand tail (from x to
* infinity) of the Chi square probability density function
* with v degrees of freedom:
*
*
*                                  inf.
*                                    -
*                        1          | |  v/2-1  -t/2
*  P( x | v )   =   -----------     |   t      e     dt
*                    v/2  -       | |
*                   2    | (v/2)   -
*                                   x
*
* where x is the Chi-square variable.
*
* The incomplete Gamma integral is used, according to the
* formula
*
*	y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).
*
*
* The arguments must both be positive.
*
*
*
* ACCURACY:
*
* See igamc().
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* chdtrc domain  x < 0 or v < 1        0.0
*/
/*							chdtri()
*
*	Inverse of complemented Chi-square distribution
*
*
*
* SYNOPSIS:
*
* double df, x, y, chdtri();
*
* x = chdtri( df, y );
*
*
*
*
* DESCRIPTION:
*
* Finds the Chi-square argument x such that the integral
* from x to infinity of the Chi-square density is equal
* to the given cumulative probability y.
*
* This is accomplished using the inverse Gamma integral
* function and the relation
*
*    x/2 = igamci( df/2, y );
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
* chdtri domain   y < 0 or y > 1        0.0
*                     v < 1
*
*/

/*                                                             chdtr() */


/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1984, 1987 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

use crate::cephes64::igam::{igam, igamc};
use crate::cephes64::igami::igamci;

pub fn chdtrc(df: f64, x: f64) -> f64 {
    //! Complemented Chi-square distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the area under the right hand tail (from x to
    //! infinity) of the Chi square probability density function
    //! with v degrees of freedom:
    //!
    #![doc=include_str!("chdtrc.svg")]
    //!
    //! where `x` is the Chi-square variable.
    //!
    //! The incomplete Gamma integral is used, according to the
    //! formula
    //!
    //! `y = chdtr( v, x ) = igamc( v/2.0, x/2.0 )`
    //!
    //! The arguments must both be positive.
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igamc`](crate::cephes64::igamc).

    if x < 0.0 {
        1.0 /* modified by T. Oliphant */
    } else {
        igamc(df / 2.0, x / 2.0)
    }	
}

// $$\mathrm{P}(x, v) = \frac{1}{2^{v/2}\,\Gamma(v/2)}\,
// \int_x^{\infty}{t^{v/2-1}\,e^{-t/2}\,dt}$$

// TODO: I think the equation in the documentation is wrong...
pub fn chdtr(df: f64, x: f64) -> f64 {
    //! Chi-square distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the area under the left hand tail (from `0` to `x`)
    //! of the Chi square probability density function with
    //! v degrees of freedom.
    //!
    #![doc=include_str!("chdtrc.svg")]
    //!
    //! where x is the Chi-square variable.
    //!
    //! The incomplete Gamma integral is used, according to the
    //! formula
    //!
    //! `y = chdtr( v, x ) = igam( v/2.0, x/2.0 )`
    //!
    //!
    //! The arguments must both be positive.
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igam`](crate::cephes64::igam).

    if x < 0.0 {		/* || (df < 1.0) ) */
        //sf_error("chdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        igam(df / 2.0, x / 2.0)
    }
}

pub fn chdtri(df: f64, y: f64) -> f64 {
    //! Inverse of complemented Chi-square distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Finds the Chi-square argument `x` such that the integral
    //! from `x` to infinity of the Chi-square density is equal
    //! to the given cumulative probability `y`.
    //!
    //! This is accomplished using the inverse Gamma integral
    //! function and the relation
    //!
    //! `x/2 = igamci( df/2, y );`
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igami`](crate::cephes64::igami).

    if !(0.0..=1.0).contains(&y) {	/* || (df < 1.0) ) */
        //sf_error("chdtri", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        2.0 * igamci(0.5 * df, y)
    }
}

// Very few tests due to simplicity of implementation
#[cfg(test)]
mod chdtr_tests {
    use super::*;

    #[test]
    fn chdtr_trivials() {
        assert_eq!(chdtrc(10.0, -1e-10), 1.0);
        assert_eq!(chdtrc(10.0, 0.0), 1.0);
        assert_eq!(chdtrc(10.0, f64::INFINITY), 0.0);

        assert_eq!(chdtr(10.0, -1e-10).is_nan(), true);
        assert_eq!(chdtr(10.0, 0.0), 0.0);
        assert_eq!(chdtr(10.0, f64::INFINITY), 1.0);

        assert_eq!(chdtri(10.0, -1e-10).is_nan(), true);
        assert_eq!(chdtri(10.0, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(chdtri(10.0, 0.0), f64::INFINITY);
        assert_eq!(chdtri(10.0, 1.0), 0.0);
    }

    #[test]
    fn chdtr_tests() {
        assert_eq!(chdtrc(1.5, 3.0), 0.1475995543647508);
        assert_eq!(chdtrc(3.0, 1.5), 0.6822703303362125);

        assert_eq!(chdtr(1.5, 3.0), 0.8524004456352492);
        assert_eq!(chdtr(3.0, 1.5), 0.31772966966378746);

        assert_eq!(chdtri(1.5, 0.5), 0.9083339566477122);
        assert_eq!(chdtri(3.0, 0.5), 2.3659738843753377);
    }
}