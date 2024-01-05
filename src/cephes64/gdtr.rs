/*                                                     gdtr.c
*
*     Gamma distribution function
*
*
*
* SYNOPSIS:
*
* double a, b, x, y, gdtr();
*
* y = gdtr( a, b, x );
*
*
*
* DESCRIPTION:
*
* Returns the integral from zero to x of the Gamma probability
* density function:
*
*
*                x
*        b       -
*       a       | |   b-1  -at
* y =  -----    |    t    e    dt
*       -     | |
*      | (b)   -
*               0
*
*  The incomplete Gamma integral is used, according to the
* relation
*
* y = igam( b, ax ).
*
*
* ACCURACY:
*
* See igam().
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* gdtr domain         x < 0            0.0
*
*/
/*							gdtrc.c
*
*	Complemented Gamma distribution function
*
*
*
* SYNOPSIS:
*
* double a, b, x, y, gdtrc();
*
* y = gdtrc( a, b, x );
*
*
*
* DESCRIPTION:
*
* Returns the integral from x to infinity of the Gamma
* probability density function:
*
*
*               inf.
*        b       -
*       a       | |   b-1  -at
* y =  -----    |    t    e    dt
*       -     | |
*      | (b)   -
*               x
*
*  The incomplete Gamma integral is used, according to the
* relation
*
* y = igamc( b, ax ).
*
*
* ACCURACY:
*
* See igamc().
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* gdtrc domain         x < 0            0.0
*
*/

/*                                                     gdtr()  */


/*
* Cephes Math Library Release 2.3:  March,1995
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

use crate::cephes64::igam::{igam, igamc};
use crate::cephes64::igami::igamci;

// $$y = \frac{a^b}{\Gamma(b)}\,\int_0^x{t^{b-1}\,e^{-a\,t}\,dt}$$

pub fn gdtr(a: f64, b: f64, x: f64) -> f64 {
    //! Gamma distribution function
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the integral from zero to x of the Gamma probability
    //! density function:
    //!
    #![doc=include_str!("gdtr.svg")]
    //!
    //! The incomplete Gamma integral is used, according to the
    //! relation
    //!
    //! `y = igam( b, ax )`
    //!
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igam`](crate::cephes64::igam).

    if x < 0.0 {
        //sf_error("gdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        igam(b, a * x)
    }
}

// $$y = \frac{a^b}{\Gamma(b)}\,\int_x^{\infty}{t^{b-1}\,e^{-a\,t}\,dt}$$

pub fn gdtrc(a: f64, b: f64, x: f64) -> f64 {
    //! Complemented Gamma distribution function
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the integral from x to infinity of the Gamma
    //! probability density function:
    //!
    #![doc=include_str!("gdtrc.svg")]
    //!
    //! The incomplete Gamma integral is used, according to the
    //! relation
    //!
    //! `y = igamc( b, ax )`
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::igamc`](crate::cephes64::igamc).

    if x < 0.0 {
        //sf_error("gdtrc", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        igamc(b, a * x)
    }
}


pub fn gdtri(a: f64, b: f64, y: f64) -> f64 {
    //! Inverse gamma distribution function

    if y < 0.0 || y > 1.0 || a <= 0.0 || b < 0.0 {
        //sf_error("gdtri", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        igamci(b, 1.0 - y) / a
    }
}

// Very few tests due to simplicity of implementation
#[cfg(test)]
mod gdtr_tests {
    use super::*;

    #[test]
    fn gdtr_trivials() {
        assert_eq!(gdtrc(0.0, 0.0, 0.0).is_nan(), true);
        assert_eq!(gdtrc(0.1, 0.0, 0.0).is_nan(), true);

        assert_eq!(gdtr(0.0, 1.0, 1.0), 0.0);
        assert_eq!(gdtr(1.0, 0.0, 1.0), 1.0); // TODO: Is this right?
        assert_eq!(gdtr(1.0, 1.0, -1e-10).is_nan(), true);

        assert_eq!(gdtri(0.0, 1.0, 1.0).is_nan(), true);
        assert_eq!(gdtri(1.0, 0.0, 1.0), f64::INFINITY);
        assert_eq!(gdtri(1.0, 1.0, 0.0), 0.0);
        assert_eq!(gdtri(1.0, 1.0, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(gdtri(1.0, 1.0, 1.0), f64::INFINITY);
    }

    #[test]
    fn gdtr_tests() {
        assert_eq!(gdtrc(1.0, 1.0, 1.0), 0.36787944117144245);
        assert_eq!(gdtrc(1.0, 1.0, 0.0), 1.0);
        assert_eq!(gdtrc(2.0, 3.0, 4.0), 0.013753967744002971);

        assert_eq!(gdtr(1.0, 1.0, 1.0), 0.6321205588285577);
        assert_eq!(gdtr(1.0, 1.0, 0.0), 0.0);
        assert_eq!(gdtr(2.0, 3.0, 4.0), 0.986246032255997);

        assert_eq!(gdtri(1.0, 1.0, 0.3), 0.35667494393873245);
        assert_eq!(gdtri(1.0, 1.0, 0.6), 0.9162907318741551);
        assert_eq!(gdtri(2.0, 3.0, 0.8), 2.139514930062667);
        assert_eq!(gdtri(2.0, 3.0, 0.9999), 6.964085309003543);
    }
}