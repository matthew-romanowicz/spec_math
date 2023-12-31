/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

use crate::cephes64::incbet::incbet;

pub fn btdtr(a: f64, b: f64, x: f64) -> f64 {
    //! Beta distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the area from zero to x under the beta density
    //! function.
    //!
    //! This function is identical to the incomplete beta
    //! integral function `incbet(a, b, x)`.
    //!
    //! The complemented function is
    //!
    //! `1 - P(1-x)  =  incbet( b, a, x );`
    //!
    //! ## ACCURACY:
    //!
    //! See [`cephes64::incbet`](crate::cephes64::incbet)

    incbet(a, b, x)
}