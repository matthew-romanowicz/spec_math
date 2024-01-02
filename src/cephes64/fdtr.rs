/*                                                     fdtr.c
*
*     F distribution
*
*
*
* SYNOPSIS:
*
* double df1, df2;
* double x, y, fdtr();
*
* y = fdtr( df1, df2, x );
*
* DESCRIPTION:
*
* Returns the area from zero to x under the F density
* function (also known as Snedcor's density or the
* variance ratio density).  This is the density
* of x = (u1/df1)/(u2/df2), where u1 and u2 are random
* variables having Chi square distributions with df1
* and df2 degrees of freedom, respectively.
*
* The incomplete beta integral is used, according to the
* formula
*
*     P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).
*
*
* The arguments a and b are greater than zero, and x is
* nonnegative.
*
* ACCURACY:
*
* Tested at random points (a,b,x).
*
*                x     a,b                     Relative error:
* arithmetic  domain  domain     # trials      peak         rms
*    IEEE      0,1    0,100       100000      9.8e-15     1.7e-15
*    IEEE      1,5    0,100       100000      6.5e-15     3.5e-16
*    IEEE      0,1    1,10000     100000      2.2e-11     3.3e-12
*    IEEE      1,5    1,10000     100000      1.1e-11     1.7e-13
* See also incbet.c.
*
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* fdtr domain     a<0, b<0, x<0         0.0
*
*/

/*                         fdtrc()
*
*  Complemented F distribution
*
*
*
* SYNOPSIS:
*
* double df1, df2;
* double x, y, fdtrc();
*
* y = fdtrc( df1, df2, x );
*
* DESCRIPTION:
*
* Returns the area from x to infinity under the F density
* function (also known as Snedcor's density or the
* variance ratio density).
*
*
*                      inf.
*                       -
*              1       | |  a-1      b-1
* 1-P(x)  =  ------    |   t    (1-t)    dt
*            B(a,b)  | |
*                     -
*                      x
*
*
* The incomplete beta integral is used, according to the
* formula
*
*  P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).
*
*
* ACCURACY:
*
* Tested at random points (a,b,x) in the indicated intervals.
*                x     a,b                     Relative error:
* arithmetic  domain  domain     # trials      peak         rms
*    IEEE      0,1    1,100       100000      3.7e-14     5.9e-16
*    IEEE      1,5    1,100       100000      8.0e-15     1.6e-15
*    IEEE      0,1    1,10000     100000      1.8e-11     3.5e-13
*    IEEE      1,5    1,10000     100000      2.0e-11     3.0e-12
* See also incbet.c.
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* fdtrc domain    a<0, b<0, x<0         0.0
*
*/

/*                         fdtri()
*
*  Inverse of F distribution
*
*
*
* SYNOPSIS:
*
* double df1, df2;
* double x, p, fdtri();
*
* x = fdtri( df1, df2, p );
*
* DESCRIPTION:
*
* Finds the F density argument x such that the integral
* from -infinity to x of the F density is equal to the
* given probability p.
*
* This is accomplished using the inverse beta integral
* function and the relations
*
*      z = incbi( df2/2, df1/2, p )
*      x = df2 (1-z) / (df1 z).
*
* Note: the following relations hold for the inverse of
* the uncomplemented F distribution:
*
*      z = incbi( df1/2, df2/2, p )
*      x = df2 z / (df1 (1-z)).
*
* ACCURACY:
*
* Tested at random points (a,b,p).
*
*              a,b                     Relative error:
* arithmetic  domain     # trials      peak         rms
*  For p between .001 and 1:
*    IEEE     1,100       100000      8.3e-15     4.7e-16
*    IEEE     1,10000     100000      2.1e-11     1.4e-13
*  For p between 10^-6 and 10^-3:
*    IEEE     1,100        50000      1.3e-12     8.4e-15
*    IEEE     1,10000      50000      3.0e-12     4.8e-14
* See also fdtrc.c.
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* fdtri domain   p <= 0 or p > 1       NaN
*                     v < 1
*
*/

/*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/


use crate::cephes64::incbet::incbet;
use crate::cephes64::incbi::incbi;

pub fn fdtrc(a: f64, b: f64, x: f64) -> f64 {
    //! Complemented F distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the area from `x` to infinity under the F density
    //! function (also known as Snedcor's density or the
    //! variance ratio density).
    //!
    #![doc=include_str!("fdtrc.svg")]
    //!
    //! The incomplete beta integral is used, according to the
    //! formula
    //!
    //! `P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).`
    //!
    //! ## ACCURACY:
    //!
    //! Tested at random points (`a`, `b`, `x`) in the indicated intervals.
    //!
    //! Relative Error:
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>x Domain</th>
    //!     <th>a, b Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 1</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>3.7e-14</td>
    //!     <td> 5.9e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>1, 5</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>8.0e-15</td>
    //!     <td>1.6e-15</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 1</td>
    //!     <td>0, 10000</td>
    //!     <td>100000</td>
    //!     <td>1.8e-11</td>
    //!     <td>3.5e-13</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>1, 5</td>
    //!     <td>0, 10000</td>
    //!     <td>100000</td>
    //!     <td>2.0e-11</td>
    //!     <td>3.0e-12</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::incbet`](crate::cephes64::incbet).

    if a <= 0.0 || b <= 0.0 || x < 0.0 {
        //sf_error("fdtrc", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let w = b / (b + a * x);
        incbet(0.5 * b, 0.5 * a, w)
    }
}


pub fn fdtr(a: f64, b: f64, x: f64) -> f64 {
    //! F distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the area from zero to x under the F density
    //! function (also known as Snedcor's density or the
    //! variance ratio density).  This is the density
    //! of x = (u1/df1)/(u2/df2), where u1 and u2 are random
    //! variables having Chi square distributions with df1
    //! and df2 degrees of freedom, respectively.
    //!
    //! The incomplete beta integral is used, according to the
    //! formula
    //!
    //! `P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) )`
    //!
    //! The arguments `a` and `b` are greater than zero, and `x` is
    //! nonnegative.
    //!
    //! ## ACCURACY:
    //!
    //! Tested at random points (`a`, `b`, `x`).
    //!
    //! Relative Error:
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>x Domain</th>
    //!     <th>a, b Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 1</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>9.8e-15</td>
    //!     <td>1.7e-15</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>1, 5</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>6.5e-15</td>
    //!     <td>3.5e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 1</td>
    //!     <td>0, 10000</td>
    //!     <td>100000</td>
    //!     <td>2.2e-11</td>
    //!     <td>3.3e-12</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>1, 5</td>
    //!     <td>0, 10000</td>
    //!     <td>100000</td>
    //!     <td>1.1e-11</td>
    //!     <td>1.7e-13</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::incbet`](crate::cephes64::incbet).

    if a <= 0.0 || b <= 0.0 || x < 0.0 {
        //sf_error("fdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let mut w = a * x;
        w = w / (b + w);
        incbet(0.5 * a, 0.5 * b, w)
    }
}


pub fn fdtri(a: f64, b: f64, y: f64) -> f64 {
    //! Inverse of F distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Finds the F density argument x such that the integral
    //! from -infinity to x of the F density is equal to the
    //! given probability p.
    //!
    //! This is accomplished using the inverse beta integral
    //! function and the relations
    //!
    //! `z = incbi( df2/2, df1/2, p )`
    //!
    //! `x = df2 (1-z) / (df1 z)`
    //!
    //! Note: the following relations hold for the inverse of
    //! the uncomplemented F distribution:
    //!
    //! `z = incbi( df1/2, df2/2, p )`
    //!
    //! `x = df2 z / (df1 (1-z))`
    //!
    //! ## ACCURACY:
    //!
    //! Tested at random points (`a`, `b`, `p`).
    //!
    //! Relative error:
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>p Domain</th>
    //!     <th>a, b Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0.001, 1</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>8.3e-15</td>
    //!     <td>4.7e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0.001, 1</td>
    //!     <td>0, 10000</td>
    //!     <td>100000</td>
    //!     <td>2.1e-11</td>
    //!     <td>1.4e-13</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>1e-6, 1e-3</td>
    //!     <td>0, 100</td>
    //!     <td>50000</td>
    //!     <td>1.3e-12</td>
    //!     <td>8.4e-15</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>1e-6, 1e-3</td>
    //!     <td>0, 10000</td>
    //!     <td>50000</td>
    //!     <td>3.0e-12</td>
    //!     <td>4.8e-14</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::fdtrc`](crate::cephes64::fdtrc).

    if a <= 0.0 || b <= 0.0 || y <= 0.0 || y > 1.0 {
        //sf_error("fdtri", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    } 

    let y = 1.0 - y;
    /* Compute probability for x = 0.5.  */
    let mut w = incbet(0.5 * b, 0.5 * a, 0.5);
    /* If that is greater than y, then the solution w < .5.
    * Otherwise, solve at 1-y to remove cancellation in (b - b*w).  */
    if w > y || y < 0.001 {
        w = incbi(0.5 * b, 0.5 * a, y);
        (b - b * w) / (a * w)
    } else {
        w = incbi(0.5 * a, 0.5 * b, 1.0 - y);
        b * w / (a * (1.0 - w))
    }
}

// Very few tests due to simplicity of implementation
#[cfg(test)]
mod fdtr_tests {
    use super::*;

    #[test]
    fn fdtr_trivials() {
        assert_eq!(fdtrc(0.0, 1.0, 1.0).is_nan(), true);
        assert_eq!(fdtrc(1.0, 0.0, 1.0).is_nan(), true);
        assert_eq!(fdtrc(1.0, 1.0, -1e-10).is_nan(), true);

        assert_eq!(fdtr(0.0, 1.0, 1.0).is_nan(), true);
        assert_eq!(fdtr(1.0, 0.0, 1.0).is_nan(), true);
        assert_eq!(fdtr(1.0, 1.0, -1e-10).is_nan(), true);

        assert_eq!(fdtri(0.0, 1.0, 1.0).is_nan(), true);
        assert_eq!(fdtri(1.0, 0.0, 1.0).is_nan(), true);
        assert_eq!(fdtri(1.0, 1.0, 0.0).is_nan(), true);
        assert_eq!(fdtri(1.0, 1.0, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(fdtri(1.0, 1.0, 1.0), f64::INFINITY);
    }

    #[test]
    fn fdtr_tests() {
        assert_eq!(fdtrc(1.0, 1.0, 1.0), 0.5000000000000001);
        assert_eq!(fdtrc(1.0, 1.0, 0.0), 1.0);
        assert_eq!(fdtrc(2.0, 3.0, 4.0), 0.14242717305466185);

        assert_eq!(fdtr(1.0, 1.0, 1.0), 0.5000000000000001);
        assert_eq!(fdtr(1.0, 1.0, 0.0), 0.0);
        assert_eq!(fdtr(2.0, 3.0, 4.0), 0.8575728269453382);

        assert_eq!(fdtri(1.0, 1.0, 0.3), 0.2596161836824997);
        assert_eq!(fdtri(1.0, 1.0, 0.6), 1.8944271909999169);
        assert_eq!(fdtri(2.0, 3.0, 0.8), 2.8860266073192995);
        assert_eq!(fdtri(2.0, 3.0, 0.9999), 694.738325041968);
    }
}