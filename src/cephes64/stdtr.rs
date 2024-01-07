/*                                                     stdtr.c
*
*     Student's t distribution
*
*
*
* SYNOPSIS:
*
* double t, stdtr();
* short k;
*
* y = stdtr( k, t );
*
*
* DESCRIPTION:
*
* Computes the integral from minus infinity to t of the Student
* t distribution with integer k > 0 degrees of freedom:
*
*                                      t
*                                      -
*                                     | |
*              -                      |         2   -(k+1)/2
*             | ( (k+1)/2 )           |  (     x   )
*       ----------------------        |  ( 1 + --- )        dx
*                     -               |  (      k  )
*       sqrt( k pi ) | ( k/2 )        |
*                                   | |
*                                    -
*                                   -inf.
*
* Relation to incomplete beta integral:
*
*        1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
* where
*        z = k/(k + t**2).
*
* For t < -2, this is the method of computation.  For higher t,
* a direct method is derived from integration by parts.
* Since the function is symmetric about t=0, the area under the
* right tail of the density is found by calling the function
* with -t instead of t.
*
* ACCURACY:
*
* Tested at random 1 <= k <= 25.  The "domain" refers to t.
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE     -100,-2      50000       5.9e-15     1.4e-15
*    IEEE     -2,100      500000       2.7e-15     4.9e-17
*/

/*                                                     stdtri.c
*
*     Functional inverse of Student's t distribution
*
*
*
* SYNOPSIS:
*
* double p, t, stdtri();
* int k;
*
* t = stdtri( k, p );
*
*
* DESCRIPTION:
*
* Given probability p, finds the argument t such that stdtr(k,t)
* is equal to p.
*
* ACCURACY:
*
* Tested at random 1 <= k <= 100.  The "domain" refers to p:
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE    .001,.999     25000       5.7e-15     8.0e-16
*    IEEE    10^-6,.001    25000       2.0e-12     2.9e-14
*/


/*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

use crate::cephes64::consts::{M_PI, MACHEP};
use crate::cephes64::incbet::incbet;
use crate::cephes64::incbi::incbi;

pub fn stdtr(k: isize, t: f64) -> f64
{
    //double x, rk, z, f, tz, p, xsqk;
    //int j;

    if (k <= 0) {
        //sf_error("stdtr", SF_ERROR_DOMAIN, NULL);
        return (f64::NAN);
    }

    if (t == 0.0) {
        return (0.5);
    }

    if (t < -2.0) {
        let rk = k as f64;
        let z = rk / (rk + t * t);
        let p = 0.5 * incbet(0.5 * rk, 0.5, z);
        return (p);
    }

    /*     compute integral from -t to + t */

    let x = if (t < 0.0) {
        -t
    } else {
        t
    };

    let rk = k as f64;			/* degrees of freedom */
    let z = 1.0 + (x * x) / rk;

    let mut p: f64;
    /* test if k is odd or even */
    if ((k & 1) != 0) {

        /*      computation for odd k   */

        let xsqk = x / rk.sqrt();
        p = xsqk.atan();
        if (k > 1) {
            let mut f = 1.0;
            let mut tz = 1.0;
            let mut j = 3;
            while ((j <= k - 2) && ((tz / f) > MACHEP)) {
                tz *= (j - 1) as f64 / (z * j as f64);
                f += tz;
                j += 2;
            }
            p += f * xsqk / z;
        }
        p *= 2.0 / M_PI;
    }


    else {

        /*      computation for even k  */

        let mut f = 1.0;
        let mut tz = 1.0;
        let mut j = 2;

        while ((j <= (k - 2)) && ((tz / f) > MACHEP)) {
            tz *= (j - 1) as f64 / (z * j as f64);
            f += tz;
            j += 2;
        }
        p = f * x / (z * rk).sqrt();
    }

    /*     common exit     */


    if (t < 0.0) {
        p = -p;			/* note destruction of relative accuracy */
    }

    p = 0.5 + 0.5 * p;
    return (p);
}

pub fn stdtri(k: isize, p: f64) -> f64
{
    // double t, rk, z;
    // int rflg;

    if (k <= 0 || p <= 0.0 || p >= 1.0) {
        //sf_error("stdtri", SF_ERROR_DOMAIN, NULL);
        return (f64::NAN);
    }

    let rk = k as f64;

    if (p > 0.25 && p < 0.75) {
        if (p == 0.5) {
            return (0.0);
        }
        let z = 1.0 - 2.0 * p;
        let z = incbi(0.5, 0.5 * rk, z.abs());
        let mut t = (rk * z / (1.0 - z)).sqrt();
        if (p < 0.5) {
            t = -t;
        }
        return (t);
    }
    
    let mut rflg = -1;
    let mut p = p;
    if (p >= 0.5) {
        p = 1.0 - p;
        rflg = 1;
    }
    let z = incbi(0.5 * rk, 0.5, 2.0 * p);

    if (f64::MAX * z < rk) {
        return (rflg as f64 * f64::INFINITY);
    } 
    let t = (rk / z - rk).sqrt();
    return (rflg as f64 * t);
}