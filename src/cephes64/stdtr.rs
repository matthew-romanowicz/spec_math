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

    if k <= 0 {
        //sf_error("stdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if t == 0.0 {
        0.5
    } else if t < -2.0 {
        let rk = k as f64;
        let z = rk / (rk + t * t);
        0.5 * incbet(0.5 * rk, 0.5, z)
    } else {
        // compute integral from -t to + t

        let x = t.abs();

        let rk = k as f64; // degrees of freedom
        let z = 1.0 + (x * x) / rk;

        let mut p: f64;

        // test if k is odd or even
        if k & 1 != 0 {
            // computation for odd k
            let xsqk = x / rk.sqrt();
            p = xsqk.atan();
            if k > 1 {
                let mut f = 1.0;
                let mut tz = 1.0;
                let mut j = 3;
                while j <= k - 2 && tz / f > MACHEP {
                    tz *= (j - 1) as f64 / (z * j as f64);
                    f += tz;
                    j += 2;
                }
                p += f * xsqk / z;
            }
            p *= 2.0 / M_PI;

        } else {
            // computation for even k

            let mut f = 1.0;
            let mut tz = 1.0;
            let mut j = 2;

            while j <= k - 2 && tz / f > MACHEP {
                tz *= (j - 1) as f64 / (z * j as f64);
                f += tz;
                j += 2;
            }
            p = f * x / (z * rk).sqrt();
        }

        // common exit


        if t < 0.0 {
            p = -p; // note destruction of relative accuracy
        }

        0.5 + 0.5 * p
    }
}

pub fn stdtri(k: isize, p: f64) -> f64
{
    // double t, rk, z;
    // int rflg;

    if k <= 0 || !(0.0..=1.0).contains(&p) {
        //sf_error("stdtri", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else if p == 0.0 { //M. Romanowicz: p = 0 and 1 should return inf, not nan
        -f64::INFINITY
    } else if p == 1.0 {
        f64::INFINITY
    } else {

        let rk = k as f64;

        if p > 0.25 && p < 0.75 {
            if p == 0.5 {
                return 0.0;
            }
            let z = 1.0 - 2.0 * p;
            let z = incbi(0.5, 0.5 * rk, z.abs());
            let t = (rk * z / (1.0 - z)).sqrt();
            if p < 0.5 {
                -t
            } else {
                t
            }
        } else {
            let mut rflg = -1;
            let mut p = p;
            if p >= 0.5 {
                p = 1.0 - p;
                rflg = 1;
            }
            let z = incbi(0.5 * rk, 0.5, 2.0 * p);

            if f64::MAX * z < rk {
                return rflg as f64 * f64::INFINITY;
            } 
            let t = (rk / z - rk).sqrt();
            rflg as f64 * t
        }
    }
}

#[cfg(test)]
mod stdtr_tests {
    use super::*;

    #[test]
    fn stdtr_trivials() {
        assert!(stdtr(0, 0.0).is_nan());
        assert!(stdtr(-1, 0.5).is_nan());
        assert_eq!(stdtr(1, 0.0), 0.5);
        assert_eq!(stdtr(1, f64::INFINITY), 1.0);
        assert_eq!(stdtr(1, -f64::INFINITY), 0.0);
    }

    #[test]
    fn stdtr_odd_k() {
        assert_eq!(stdtr(1, -2.0), 0.14758361765043326);
        assert_eq!(stdtr(1, -1e-10), 0.499999999968169);
        assert_eq!(stdtr(1, 1e-10), 0.500000000031831);
        assert_eq!(stdtr(1, 1.5), 0.8128329581890013);
        assert_eq!(stdtr(1, 10.0), 0.9682744825694465);
        assert_eq!(stdtr(1, 1e5), 0.9999968169011383);
        assert_eq!(stdtr(1, 1e10), 0.999999999968169);

        assert_eq!(stdtr(5, -2.0), 0.05096973941492916);
        assert_eq!(stdtr(5, -1e-10), 0.4999999999620393);
        assert_eq!(stdtr(5, 1e-10), 0.5000000000379606);
        assert_eq!(stdtr(5, 1.5), 0.9030481598787634);
        assert_eq!(stdtr(5, 10.0), 0.9999145262121285);
        assert_eq!(stdtr(5, 5e2), 0.9999999999996964);

        assert_eq!(stdtr(101, -2.0), 0.02409260560369547);
        assert_eq!(stdtr(101, -1e-10), 0.4999999999602044);
        assert_eq!(stdtr(101, 1e-10), 0.5000000000397956);
        assert_eq!(stdtr(101, 1.5), 0.9316330369201067);
        assert_eq!(stdtr(101, 9.0), 0.9999999999999928);

        // TODO: This takes forever and is inaccurate
        //assert_eq!(stdtr(2147483647, -2.0), 0.022750132011032955);
    }

    #[test]
    fn stdtr_even_k() {
        assert_eq!(stdtr(2, -2.0), 0.09175170953613693);
        assert_eq!(stdtr(2, -1e-10), 0.49999999996464467);
        assert_eq!(stdtr(2, 1e-10), 0.5000000000353554);
        assert_eq!(stdtr(2, 1.5), 0.8638034375544994);
        assert_eq!(stdtr(2, 10.0), 0.9950737714883372);
        assert_eq!(stdtr(2, 1e5), 0.99999999995);
        assert_eq!(stdtr(2, 1e6), 0.9999999999995);

        assert_eq!(stdtr(6, -2.0), 0.046213155765837566);
        assert_eq!(stdtr(6, -1e-10), 0.4999999999617267);
        assert_eq!(stdtr(6, 1e-10), 0.5000000000382733);
        assert_eq!(stdtr(6, 1.5), 0.9078596319292591);
        assert_eq!(stdtr(6, 10.0), 0.9999710400862252);
        assert_eq!(stdtr(6, 5e2), 0.9999999999999979);

        assert_eq!(stdtr(100, -2.0), 0.024106089365566685);
        assert_eq!(stdtr(100, -1e-10), 0.4999999999602054);
        assert_eq!(stdtr(100, 1e-10), 0.5000000000397946);
        assert_eq!(stdtr(100, 1.5), 0.9316174709376557);
        assert_eq!(stdtr(100, 9.0), 0.9999999999999927);

    }

    #[test]
    fn stdtr_t_lt_neg2() {
        assert_eq!(stdtr(1, -1e150), 3.1830988618379074e-151);
        assert_eq!(stdtr(1, -10.0), 0.031725517430553574);
        assert_eq!(stdtr(1, -2.0 - 1e-10), 0.1475836176440671);

        assert_eq!(stdtr(10, -1e20), 1.2304687499999997e-196);
        assert_eq!(stdtr(10, -10.0), 7.947765877982062e-7);
        assert_eq!(stdtr(10, -2.0 - 1e-10), 0.03669401737925559);

        assert_eq!(stdtr(100, -20.0), 4.997133930668492e-37);
        assert_eq!(stdtr(100, -10.0), 4.95084449229707e-17);
        assert_eq!(stdtr(100, -2.0 - 1e-10), 0.024106089360076077);

        // TODO: This is inaccurate
        // assert_eq!(stdtr(2147483647, -20.0), 2.753675665465448e-89);
        // assert_eq!(stdtr(2147483647, -10.0), 7.61986207143429e-24);
        // assert_eq!(stdtr(2147483647, -2.0 - 1e-10), 0.022750132005633864);
        // assert_eq!(stdtr(2147483647, -2.0), 0.022750132011032955);
    }
}

#[cfg(test)]
mod stdtri_tests {
    use super::*;

    #[test]
    fn stdtri_trivials() {
        assert!(stdtri(0, 0.5).is_nan());
        assert!(stdtri(-1, 0.5).is_nan());
        assert!(stdtri(1, -1e-10).is_nan());
        assert!(stdtri(1, 1.0 + 1e10).is_nan());
        assert_eq!(stdtri(1, 0.0), -f64::INFINITY);
        assert_eq!(stdtri(1, 1.0), f64::INFINITY);
        assert_eq!(stdtri(100, 0.5), 0.0);
    }

    #[test]
    fn stdtri_p_medium() { // 0.25 < p < 0.75

        // Note: These results do not match scipy's stdtrit, but they appear
        // to be more accurate
        assert_eq!(stdtri(1, 0.25 + 1e-10), -0.9999999993716808);
        assert_eq!(stdtri(1, 0.35), -0.5095254494944288);
        assert_eq!(stdtri(1, 0.5 - 1e-10), -3.141592913526334e-10);
        assert_eq!(stdtri(1, 0.5 + 1e-10), 3.141592913526334e-10);
        assert_eq!(stdtri(1, 0.65), 0.5095254494944288);
        assert_eq!(stdtri(1, 0.75 - 1e-10), 0.9999999993716808);

        assert_eq!(stdtri(10, 0.25 + 1e-10), -0.6998120609781328);
        assert_eq!(stdtri(10, 0.35), -0.39659149375562175);
        assert_eq!(stdtri(10, 0.5 - 1e-10), -2.5699782475714284e-10);
        assert_eq!(stdtri(10, 0.5 + 1e-10), 2.5699782475714284e-10);
        assert_eq!(stdtri(10, 0.65), 0.39659149375562175);
        assert_eq!(stdtri(10, 0.75 - 1e-10), 0.6998120609781328);

        assert_eq!(stdtri(100, 0.25 + 1e-10), -0.6769510426949146);
        assert_eq!(stdtri(100, 0.35), -0.38642898040767143);
        assert_eq!(stdtri(100, 0.5 - 1e-10), -2.5129027882894715e-10);
        assert_eq!(stdtri(100, 0.5 + 1e-10), 2.5129027882894715e-10);
        assert_eq!(stdtri(100, 0.65), 0.38642898040767143);
        assert_eq!(stdtri(100, 0.75 - 1e-10), 0.6769510426949146);
    }

    #[test]
    fn stdtri_p_large() { // p <= 0.25 or 0.75 <= p

        // Note: These results do not match scipy's stdtrit, but they appear
        // to be more accurate
        assert_eq!(stdtri(1, 1e-10), -3183098861.8379073);
        assert_eq!(stdtri(1, 0.15), -1.9626105055051508);
        assert_eq!(stdtri(1, 0.25), -1.0000000000000007);
        assert_eq!(stdtri(1, 0.75), 1.0000000000000007);
        assert_eq!(stdtri(1, 0.85), 1.9626105055051506);
        assert_eq!(stdtri(1, 1.0 - 1e-10), 3183098598.467149);

        assert_eq!(stdtri(10, 1e-10), -25.46600802169773);
        assert_eq!(stdtri(10, 0.15), -1.0930580735905255);
        assert_eq!(stdtri(10, 0.25), -0.6998120613124323);
        assert_eq!(stdtri(10, 0.75), 0.6998120613124323);
        assert_eq!(stdtri(10, 0.85), 1.0930580735905255);
        assert_eq!(stdtri(10, 1.0 - 1e-10), 25.466007808016013);

        assert_eq!(stdtri(100, 1e-10), -7.083375481400722);
        assert_eq!(stdtri(100, 0.15), -1.0418359009083464);
        assert_eq!(stdtri(100, 0.25), -0.6769510430114689);
        assert_eq!(stdtri(100, 0.75), 0.6769510430114689);
        assert_eq!(stdtri(100, 0.85), 1.0418359009083464);
        assert_eq!(stdtri(100, 1.0 - 1e-10), 7.083375464183688);
    }
}