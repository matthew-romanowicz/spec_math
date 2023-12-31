/*                                                     bdtr.c
*
*     Binomial distribution
*
*
*
* SYNOPSIS:
*
* int k, n;
* double p, y, bdtr();
*
* y = bdtr( k, n, p );
*
* DESCRIPTION:
*
* Returns the sum of the terms 0 through k of the Binomial
* probability density:
*
*   k
*   --  ( n )   j      n-j
*   >   (   )  p  (1-p)
*   --  ( j )
*  j=0
*
* The terms are not summed directly; instead the incomplete
* beta integral is employed, according to the formula
*
* y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
*
* The arguments must be positive, with p ranging from 0 to 1.
*
* ACCURACY:
*
* Tested at random points (a,b,p), with p between 0 and 1.
*
*               a,b                     Relative error:
* arithmetic  domain     # trials      peak         rms
*  For p between 0.001 and 1:
*    IEEE     0,100       100000      4.3e-15     2.6e-16
* See also incbet.c.
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* bdtr domain         k < 0            0.0
*                     n < k
*                     x < 0, x > 1
*/
/*							bdtrc()
*
*	Complemented binomial distribution
*
*
*
* SYNOPSIS:
*
* int k, n;
* double p, y, bdtrc();
*
* y = bdtrc( k, n, p );
*
* DESCRIPTION:
*
* Returns the sum of the terms k+1 through n of the Binomial
* probability density:
*
*   n
*   --  ( n )   j      n-j
*   >   (   )  p  (1-p)
*   --  ( j )
*  j=k+1
*
* The terms are not summed directly; instead the incomplete
* beta integral is employed, according to the formula
*
* y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
*
* The arguments must be positive, with p ranging from 0 to 1.
*
* ACCURACY:
*
* Tested at random points (a,b,p).
*
*               a,b                     Relative error:
* arithmetic  domain     # trials      peak         rms
*  For p between 0.001 and 1:
*    IEEE     0,100       100000      6.7e-15     8.2e-16
*  For p between 0 and .001:
*    IEEE     0,100       100000      1.5e-13     2.7e-15
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* bdtrc domain      x<0, x>1, n<k       0.0
*/
/*							bdtri()
*
*	Inverse binomial distribution
*
*
*
* SYNOPSIS:
*
* int k, n;
* double p, y, bdtri();
*
* p = bdtri( k, n, y );
*
* DESCRIPTION:
*
* Finds the event probability p such that the sum of the
* terms 0 through k of the Binomial probability density
* is equal to the given cumulative probability y.
*
* This is accomplished using the inverse beta integral
* function and the relation
*
* 1 - p = incbi( n-k, k+1, y ).
*
* ACCURACY:
*
* Tested at random points (a,b,p).
*
*               a,b                     Relative error:
* arithmetic  domain     # trials      peak         rms
*  For p between 0.001 and 1:
*    IEEE     0,100       100000      2.3e-14     6.4e-16
*    IEEE     0,10000     100000      6.6e-12     1.2e-13
*  For p between 10^-6 and 0.001:
*    IEEE     0,100       100000      2.0e-12     1.3e-14
*    IEEE     0,10000     100000      1.5e-12     3.2e-14
* See also incbi.c.
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* bdtri domain     k < 0, n <= k         0.0
*                  x < 0, x > 1
*/

/*                                                             bdtr() */


/*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

use crate::cephes64::unity::{expm1, log1p};
use crate::cephes64::incbet::incbet;
use crate::cephes64::incbi::incbi;

pub fn bdtrc(k: f64, n: isize, p: f64) -> f64
{
    //double dk, dn;
    let fk = k.floor();

    if p.is_nan() || k.is_nan() {
        return f64::NAN;
    }

    if !(0.0..=1.0).contains(&p) || (n as f64) < fk {
        //sf_error("bdtrc", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    }

    if fk < 0.0 {
        return 1.0;
    }

    if fk == n as f64 {
        return 0.0;
    }

    let dn = n as f64 - fk;
    if k == 0.0 { // TODO: Should this be fk like in bdtr?
        if p < 0.01 {
            -expm1(dn * log1p(-p))
        } else {
            1.0 - (1.0 - p).powf(dn)
        }
    } else {
        incbet(fk + 1.0, dn, p)
    }
}



pub fn bdtr(k: f64, n: isize, p: f64) -> f64 {
    //double dk, dn;
    let fk = k.floor();

    if p.is_nan() || k.is_nan() {
        return f64::NAN;
    }

    if !(0.0..=1.0).contains(&p) || fk < 0.0 || (n as f64) < fk {
        //sf_error("bdtr", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    }

    if fk == n as f64 {
        return 1.0;
    }

    let dn = n as f64 - fk;
    if fk == 0.0 {
        (1.0 - p).powf(dn)
    }
    else {
        incbet(dn, fk + 1.0, 1.0 - p)
    }
}


pub fn bdtri(k: f64, n: isize, y: f64) -> f64 {
    //double p, dn, dk;
    let fk = k.floor();

    if k.is_nan() {
        return f64::NAN;
    }

    if !(0.0..=1.0).contains(&y) || fk < 0.0 || (n as f64) <= fk {
        //sf_error("bdtri", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    }

    let dn = n as f64 - fk;

    if fk == 0.0 {
        if y > 0.8 {
            -expm1(log1p(y - 1.0) / dn)
        } else {
            1.0 - y.powf(1.0 / dn)
        }
    } else {
        let dk = fk + 1.0;
        let p = incbet(dn, dk, 0.5);
        if p > 0.5 {
            incbi(dk, dn, 1.0 - y)
        } else {
            1.0 - incbi(dn, dk, y)
        }
    }
}

#[cfg(test)]
mod bdtrc_tests {
    use super::*;

    #[test]
    fn bdtrc_trivials() {
        assert_eq!(bdtrc(f64::NAN, 10, 0.5).is_nan(), true);
        assert_eq!(bdtrc(5.0, 10, f64::NAN).is_nan(), true);
        assert_eq!(bdtrc(5.0, 10, -1e-10).is_nan(), true);
        assert_eq!(bdtrc(5.0, 10, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(bdtrc(11.0, 10, 0.5).is_nan(), true);
        assert_eq!(bdtrc(10.0, 10, 0.5), 0.0);
        assert_eq!(bdtrc(10.9, 10, 0.5), 0.0);
        assert_eq!(bdtrc(-1.0, 10, 0.5), 1.0); 
    }

    #[test]
    fn bdtrc_k_0() {
        assert_eq!(bdtrc(0.0, 0, 1e-15), 0.0);
        assert_eq!(bdtrc(0.0, 1, 1e-15), 1.0000000000000003e-15);
        assert_eq!(bdtrc(0.0, 10, 1e-15), 9.999999999999956e-15);
        assert_eq!(bdtrc(0.0, 2147483647, 1e-15), 2.1474813411586448e-06);

        assert_eq!(bdtrc(0.0, 0, 0.5), 0.0);
        assert_eq!(bdtrc(0.0, 1, 0.4), 0.4);
        assert_eq!(bdtrc(0.0, 10, 0.6), 0.9998951424);
        assert_eq!(bdtrc(0.0, 2147483647, 0.8), 1.0);

        assert_eq!(bdtrc(0.0, 0, 1.0 - 1e-10), 0.0);
        assert_eq!(bdtrc(0.0, 1, 1.0 - 1e-10), 0.9999999999);
        assert_eq!(bdtrc(0.0, 10, 1.0 - 1e-10), 1.0);
        assert_eq!(bdtrc(0.0, 2147483647, 1.0 - 1e-10), 1.0);
    }

    #[test]
    fn bdtrc_nominal() {
        assert_eq!(bdtrc(0.9, 1, 1e-15), 1.0e-15);
        assert_eq!(bdtrc(0.9, 2147483647, 1e-15), 2.1474813411586494e-6);
        assert_eq!(bdtrc(0.9, 0, 0.5), 0.0);
        assert_eq!(bdtrc(0.9, 10, 0.6), 0.9998951424);
        assert_eq!(bdtrc(0.9, 1, 1.0 - 1e-10), 0.9999999999);
        assert_eq!(bdtrc(0.9, 10, 1.0 - 1e-10), 0.9999999999999999);

        assert_eq!(bdtrc(1.9, 2, 1e-15), 1e-30);
        assert_eq!(bdtrc(5.6, 10, 1e-15), 2.0999999999999935e-88);
        assert_eq!(bdtrc(2147483646.0, 2147483647, 1e-15), 0.0);

        assert_eq!(bdtrc(1.9, 2, 0.5), 0.25);
        assert_eq!(bdtrc(5.6, 10, 0.4), 0.1662386176);
        assert_eq!(bdtrc(1718000000.0, 2147483647, 0.8), 0.24015971937437247);

        assert_eq!(bdtrc(1.9, 2, 0.9), 0.81);
        assert_eq!(bdtrc(8.6, 10, 0.2), 4.198400000000002e-06);
        assert_eq!(bdtrc(644000000.0, 2147483647, 0.3), 0.9999999999999999);
        assert_eq!(bdtrc(645000000.0, 2147483647, 0.3), 5.085779692837519e-277);

        assert_eq!(bdtrc(1.9, 2, 1.0 - 1e-10), 0.9999999998);
        assert_eq!(bdtrc(5.6, 10, 1.0 - 1e-10), 0.9999999999999999);
        assert_eq!(bdtrc(1718000000.0, 2147483647, 1.0 - 1e-10), 0.9999999999999999);
    }
}

#[cfg(test)]
mod bdtr_tests {
    use super::*;

    #[test]
    fn bdtr_trivials() {
        assert_eq!(bdtr(f64::NAN, 10, 0.5).is_nan(), true);
        assert_eq!(bdtr(5.0, 10, f64::NAN).is_nan(), true);
        assert_eq!(bdtr(5.0, 10, -1e-10).is_nan(), true);
        assert_eq!(bdtr(5.0, 10, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(bdtr(11.0, 10, 0.5).is_nan(), true);
        assert_eq!(bdtr(10.0, 10, 0.5), 1.0);
        assert_eq!(bdtr(10.9, 10, 0.5), 1.0);

        // TODO: Should this be NAN? k rounds to 0
        assert_eq!(bdtr(-1e-10, 10, 0.5).is_nan(), true); 
    }

    #[test]
    fn bdtr_fk_0() {
        assert_eq!(bdtr(0.0, 0, 1e-15), 1.0);
        assert_eq!(bdtr(0.9, 1, 1e-15), 0.999999999999999);
        assert_eq!(bdtr(0.0, 10, 1e-15), 0.99999999999999);
        assert_eq!(bdtr(0.9, 2147483647, 1e-15), 0.9999978542350912);

        assert_eq!(bdtr(0.9, 0, 0.5), 1.0);
        assert_eq!(bdtr(0.0, 1, 0.4), 0.6);
        assert_eq!(bdtr(0.9, 10, 0.6), 0.00010485760000000006);
        assert_eq!(bdtr(0.0, 2147483647, 0.8), 0.0);

        assert_eq!(bdtr(0.0, 0, 1.0 - 1e-10), 1.0);
        assert_eq!(bdtr(0.9, 1, 1.0 - 1e-10), 1.000000082740371e-10);
        assert_eq!(bdtr(0.9, 10, 1.0 - 1e-10), 1.000000827404018e-100);
        assert_eq!(bdtr(0.0, 2147483647, 1.0 - 1e-10), 0.0);
    }

    #[test]
    fn bdtr_nominal() {
        assert_eq!(bdtr(1.9, 2, 1e-15), 0.9999999999999999);
        assert_eq!(bdtr(5.6, 10, 1e-15), 0.9999999999999999);
        assert_eq!(bdtr(2147483646.0, 2147483647, 1e-15), 0.9999999999999999);

        assert_eq!(bdtr(1.9, 2, 0.5), 0.75);
        assert_eq!(bdtr(5.6, 10, 0.4), 0.8337613824);
        assert_eq!(bdtr(1718000000.0, 2147483647, 0.8), 0.7598402806256275);

        assert_eq!(bdtr(1.9, 2, 0.9), 0.18999999999999995);
        assert_eq!(bdtr(8.6, 10, 0.2), 0.9999958016);
        assert_eq!(bdtr(644000000.0, 2147483647, 0.3), 4.06105281232279e-31);

        assert_eq!(bdtr(1.9, 2, 1.0 - 1e-10), 2.000000165380742e-10);
        assert_eq!(bdtr(5.6, 10, 1.0 - 1e-10), 2.5200010414788468e-48);
        assert_eq!(bdtr(1718000000.0, 2147483647, 1.0 - 1e-10), 0.0);
    }
}

#[cfg(test)]
mod bdtri_tests {
    use super::*;

    #[test]
    fn bdtri_trivials() {
        assert_eq!(bdtri(f64::NAN, 10, 0.5).is_nan(), true);
        assert_eq!(bdtri(5.0, 10, f64::NAN).is_nan(), true);
        assert_eq!(bdtri(5.0, 10, -1e-10).is_nan(), true);
        assert_eq!(bdtri(5.0, 10, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(bdtri(11.0, 10, 0.5).is_nan(), true);
        assert_eq!(bdtri(10.0, 10, 0.5).is_nan(), true);
        assert_eq!(bdtri(10.9, 10, 0.5).is_nan(), true);

        // TODO: Should this be NAN? k rounds to 0
        assert_eq!(bdtri(-1e-10, 10, 0.5).is_nan(), true); 
    }

    #[test]
    fn bdtri_fk_0() {
        assert_eq!(bdtri(0.9, 1, 1e-15), 0.999999999999999);
        assert_eq!(bdtri(0.0, 10, 1e-15), 0.9683772233983162);
        assert_eq!(bdtri(0.9, 2147483647, 1e-15), 1.6083370968367205e-08);

        assert_eq!(bdtri(0.0, 1, 0.4), 0.6);
        assert_eq!(bdtri(0.9, 10, 0.6), 0.04979978349432357);
        assert_eq!(bdtri(0.0, 2147483647, 0.8), 1.039093255883472e-10);

        assert_eq!(bdtri(0.9, 1, 1.0 - 1e-10), 1.000000082740371e-10);
        assert_eq!(bdtri(0.9, 10, 1.0 - 1e-10), 1.0000000827853709e-11);
        assert_eq!(bdtri(0.0, 2147483647, 1.0 - 1e-10), 4.6566132607685043e-20);
    }

    #[test]
    fn bdtri_nominal() {
        assert_eq!(bdtri(1.9, 2, 1e-15), 0.9999999999999996);
        assert_eq!(bdtri(5.6, 10, 1e-15), 0.9996690440210801);
        assert_eq!(bdtri(2147483646.0, 2147483647, 1e-15), 0.9999999999999999);

        assert_eq!(bdtri(1.9, 2, 0.5), 0.7071067811865475);
        assert_eq!(bdtri(5.6, 10, 0.4), 0.5869309013035198);
        assert_eq!(bdtri(1718000000.0, 2147483647, 0.8), 0.7999988275359662);

        assert_eq!(bdtri(1.9, 2, 0.9), 0.3162277660168379);
        assert_eq!(bdtri(8.6, 10, 0.2), 0.916740063444129);
        assert_eq!(bdtri(644000000.0, 2147483647, 0.3), 0.2998910545694471);

        assert_eq!(bdtri(1.9, 2, 1.0 - 1e-10), 1.0000000413701846e-05);
        assert_eq!(bdtri(5.6, 10, 1.0 - 1e-10), 0.008881962019336944);
        assert_eq!(bdtri(1718000000.0, 2147483647, 1.0 - 1e-10), 0.7999511799804422);
    }
}