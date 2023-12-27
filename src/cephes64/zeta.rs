/*                                                     zeta.c
*
*     Riemann zeta function of two arguments
*
*
*
* SYNOPSIS:
*
* double x, q, y, zeta();
*
* y = zeta( x, q );
*
*
*
* DESCRIPTION:
*
*
*
*                 inf.
*                  -        -x
*   zeta(x,q)  =   >   (k+q)
*                  -
*                 k=0
*
* where x > 1 and q is not a negative integer or zero.
* The Euler-Maclaurin summation formula is used to obtain
* the expansion
*
*                n
*                -       -x
* zeta(x,q)  =   >  (k+q)
*                -
*               k=1
*
*           1-x                 inf.  B   x(x+1)...(x+2j)
*      (n+q)           1         -     2j
*  +  ---------  -  -------  +   >    --------------------
*        x-1              x      -                   x+2j+1
*                   2(n+q)      j=1       (2j)! (n+q)
*
* where the B2j are Bernoulli numbers.  Note that (see zetac.c)
* zeta(x,1) = zetac(x) + 1.
*
*
*
* ACCURACY:
*
*
*
* REFERENCE:
*
* Gradshteyn, I. S., and I. M. Ryzhik, Tables of Integrals,
* Series, and Products, p. 1073; Academic Press, 1980.
*
*/

/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1984, 1987 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

use crate::cephes64::consts::MACHEP;

/* Expansion coefficients
* for Euler-Maclaurin summation formula
* (2k)! / B2k
* where B2k are Bernoulli numbers
*/
const A: [f64; 12] = [
    12.0,
    -720.0,
    30240.0,
    -1209600.0,
    47900160.0,
    -1.8924375803183791606e9,	/*1.307674368e12/691 */
    7.47242496e10,
    -2.950130727918164224e12,	/*1.067062284288e16/3617 */
    1.1646782814350067249e14,	/*5.109094217170944e18/43867 */
    -4.5979787224074726105e15,	/*8.028576626982912e20/174611 */
    1.8152105401943546773e17,	/*1.5511210043330985984e23/854513 */
    -7.1661652561756670113e18	/*1.6938241367317436694528e27/236364091 */
];

/* 30 Nov 86 -- error in third coefficient fixed */


pub fn zeta(x: f64, q: f64) -> f64 {

    if x == 1.0 {
        return f64::INFINITY;
    } else if (x < 1.0) {
        //sf_error("zeta", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    } else if (q <= 0.0) {
        if q == q.floor() {
            //sf_error("zeta", SF_ERROR_SINGULAR, NULL);
            return f64::INFINITY;
        } else if (x != x.floor()) {
            return f64::NAN; /* because q^-x not defined */
        }
    }

    /* Asymptotic expansion
    * https://dlmf.nist.gov/25.11#E43
    */
    if (q > 1e8) {
        return (1.0 / (x - 1.0) + 1.0 / (2.0 * q)) * q.powf(1.0 - x);
    }

    /* Euler-Maclaurin summation formula */

    /* Permit negative q but continue sum until n+q > +9 .
    * This case should be handled by a reflection formula.
    * If q<0 and x is an integer, there is a relation to
    * the polyGamma function.
    */
    let mut s = q.powf(-x);
    let mut a = q;
    let mut i = 0;
    let mut b = 0.0;
    while (i < 9) || (a <= 9.0) {
        i += 1;
        a += 1.0;
        b = a.powf(-x);
        s += b;
        if (b / s).abs() < MACHEP {
            return s;
        }
    }

    let w = a;
    s += b * w / (x - 1.0);
    s -= 0.5 * b;
    a = 1.0;
    let mut k = 0.0;
    for i in 0..12 {
        a *= x + k;
        b /= w;
        let t = a * b / A[i];
        s = s + t;
        if (t / s).abs() < MACHEP {
            return s;
        }
        k += 1.0;
        a *= x + k;
        b /= w;
        k += 1.0;
    }

    // TODO: Find a way to test this case
    return s;
}

#[cfg(test)]
mod zeta_tests {
    use super::*;

    #[test]
    fn zeta_trivials() {
        assert_eq!(zeta(1.0, 1.0), f64::INFINITY);
        assert_eq!(zeta(1.0 - 1e-10, 1.0).is_nan(), true);
        assert_eq!(zeta(3.0, 0.0), f64::INFINITY);
        assert_eq!(zeta(3.0, -10.0), f64::INFINITY);
        assert_eq!(zeta(3.1, -10.5).is_nan(), true);
    }

    #[test]
    fn zeta_large_q() {
        assert_eq!(zeta(1.5, 1.1e8), 0.00019069251828251056);
        assert_eq!(zeta(10.0, 1.1e8), 4.712195952465921e-74);
        assert_eq!(zeta(1e10, 1.1e10), 0.0);
        assert_eq!(zeta(10.0, 1.1e15), 4.712195759694296e-137);
        assert_eq!(zeta(10.0, 1.1e25), 4.712195759694275e-227);
        assert_eq!(zeta(10.0, 1.1e35), 4.7121956e-317);
    }

    #[test]
    fn zeta_neg_q() {
        assert_eq!(zeta(25.0, -0.9), 1.0000000000000055e+25);
        assert_eq!(zeta(20.0, -1.5), 2097152.0006014686);
        assert_eq!(zeta(100.0, -20.1), 9.99999999998579e+99);
        assert_eq!(zeta(10.0, -1e5 + 0.1), 9999999997.435017);
        assert_eq!(zeta(100.0, -1e5 + 0.1), 9.999999941792339e+99);
    }

    #[test]
    fn zeta_nominal() {
        assert_eq!(zeta(1.0 + 1e-15, 1e-20), 1.0000000000000511e+20);
        assert_eq!(zeta(1.001, 0.002), 1503.6909814191445);
        assert_eq!(zeta(1.0 + 1e-10, 1e-15), 1000000003453880.9);
        assert_eq!(zeta(1.01, 0.01), 205.27439873109253);
        assert_eq!(zeta(1.01, 1e5), 89.12509826963665);
        assert_eq!(zeta(1.5, 0.9e7), 0.0006666666851851857);
        assert_eq!(zeta(25.0, 0.9e7), 5.2235903489005884e-169);
        assert_eq!(zeta(25.0, 0.01), 9.999999999999995e+49);
        assert_eq!(zeta(100.0, 0.01), 9.99999999999998e+199);
    }
}