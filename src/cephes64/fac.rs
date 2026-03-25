/*							fac.c
 *
 *	Factorial function
 *
 *
 *
 * SYNOPSIS:
 *
 * double y, fac();
 * int i;
 *
 * y = fac( i );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns factorial of i  =  1 * 2 * 3 * ... * i.
 * fac(0) = 1.0.
 *
 * Due to machine arithmetic bounds the largest value of
 * i accepted is 33 in DEC arithmetic or 170 in IEEE
 * arithmetic.  Greater values, or negative ones,
 * produce an error message and return MAXNUM.
 *
 *
 *
 * ACCURACY:
 *
 * For i < 34 the values are simply tabulated, and have
 * full machine accuracy.  If i > 55, fac(i) = gamma(i+1);
 * see gamma.c.
 *
 *                      Relative error:
 * arithmetic   domain      peak
 *    IEEE      0, 170    1.4e-15
 *    DEC       0, 33      1.4e-17
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Factorials of integers from 0 through 33 */
const FACTBL: [f64; 34] = [
  1.00000000000000000000E0,
  1.00000000000000000000E0,
  2.00000000000000000000E0,
  6.00000000000000000000E0,
  2.40000000000000000000E1,
  1.20000000000000000000E2,
  7.20000000000000000000E2,
  5.04000000000000000000E3,
  4.03200000000000000000E4,
  3.62880000000000000000E5,
  3.62880000000000000000E6,
  3.99168000000000000000E7,
  4.79001600000000000000E8,
  6.22702080000000000000E9,
  8.71782912000000000000E10,
  1.30767436800000000000E12,
  2.09227898880000000000E13,
  3.55687428096000000000E14,
  6.40237370572800000000E15,
  1.21645100408832000000E17,
  2.43290200817664000000E18,
  5.10909421717094400000E19,
  1.12400072777760768000E21,
  2.58520167388849766400E22,
  6.20448401733239439360E23,
  1.55112100433309859840E25,
  4.03291461126605635584E26,
  1.0888869450418352160768E28,
  3.04888344611713860501504E29,
  8.841761993739701954543616E30,
  2.6525285981219105863630848E32,
  8.22283865417792281772556288E33,
  2.6313083693369353016721801216E35,
  8.68331761881188649551819440128E36
];
const MAXFAC: usize = 170;


use crate::cephes64::gamma::gamma;

pub fn fac(i: usize) -> f64 {
    //! Factorial function
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns factorial of 'i  =  1 * 2 * 3 * ... * i`.
    //!
    //! `fac(0) = 1.0`.
    //!
    //! Due to machine arithmetic bounds the largest value of
    //! `i` accepted is `170`. Greater values  return `f64::INFINITY`.
    //!
    //! ## ACCURACY:
    //!
    //! For `i < 34` the values are simply tabulated, and have
    //! full machine accuracy.  If `i > 55`, `fac(i) = gamma(i+1)`.
    //!
    //! See [`cephes64::gamma`](crate::cephes64::gamma).
    //!
    //! ## ACCURACY:
    //!
    //! Relative error.
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>Domain</th>
    //!     <th>Peak</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 170</td>
    //!     <td>1.4e-15</td>
    //! </tr>
    //!</table>
 

    if ( i > MAXFAC ) {
        // mtherr( "fac", OVERFLOW );
        // return( MAXNUM );
        return f64::INFINITY;
	}

    /* Get answer from table for small i. */
    if ( i < 34 )
	{
	    return FACTBL[i];
	}
    /* Use gamma function for large i. */
    if ( i > 55 )
    {
        let x = i + 1;
        return gamma(x as f64);
    }
    /* Compute directly for intermediate i. */
    let mut n = 34.0;
    let mut f = 34.0;
    for j in 35..=i //( j=35; j<=i; j++ )
    {
        n += 1.0;
        f *= n;
    }
    f *= FACTBL[33];
    return f;
}

#[cfg(test)]
mod fac_tests {
    use super::*;

    #[test]
    fn fac_tests() {
        assert_eq!(fac(0), 1.0);
        assert_eq!(fac(1), 1.0);
        assert_eq!(fac(16), 20922789888000.0);
        assert_eq!(fac(33), 8.683317618811886e+36);
        assert_eq!(fac(34), 2.9523279903960412e+38);
        assert_eq!(fac(35), 1.0333147966386144e+40);
        assert_eq!(fac(55), 1.269640335365828e+73);
        assert_eq!(fac(56), 7.109985878048636e+74);
        assert_eq!(fac(170), 7.257415615308e+306);
        assert_eq!(fac(171), f64::INFINITY);
    }
}