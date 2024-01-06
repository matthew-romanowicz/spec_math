/*                                                     nbdtr.c
*
*     Negative binomial distribution
*
*
*
* SYNOPSIS:
*
* int k, n;
* double p, y, nbdtr();
*
* y = nbdtr( k, n, p );
*
* DESCRIPTION:
*
* Returns the sum of the terms 0 through k of the negative
* binomial distribution:
*
*   k
*   --  ( n+j-1 )   n      j
*   >   (       )  p  (1-p)
*   --  (   j   )
*  j=0
*
* In a sequence of Bernoulli trials, this is the probability
* that k or fewer failures precede the nth success.
*
* The terms are not computed individually; instead the incomplete
* beta integral is employed, according to the formula
*
* y = nbdtr( k, n, p ) = incbet( n, k+1, p ).
*
* The arguments must be positive, with p ranging from 0 to 1.
*
* ACCURACY:
*
* Tested at random points (a,b,p), with p between 0 and 1.
*
*               a,b                     Relative error:
* arithmetic  domain     # trials      peak         rms
*    IEEE     0,100       100000      1.7e-13     8.8e-15
* See also incbet.c.
*
*/
/*
* Cephes Math Library Release 2.3:  March, 1995
* Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

use crate::cephes64::incbet::incbet;
use crate::cephes64::incbi::incbi;

pub fn nbdtrc(k: isize, n: isize, p: f64) -> f64 {
    //! Complemented negative binomial distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the sum of the terms k+1 to infinity of the negative
    //! binomial distribution:
    //!
    #![doc=include_str!("nbdtrc.svg")]
    //!
    //! The terms are not computed individually; instead the incomplete
    //! beta integral is employed, according to the formula
    //!
    //! `y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p )`
    //!
    //! The arguments must be positive, with `p` ranging from 0 to 1.
    //!
    //! ## ACCURACY:
    //!
    //! Tested at random points (`a`, `b`, `p`), with `p` between 0 and 1.
    //!
    //! Relative error
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>a, b Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>1.7e-13</td>
    //!     <td>8.8e-15</td>
    //! </tr>
    //!</table>
    //!
    //! See also [`cephes64::incbet`](crate::cephes64::incbet).

    if !(0.0..=1.0).contains(&p) || k < 0 {
        //sf_error("nbdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let dk = k + 1;
        let dn = n;
        incbet(dk as f64, dn as f64, 1.0 - p)
    }
}



pub fn nbdtr(k: isize, n: isize, p: f64) -> f64 {
    //! Negative binomial distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the sum of the terms 0 through k of the negative
    //! binomial distribution:
    //!
    #![doc=include_str!("nbdtr.svg")]
    //!
    //! In a sequence of Bernoulli trials, this is the probability
    //! that k or fewer failures precede the nth success.
    //!
    //! The terms are not computed individually; instead the incomplete
    //! beta integral is employed, according to the formula
    //!
    //! `y = nbdtr( k, n, p ) = incbet( n, k+1, p )`
    //!
    //! The arguments must be positive, with p ranging from 0 to 1.
    //!
    //! ## ACCURACY:
    //!
    //! Tested at random points (`a`, `b`, `p`), with `p` between 0 and 1.
    //!
    //! Relative error
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>a, b Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>1.7e-13</td>
    //!     <td>8.8e-15</td>
    //! </tr>
    //!</table>
    //!
    //! See also [`cephes64::incbet`](crate::cephes64::incbet).

    if !(0.0..=1.0).contains(&p) || k < 0 {
        //sf_error("nbdtr", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let dk = k + 1;
        let dn = n;
        incbet(dn as f64, dk as f64, p)
    }
}



pub fn nbdtri(k: isize, n: isize, p: f64) -> f64 {
    //! Functional inverse of negative binomial distribution
    //!
    //! ## DESCRIPTION:
    //!
    //! Finds the argument `p` such that `nbdtr(k,n,p)` is equal to `y`.
    //!
    //! ## ACCURACY:
    //!
    //! Tested at random points (`a`, `b`, `y`), with `y` between 0 and 1.
    //!
    //! Relative error:
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>a, b Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 100</td>
    //!     <td>100000</td>
    //!     <td>1.5e-14</td>
    //!     <td>8.5e-16</td>
    //! </tr>
    //!</table>
    //!
    //! See also [`cephes64::incbi`](crate::cephes64::incbi).

    if !(0.0..=1.0).contains(&p) || k < 0 {
        //sf_error("nbdtri", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {
        let dk = k + 1;
        let dn = n;
        incbi(dn as f64, dk as f64, p)
    }
}

// Very few tests due to simplicity of implementation
#[cfg(test)]
mod nbdtr_tests {
    use super::*;

    #[test]
    fn nbdtr_trivials() {
        assert_eq!(nbdtrc(-1, 10, 0.5).is_nan(), true);
        assert_eq!(nbdtrc(1, 10, -1e-15).is_nan(), true);
        assert_eq!(nbdtrc(1, 10, 1.0 + 1e-15).is_nan(), true);

        assert_eq!(nbdtr(-1, 10, 0.5).is_nan(), true);
        assert_eq!(nbdtr(1, 10, -1e-15).is_nan(), true);
        assert_eq!(nbdtr(1, 10, 1.0 + 1e-15).is_nan(), true);

        assert_eq!(nbdtri(-1, 10, 0.5).is_nan(), true);
        assert_eq!(nbdtri(1, 10, -1e-15).is_nan(), true);
        assert_eq!(nbdtri(1, 10, 1.0 + 1e-15).is_nan(), true);
    }

    #[test]
    fn nbdtr_tests() {
        assert_eq!(nbdtrc(5, 10, 0.5), 0.8491210937499999);
        assert_eq!(nbdtrc(5, 10, 0.9), 0.0022496700850479973);
        assert_eq!(nbdtrc(15, 10, 0.5), 0.11476147174835194);
        assert_eq!(nbdtrc(15, 15, 0.1), 0.9999999644052073);
        assert_eq!(nbdtrc(15, 10, 0.0), 1.0);
        assert_eq!(nbdtrc(15, 10, 1.0), 0.0);

        assert_eq!(nbdtr(5, 10, 0.5), 0.15087890625000008);
        assert_eq!(nbdtr(5, 10, 0.9), 0.997750329914952);
        assert_eq!(nbdtr(15, 10, 0.5), 0.8852385282516481);
        assert_eq!(nbdtr(15, 15, 0.1), 3.559479269444176e-08);
        assert_eq!(nbdtr(15, 10, 0.0), 0.0);
        assert_eq!(nbdtr(15, 10, 1.0), 1.0);

        assert_eq!(nbdtri(5, 10, 0.5), 0.6303295486548143);
        assert_eq!(nbdtri(5, 10, 0.9), 0.7744087263163602);
        assert_eq!(nbdtri(15, 10, 0.5), 0.3816135642524775);
        assert_eq!(nbdtri(15, 15, 0.1), 0.36969603948477525);
        assert_eq!(nbdtri(15, 10, 0.0), 0.0);
        assert_eq!(nbdtri(15, 10, 1.0), 1.0);
    }
}