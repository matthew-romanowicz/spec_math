/*                                                     yn.c
*
*     Bessel function of second kind of integer order
*
*
*
* SYNOPSIS:
*
* double x, y, yn();
* int n;
*
* y = yn( n, x );
*
*
*
* DESCRIPTION:
*
* Returns Bessel function of order n, where n is a
* (possibly negative) integer.
*
* The function is evaluated by forward recurrence on
* n, starting with values computed by the routines
* y0() and y1().
*
* If n = 0 or 1 the routine for y0 or y1 is called
* directly.
*
*
*
* ACCURACY:
*
*
*                      Absolute error, except relative
*                      when y > 1:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       30000       3.4e-15     4.3e-16
*
*
* ERROR MESSAGES:
*
*   message         condition      value returned
* yn singularity   x = 0              INFINITY
* yn overflow                         INFINITY
*
* Spot checked against tables for x, n between 0 and 100.
*
*/

/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

use crate::cephes64::j0::y0;
use crate::cephes64::j1::y1;

pub fn yn(n: isize, x: f64) -> f64 {
    //! Bessel function of second kind of integer order
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns Bessel function of order `n`, where `n` is a
    //! (possibly negative) integer.
    //!
    //! The function is evaluated by forward recurrence on
    //! `n`, starting with values computed by the routines
    //! [`cephes64::y0`](crate::cephes64::y0) and [`cephes64::y1`](crate::cephes64::y01).
    //!
    //! If `n = 0` or `1` the routine for [`cephes64::y0`](crate::cephes64::y0) or [`cephes64::y1`](crate::cephes64::y1) is called
    //! directly.
    //!
    //! ## ACCURACY:
    //!
    //! Absolute error, except relative
    //!                      when y > 1:
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 30</td>
    //!     <td>30000</td>
    //!     <td>3.4e-15</td>
    //!     <td>4.3e-16</td>
    //! </tr>
    //!</table>

    let mut n = n;
    let sign = if n < 0 {
        n = -n;
        if (n & 1) == 0 {	/* -1**n */
            1.0
        } else {
            -1.0
        }
    } else {
        1.0
    };
    


    if n == 0 {
        sign * y0(x)
    } else if n == 1 {
        sign * y1(x)
    } else if x == 0.0 { /* test for overflow */
        //sf_error("yn", SF_ERROR_SINGULAR, NULL);
        -f64::INFINITY * sign
    } else if x < 0.0 {
        //sf_error("yn", SF_ERROR_DOMAIN, NULL);
        f64::NAN
    } else {

        /* forward recurrence on n */

        let mut anm2 = y0(x);
        let mut anm1 = y1(x);
        let mut r = 2.0;
        let mut an = 0.0; // Value doesn't matter, just need to declare it
        for _ in 1..n {
            an = r * anm1 / x - anm2;
            anm2 = anm1;
            anm1 = an;
            r += 2.0;
        }
        sign * an
    }
}

#[cfg(test)]
mod yn_tests {
    use super::*;

    #[test]
    fn yn_pos_tests() {
        assert_eq!(yn(3, 1e-20), -5.092958178940651e+60);
        assert_eq!(yn(3, 1.0), -5.821517605964729);
        assert_eq!(yn(3, 10.0), -0.25136265718383743);
        assert_eq!(yn(3, 1e4), -0.007097801307053426);

        assert_eq!(yn(10, 2e-20), -1.1550829149837402e+205);
        assert_eq!(yn(10, 2.0), -129184.54220803932);
        assert_eq!(yn(10, 20.0), -0.043894653515658466);
        assert_eq!(yn(10, 2e4), 0.0009085311471466713);

        assert_eq!(yn(30, 2e-5), -2.8144202539011476e+180);
        assert_eq!(yn(30, 2.0), -2.9132238482189055e+30);
        assert_eq!(yn(30, 20.0), -114.97814626308332);
        assert_eq!(yn(30, 2e4), 0.000796991966969164);
    }

    #[test]
    fn yn_neg_tests() {
        assert_eq!(yn(-1, 1e-20), 6.3661977236758135e+19);
        assert_eq!(yn(-1, 1.0), 0.7812128213002888);
        assert_eq!(yn(-1, 10.0), -0.24901542420695388);
        assert_eq!(yn(-1, 1e4), -0.00709634275253725);

        assert_eq!(yn(-2, 2e-20), -3.183098861837907e+39);
        assert_eq!(yn(-2, 2.0), -0.6174081041906827);
        assert_eq!(yn(-2, 20.0), -0.07919175824563582);
        assert_eq!(yn(-2, 2e4), 0.0009218922962568055);

        assert_eq!(yn(-30, 2e-5), -2.8144202539011476e+180);
        assert_eq!(yn(-30, 2.0), -2.9132238482189055e+30);
        assert_eq!(yn(-30, 20.0), -114.97814626308332);
        assert_eq!(yn(-30, 2e4), 0.000796991966969164);
    }

    #[test]
    fn yn_2_tests() {
        assert_eq!(yn(2, 1e-20), -1.2732395447351627e+40);
        assert_eq!(yn(2, 1e-10), -1.2732395447351627e+20);
        assert_eq!(yn(2, 1e-5), -12732395447.669935);
        assert_eq!(yn(2, 0.1), -127.64478324269018);
        assert_eq!(yn(2, 1.0), -1.6506826068162546);
        assert_eq!(yn(2, 10.0), -0.005868082442208836);
        assert_eq!(yn(2, 100.0), 0.07683686712502791);
        assert_eq!(yn(2, 1000.0), -0.00476548664020829);
        assert_eq!(yn(2, 1e5), -0.0018467317746558585); // slightly off from scipy
        // TODO: Accuracy degrades
        // assert_eq!(yn(2, 1e10), 7.676508871255354e-06);
        // assert_eq!(yn(2, 1e20), 5.828155155322492e-11);
    }

    #[test]
    fn yn_0_trivials() {
        assert_eq!(yn(0, 0.0), -f64::INFINITY);
        assert_eq!(yn(0, -1e-20).is_nan(), true);
        assert_eq!(yn(0, -10.0).is_nan(), true);
        assert_eq!(yn(0, -f64::INFINITY).is_nan(), true);

        assert_eq!(yn(10, 0.0), -f64::INFINITY);
        assert_eq!(yn(10, -1e-20).is_nan(), true);
        assert_eq!(yn(10, -10.0).is_nan(), true);
        assert_eq!(yn(10, -f64::INFINITY).is_nan(), true);

        assert_eq!(yn(-10, 0.0), -f64::INFINITY);
        assert_eq!(yn(-10, -1e-20).is_nan(), true);
        assert_eq!(yn(-10, -10.0).is_nan(), true);
        assert_eq!(yn(-10, -f64::INFINITY).is_nan(), true);
    }

    #[test]
    fn yn_0_small_x() {
        assert_eq!(yn(0, 1e-20), -29.3912282502858);
        assert_eq!(yn(0, 1e-10), -14.732516272697243);
        assert_eq!(yn(0, 0.1), -1.5342386513503667);
        assert_eq!(yn(0, 1.0), 0.08825696421567697);
        assert_eq!(yn(0, 2.0), 0.5103756726497451);
        assert_eq!(yn(0, 3.0), 0.3768500100127906);
        assert_eq!(yn(0, 4.0), -0.016940739325064846);
        assert_eq!(yn(0, 5.0), -0.30851762524903303);
    }

    #[test]
    fn yn_0_large_x() {
        assert_eq!(yn(0, 5.0 + 1e-16), -0.30851762524903303);
        assert_eq!(yn(0, 6.0), -0.28819468398157916);
        assert_eq!(yn(0, 10.0), 0.05567116728359961);
        assert_eq!(yn(0, 100.0), -0.0772443133650831);
        assert_eq!(yn(0, 1000.0), 0.004715917977623586);
        // TODO: Accuracy reduces for very large x
        //assert_eq!(yn(0, 1e10), -7.676508871690473e-06);
        //assert_eq!(yn(0, 1e20), -5.828155155322492e-11);
    }

    #[test]
    fn yn_1_trivials() {
        assert_eq!(yn(1, 0.0), -f64::INFINITY);
        assert_eq!(yn(1, -1e-20).is_nan(), true);
        assert_eq!(yn(1, -10.0).is_nan(), true);
        assert_eq!(yn(1, -f64::INFINITY).is_nan(), true);
    }

    #[test]
    fn yn_1_small_x() {
        assert_eq!(yn(1, 1e-20), -6.3661977236758135e+19);
        assert_eq!(yn(1, 1e-10), -6366197723.675814);
        assert_eq!(yn(1, 0.1), -6.458951094702027);
        assert_eq!(yn(1, 1.0), -0.7812128213002888);
        assert_eq!(yn(1, 2.0), -0.10703243154093756);
        assert_eq!(yn(1, 3.0), 0.3246744247918001);
        assert_eq!(yn(1, 4.0), 0.3979257105571003);
        assert_eq!(yn(1, 5.0), 0.14786314339122691);
    }

    #[test]
    fn yn_1_large_x() {
        assert_eq!(yn(1, 5.0 + 1e-16), 0.14786314339122691);
        assert_eq!(yn(1, 6.0), -0.17501034430039827);
        assert_eq!(yn(1, 10.0), 0.24901542420695388);
        assert_eq!(yn(1, 100.0), -0.02037231200275932);
        assert_eq!(yn(1, 1000.0), -0.024784331292351868);
        // TODO: Accuracy reduces for very large x
        //assert_eq!(yn(1, 1e10), -2.1755990258465346e-06);
        //assert_eq!(yn(1, 1e20), -5.828155155322492e-11);
    }
}