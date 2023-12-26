/*                                                     kn.c
*
*     Modified Bessel function, third kind, integer order
*
*
*
* SYNOPSIS:
*
* double x, y, kn();
* int n;
*
* y = kn( n, x );
*
*
*
* DESCRIPTION:
*
* Returns modified Bessel function of the third kind
* of order n of the argument.
*
* The range is partitioned into the two intervals [0,9.55] and
* (9.55, infinity).  An ascending power series is used in the
* low range, and an asymptotic expansion in the high range.
*
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0,30        90000       1.8e-8      3.0e-10
*
*  Error is high only near the crossover point x = 9.55
* between the two expansions used.
*/


/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
*/


/*
* Algorithm for Kn.
*                        n-1
*                    -n   -  (n-k-1)!    2   k
* K (x)  =  0.5 (x/2)     >  -------- (-x /4)
*  n                      -     k!
*                        k=0
*
*                     inf.                                   2   k
*        n         n   -                                   (x /4)
*  + (-1)  0.5(x/2)    >  {p(k+1) + p(n+k+1) - 2log(x/2)} ---------
*                      -                                  k! (n+k)!
*                     k=0
*
* where  p(m) is the psi function: p(1) = -EUL and
*
*                       m-1
*                        -
*       p(m)  =  -EUL +  >  1/k
*                        -
*                       k=1
*
* For large x,
*                                          2        2     2
*                                       u-1     (u-1 )(u-3 )
* K (z)  =  sqrt(pi/2z) exp(-z) { 1 + ------- + ------------ + ...}
*  v                                        1            2
*                                     1! (8z)     2! (8z)
* asymptotically, where
*
*            2
*     u = 4 v .
*
*/

const EUL: f64 =  5.772156649015328606065e-1;
const MAXFAC: isize = 31;

use crate::cephes64::consts::{M_PI, MACHEP, MAXLOG};

pub fn kn(nn: isize, x: f64) -> f64 {
    //! Modified Bessel function, third kind, integer order
    //! 
    //! ## DESCRIPTION:
    //! 
    //! Returns modified Bessel function of the third kind
    //! of order n of the argument.
    //!
    //! The range is partitioned into the two intervals [0,9.55] and
    //! (9.55, infinity).  An ascending power series is used in the
    //! low range, and an asymptotic expansion in the high range.
    //!
    //! ## ACCURACY:
    //!
    //! Relative error:
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
    //!     <td>90000</td>
    //!     <td>1.8e-8</td>
    //!     <td>3.0e-10</td>
    //! </tr>
    //!</table>
    //!
    //! Error is high only near the crossover point x = 9.55
    //! between the two expansions used.

    let n = nn.abs();

    if n > MAXFAC {
        //overf:
        //sf_error("kn", SF_ERROR_OVERFLOW, NULL);
        return f64::INFINITY;
    } else if x <= 0.0 {
        if x < 0.0 {
            //sf_error("kn", SF_ERROR_DOMAIN, NULL);
            return f64::NAN;
        }
        else {
            //sf_error("kn", SF_ERROR_SINGULAR, NULL);
            return f64::INFINITY;
        }
    }


    if x > 9.55 {
        /* Asymptotic expansion for Kn(x) */
        /* Converges to 1.4e-17 for x > 18.4 */
        if x > MAXLOG {
            //sf_error("kn", SF_ERROR_UNDERFLOW, NULL);
            return 0.0;
        }
        let k = n;
        let pn = 4.0 * (k as f64 * k as f64);
        let mut pk = 1.0;
        let z0 = 8.0 * x;
        let mut fn_ = 1.0;
        let mut t: f64 = 1.0;
        let mut s = t;
        let mut nkf = f64::INFINITY;
        let mut i = 0;
        loop {
            let z = pn - pk * pk;
            t = t * z / (fn_ * z0);
            let nk1f = (t).abs();
            if (i >= n) && (nk1f > nkf) {
                break;
            }
            nkf = nk1f;
            s += t;
            fn_ += 1.0;
            pk += 2.0;
            i += 1;
            if !((t / s).abs() > MACHEP) {
                break;
            }
        }
        let ans = (-x).exp() * (M_PI / (2.0 * x)).sqrt() * s;
        return ans;
    }
    

    let mut ans = 0.0;
    let z0 = 0.25 * x * x;
    let mut fn_ = 1.0;
    let mut pn = 0.0;
    let mut zmn = 1.0;
    let tox = 2.0 / x;

    if n > 0 {
        /* compute factorial of n and psi(n) */
        pn = -EUL;
        let mut k = 1.0;
        for _ in 1..n {
            pn += 1.0 / k;
            k += 1.0;
            fn_ *= k;
        }

        zmn = tox;

        if n == 1 {
            ans = 1.0 / x;
        } else {
            let mut nk1f = fn_ / (n as f64);
            let mut kf = 1.0;
            let mut s = nk1f;
            let z = -z0;
            let mut zn = 1.0;
            let mut t;
            for i in 1..n {
                nk1f = nk1f / (n - i) as f64;
                kf = kf * (i as f64);
                zn *= z;
                t = nk1f * zn / kf;
                s += t;
                if (f64::MAX - t.abs()) < s.abs() {
                    return f64::INFINITY;
                }
                if (tox > 1.0) && ((f64::MAX / tox) < zmn) {
                    return f64::INFINITY;
                }
                zmn *= tox;
            }
            s *= 0.5;
            t = s.abs();
            if (zmn > 1.0) && ((f64::MAX / zmn) < t){
                return f64::INFINITY;
            }
            if (t > 1.0) && ((f64::MAX / t) < zmn) {
                return f64::INFINITY;
            }
            ans = s * zmn;
        }
    }


    let tlg = 2.0 * (0.5 * x).ln();
    let mut pk = -EUL;
    let mut t = if n == 0 {
        pn = pk;
        1.0
    } else {
        pn = pn + 1.0 / (n as f64);
        1.0 / fn_
    };
    let mut s = (pk + pn - tlg) * t;
    let mut k = 1.0;
    loop {
        t *= z0 / (k * (k + n as f64));
        pk += 1.0 / k;
        pn += 1.0 / (k + n as f64);
        s += (pk + pn - tlg) * t;
        k += 1.0;
        if !((t / s).abs() > MACHEP) {
            break;
        }
    }

    s = 0.5 * s / zmn;
    if n & 1 != 0 {
        s = -s;
    }
    
    ans += s;

    return ans;
}

// TODO: This doesn't line up exactly with output from the C CEHPES library
#[cfg(test)]
mod kn_tests {
    use super::*;

    #[test]
    fn yn_pos_tests() {
        assert_eq!(kn(3, 1e-20), 8e+60);
        assert_eq!(kn(3, 1.0), 7.1012628247379448);
        assert_eq!(kn(3, 10.0), 0.000027252700261657538);
        assert_eq!(kn(3, 1e4), 0.0);

        assert_eq!(kn(10, 2e-20), 1.8143999999999999e+205);
        assert_eq!(kn(10, 2.0), 162482.40397955917);
        assert_eq!(kn(10, 20.0), 6.316214528321577e-09);
        assert_eq!(kn(10, 2e4), 0.0);

        assert_eq!(kn(30, 2e-5), 4.420880996854586e+180);
        assert_eq!(kn(30, 2.0), 4.271125754887687e+30);
        assert_eq!(kn(30, 20.0), 0.16883087719470813);
        assert_eq!(kn(30, 2e4), 0.0);
    }

    #[test]
    fn yn_neg_tests() {
        assert_eq!(kn(-3, 1e-20), 8e+60);
        assert_eq!(kn(-3, 1.0), 7.1012628247379448);
        assert_eq!(kn(-3, 10.0), 0.000027252700261657538);
        assert_eq!(kn(-3, 1e4), 0.0);

        assert_eq!(kn(-10, 2e-20), 1.8143999999999999e+205);
        assert_eq!(kn(-10, 2.0), 162482.40397955917);
        assert_eq!(kn(-10, 20.0), 6.316214528321577e-09);
        assert_eq!(kn(-10, 2e4), 0.0);

        assert_eq!(kn(-30, 2e-5), 4.420880996854586e+180);
        assert_eq!(kn(-30, 2.0), 4.271125754887687e+30);
        assert_eq!(kn(-30, 20.0), 0.16883087719470813);
        assert_eq!(kn(-30, 2e4), 0.0);
    }

    #[test]
    fn kn_1_trivials() {
        assert_eq!(kn(1, f64::INFINITY), 0.0);
        assert_eq!(kn(1, -f64::INFINITY).is_nan(), true);
        assert_eq!(kn(1, f64::NAN).is_nan(), true);
        assert_eq!(kn(1, -1e-20).is_nan(), true);
        assert_eq!(kn(1, -1.0).is_nan(), true);
        assert_eq!(kn(1, 0.0), f64::INFINITY);
    }

    #[test]
    fn kn_1_small_x() {
        assert_eq!(kn(1, 1e-20), 1e+20);

        assert_eq!(kn(1, 1e-10), 10000000000.0);
        assert_eq!(kn(1, 0.1), 9.853844780870606);
        assert_eq!(kn(1, 0.5), 1.6564411200033011);
        assert_eq!(kn(1, 1.0), 0.6019072301972346);
        assert_eq!(kn(1, 1.5), 0.27738780045684375);
        assert_eq!(kn(1, 2.0), 0.13986588181652243);
    }

    #[test]
    fn kn_1_large_x() {
        assert_eq!(kn(1, 2.0 + 1e-16), 0.13986588181652243);
        assert_eq!(kn(1, 3.0), 0.04015643112819417);
        assert_eq!(kn(1, 10.0), 1.8648773457308253e-5);
        assert_eq!(kn(1, 20.0), 5.883057969557036e-10);
        assert_eq!(kn(1, 100.0), 4.679853735636909e-45);
        assert_eq!(kn(1, 1000.0), 0.0);
    }

    #[test]
    fn kn_0_trivials() {
        assert_eq!(kn(0, f64::INFINITY), 0.0);
        assert_eq!(kn(0, -f64::INFINITY).is_nan(), true);
        assert_eq!(kn(0, f64::NAN).is_nan(), true);
        assert_eq!(kn(0, -1e-20).is_nan(), true);
        assert_eq!(kn(0, -1.0).is_nan(), true);
        assert_eq!(kn(0, 0.0), f64::INFINITY);
    }

    #[test]
    fn kn_0_small_x() {
        assert_eq!(kn(0, 1e-20), 46.16763337553933);

        assert_eq!(kn(0, 1e-10), 23.141782445598867);
        assert_eq!(kn(0, 0.1), 2.427069024702017);
        assert_eq!(kn(0, 0.5), 0.9244190712276656);
        assert_eq!(kn(0, 1.0), 0.42102443824070834);
        assert_eq!(kn(0, 1.5), 0.21380556264752568);
        assert_eq!(kn(0, 2.0), 0.11389387274953341);
    }

    #[test]
    fn kn_0_large_x() {
        assert_eq!(kn(0, 2.0 + 1e-16), 0.11389387274953341);
        assert_eq!(kn(0, 3.0), 0.03473950438627915);
        assert_eq!(kn(0, 10.0), 1.7780062319409137e-5);
        assert_eq!(kn(0, 20.0), 5.7412378153365248e-10);
        assert_eq!(kn(0, 100.0), 4.656628229175901e-45);
        assert_eq!(kn(0, 1000.0), 0.0);
    }
}