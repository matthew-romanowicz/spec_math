/*                                                     ellik.c
 *
 *     Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellik();
 *
 * y = ellik( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi | m) =   |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with m in [0, 1] and phi as indicated.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10       200000      7.4e-16     1.0e-16
 *
 *
 */


/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
/* Copyright 2014, Eric W. Moore */

/*     Incomplete elliptic integral of first kind      */

#![allow(clippy::excessive_precision)]

use super::consts::{M_PI, M_PI_2, MACHEP};
use super::ellpk::ellpk;

//$$F(\varphi, m) = \int_{0}^{\varphi} \frac{d\theta}{\sqrt{1 - m\,\sin^2(\theta)}}$$
// TODO: Expand domain to match wolfram|alpha
pub fn ellik(phi: f64,  m: f64) -> f64
{
//! Incomplete elliptic integral of the first kind
//!
//! ## Description
//! Approximates the integral 
//!
#![doc=include_str!("ellik.svg")]
//!
//! of amplitude `phi` and modulus `m`, using the arithmetic-geometric mean algorithm.
//!
//! ## Accuracy
//!
//! Tested at random points with `m` in `[0, 1]` and `phi` as indicated.
//!
//! Relative Error:
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
//!     <td>-10, 10</td>
//!     <td>200000</td>
//!     <td>7.4e-16</td>
//!     <td>1.0e-16</td>
//! </tr>
//!</table>

    if phi.is_nan() || m.is_nan() {
        return f64::NAN;
    }
    if m > 1.0 {
        return f64::NAN;
    }
    if phi.is_infinite() || m.is_infinite()
    {
        if m.is_infinite() && !phi.is_infinite() {
            return 0.0;
        }
        else if phi.is_infinite() && !m.is_infinite() {
            return phi;
        }
        else {
            return f64::NAN;
        }
    }
    if m == 0.0 {
        return phi;
    }
    let a = 1.0 - m;
    if a == 0.0 {
        if phi.abs() >= M_PI_2 {
            //sf_error("ellik", SF_ERROR_SINGULAR, NULL);
            return f64::INFINITY;
        }
        /* DLMF 19.6.8, and 4.23.42 */
       return phi.tan().asinh();
    }
    let mut phi = phi;
    let mut npio2 = (phi / M_PI_2).floor();
    if npio2.abs() % 2.0 == 1.0 {
        npio2 += 1.0;
    }
    let mut k = if npio2 != 0.0 {
	    phi -= npio2 * M_PI_2;
        ellpk(a)
    }
    else {
        0.0
    };
    let sign = if phi < 0.0 {
        phi = -phi;
        -1
    }
    else {
        0
    };
    let mut temp = if a > 1.0 {
        ellik_neg_m(phi, m)
    } else {
        let mut b = a.sqrt();
        let mut t = phi.tan();
        if t.abs() > 10.0 {
            /* Transform the amplitude */
            let e = 1.0 / (b * t);
            /* ... but avoid multiple recursions.  */
            if e.abs() < 10.0 {
                let e = e.atan();
                if npio2 == 0.0 {
                    k = ellpk(a);
                }
                let mut temp2 = k - ellik(e, m);
                if sign < 0 {
                    temp2 = -temp2;
                }
                temp2 += npio2 * k;
                return temp2;
            }
        }
        let mut a = 1.0;
        let mut c = m.sqrt();
        let mut d = 1;
        let mut mod_ = 0;

        while (c / a).abs() > MACHEP {
            let temp2 = b / a;
            phi = phi + (t * temp2).atan() + (mod_ as f64) * M_PI;
            let denom = 1.0 - temp2 * t * t;
            if denom.abs() > 10.0 * MACHEP {
                t = t * (1.0 + temp2) / denom;
                mod_ = ((phi + M_PI_2) / M_PI) as isize;
            }
            else {
                t = phi.tan();
                mod_ = ((phi - t.atan())/M_PI).floor() as isize;
            }
            c = (a - b) / 2.0;
            let temp2 = (a * b).sqrt();
            a = (a + b) / 2.0;
            b = temp2;
            d += d;
        }

        (t.atan() + (mod_ as f64) * M_PI) / (d as f64 * a)
    };

    if sign < 0 {
        temp = -temp;
    }
    temp + npio2 * k
}

fn max3(a: f64, b: f64, c: f64) -> f64 {
    a.max(b.max(c))
}

/* To calculate legendre's incomplete elliptical integral of the first kind for
 * negative m, we use a power series in phi for small m*phi*phi, an asymptotic
 * series in m for large m*phi*phi* and the relation to Carlson's symmetric
 * integral of the first kind.
 *
 * F(phi, m) = sin(phi) * R_F(cos(phi)^2, 1 - m * sin(phi)^2, 1.0)
 *           = R_F(c-1, c-m, c)
 *
 * where c = csc(phi)^2. We use the second form of this for (approximately)
 * phi > 1/(sqrt(DBL_MAX) ~ 1e-154, where csc(phi)^2 overflows. Elsewhere we
 * use the first form, accounting for the smallness of phi.
 *
 * The algorithm used is described in Carlson, B. C. Numerical computation of
 * real or complex elliptic integrals. (1994) https://arxiv.org/abs/math/9409227
 * Most variable names reflect Carlson's usage.
 *
 * In this routine, we assume m < 0 and  0 > phi > pi/2.
 */
fn ellik_neg_m(phi: f64, m: f64) -> f64
{
    //double x, y, z, x1, y1, z1, a0, A, Q, X, Y, Z, e2, e3, scale;
    let mut n: isize = 0;
    let mpp = (m*phi)*phi;

    if -mpp < 1e-6 && phi < -m {
        return phi + (-mpp*phi*phi/30.0  + 3.0*mpp*mpp/40.0 + mpp/6.0)*phi;
    }

    if -mpp > 4e7 {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = (4.0 * sp * sm / (1.0 + cp)).ln();
        let b = -(1.0 + cp/sp/sp - a) / 4.0 / m;
        return (a + b) / sm;
    }

    let (scale, x, y, z) = if phi > 1e-153 && m > -1e305 {
        let s = phi.sin();
        let csc2 = 1.0 / (s*s);
        (1.0, 1.0 / (phi.tan() * phi.tan()), csc2 - m, csc2)
    }
    else {
        (phi, 1.0, 1.0 - m*phi*phi, 1.0)
    };

    if x == y && x == z {
        return scale / x.sqrt();
    }

    let a0 = (x + y + z) / 3.0;
    let mut a = a0;
    let mut x1 = x; 
    let mut y1 = y; 
    let mut z1 = z;
    /* Carlson gives 1/pow(3*r, 1.0/6.0) for this constant. if r == eps,
     * it is ~338.38. */
    let mut q = 400.0 * max3((a0-x).abs(), (a0-y).abs(), (a0-z).abs());

    while q > a.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx*sy + sx*sz + sy*sz;
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        a = (x1 + y1 + z1) / 3.0;
        n += 1;
        q /= 4.0;
    }
    let x = (a0 - x) / a / ((1 << (2 * n)) as f64);
    let y = (a0 - y) / a / ((1 << (2 * n)) as f64);
    let z = -(x + y);

    let e2 = x*y - z*z;
    let e3 = x*y*z;

    scale * (1.0 - e2/10.0 + e3/14.0 + e2*e2/24.0 - 3.0*e2*e3/44.0) / a.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ellik_trivials() {
        assert_eq!(ellik(f64::NAN, 0.5).is_nan(), true);
        assert_eq!(ellik(0.0, f64::NAN).is_nan(), true);
        assert_eq!(ellik(0.0, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(ellik(0.0, 10.0).is_nan(), true);
        assert_eq!(ellik(0.0, f64::INFINITY).is_nan(), true);
        assert_eq!(ellik(f64::INFINITY, -f64::INFINITY).is_nan(), true);
        assert_eq!(ellik(0.5, -f64::INFINITY), 0.0);
        assert_eq!(ellik(f64::INFINITY, 0.5), f64::INFINITY);
        assert_eq!(ellik(-f64::INFINITY, 0.5), -f64::INFINITY);
        assert_eq!(ellik(0.0, 0.0), 0.0);
        assert_eq!(ellik(10.0, 0.0), 10.0);
        assert_eq!(ellik(-10.0, 0.0), -10.0);
        assert_eq!(ellik(M_PI_2 + 1e-20, 1.0), f64::INFINITY);
        assert_eq!(ellik(-M_PI_2 - 1e-20, 1.0), f64::INFINITY);
    }

    #[test]
    fn ellik_m_1() {
        assert_eq!(ellik(0.0, 1.0), 0.0);

        assert_eq!(ellik(1e-10, 1.0), 1e-10);
        assert_eq!(ellik(0.1, 1.0), 0.10016708454748019);
        assert_eq!(ellik(1.0, 1.0), 1.2261911708835171);
        assert_eq!(ellik(M_PI_2 - 0.1, 1.0), 2.9948984537675716);
        assert_eq!(ellik(M_PI_2 - 1e-10, 1.0), 23.718997415436874);
        assert_eq!(ellik(M_PI_2 - 1e-15, 1.0), 35.07367669831271);

        assert_eq!(ellik(-1e-10, 1.0), -1e-10);
        assert_eq!(ellik(-0.1, 1.0), -0.10016708454748019);
        assert_eq!(ellik(-1.0, 1.0), -1.2261911708835171);
        assert_eq!(ellik(-M_PI_2 + 0.1, 1.0), -2.9948984537675716);
        assert_eq!(ellik(-M_PI_2 + 1e-10, 1.0), -23.718997415436874);
        assert_eq!(ellik(-M_PI_2 + 1e-15, 1.0), -35.07367669831271);
    }

    #[test]
    fn ellik_m_neg() {
        assert_eq!(ellik(0.0, -0.1), 0.0);

        assert_eq!(ellik(1e-10, -0.1), 1e-10);
        assert_eq!(ellik(0.1, -0.1), 0.0999833740948805);
        assert_eq!(ellik(1.0, -0.1), 0.9868118487668738);
        assert_eq!(ellik(10.0, -0.1), 9.773849855431353);
        assert_eq!(ellik(1e10, -0.1), 9763155117.895386);
        assert_eq!(ellik(1e20, -0.1), 9.763155117905381e+19);

        assert_eq!(ellik(-1e-10, -0.1), -1e-10);
        assert_eq!(ellik(-0.1, -0.1), -0.0999833740948805);
        assert_eq!(ellik(-1.0, -0.1), -0.9868118487668738);
        assert_eq!(ellik(-10.0, -0.1), -9.773849855431353);
        assert_eq!(ellik(-1e10, -0.1), -9763155117.895386);
        assert_eq!(ellik(-1e20, -0.1), -9.763155117905381e+19);

        assert_eq!(ellik(1e-10, -1.0), 1e-10);
        assert_eq!(ellik(0.1, -1.0), 0.09983440838641279);
        assert_eq!(ellik(1.0, -1.0), 0.8963937894628946);
        assert_eq!(ellik(10.0, -1.0), 8.415142198945016);
        assert_eq!(ellik(1e10, -1.0), 8346268416.675423);
        assert_eq!(ellik(1e20, -1.0), 8.346268416740732e+19);

        assert_eq!(ellik(-1e-10, -1.0), -1e-10);
        assert_eq!(ellik(-0.1, -1.0), -0.09983440838641279);
        assert_eq!(ellik(-1.0, -1.0), -0.8963937894628946);
        assert_eq!(ellik(-10.0, -1.0), -8.415142198945016);
        assert_eq!(ellik(-1e10, -1.0), -8346268416.675423);
        assert_eq!(ellik(-1e20, -1.0), -8.346268416740732e+19);
    }

    #[test]
    fn ellik_nominal() {
        assert_eq!(ellik(M_PI, 1e-20), 3.141592653589793);
        assert_eq!(ellik(100.0, 1e-10), 100.00000000251092);
        assert_eq!(ellik(-10.0, 1e-5), -10.000023858951376);
        assert_eq!(ellik(0.1351, 0.1), 0.1351409811799065);
        assert_eq!(ellik(-0.1351, 0.1), -0.1351409811799065);
        assert_eq!(ellik(-M_PI, 0.5), -3.7081493546027438);
        assert_eq!(ellik(-100.0, 0.9), -164.44309769019648);
        assert_eq!(ellik(30.0 * M_PI + 0.1, 1.0 - 1e-10), 774.0533541855924);
        assert_eq!(ellik(-30.0 * M_PI + 0.1, 1.0 - 1e-15), -1119.2647743528726);

        assert_eq!(ellik(-M_PI_2, 0.9), -2.5780921133481733);
        assert_eq!(ellik(-3.0 * M_PI_2 - 0.05, 1e-10), -4.762388980504998);
        assert_eq!(ellik(31.0 * M_PI_2, 1e-20), 48.69468613064179);
        assert_eq!(ellik(15.0 * M_PI_2 + 0.05, 1.0 - 1e-10), 202.69884547978893);
        assert_eq!(ellik(-15.0 * M_PI_2, 1.0 - 1e-16), -296.3204195149779);
        assert_eq!(ellik(M_PI_2, 0.1), 1.6124413487202192);
        assert_eq!(ellik(M_PI_2 - 0.01, 0.1), 1.601900442706069);
        assert_eq!(ellik(M_PI_2 + 0.01, 0.1), 1.6229822547343693);
    }
}