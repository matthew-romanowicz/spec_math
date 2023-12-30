/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1984, 1987, 1993 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
/* Copyright 2014, Eric W. Moore */
 

#![allow(clippy::excessive_precision)]

use super::consts::{M_PI, M_PI_2, MACHEP};
use super::unity::cosm1;
use super::ellpk::ellpk;
use super::ellpe::ellpe;

 // from https://github.com/scipy/scipy/blob/c4ce0c4560bc635867512c4d2ea6db6f666d3eeb/scipy/special/cephes/unity.c#L144


//$$E(\varphi, m) = \int_{0}^{\varphi} \sqrt{1 - m\,\sin^2(\theta)} \,d\theta$$

pub fn ellie(phi: f64, m: f64) -> f64 {
    //! Incomplete elliptic integral of the second kind
    //!
    //! ## Description
    //! Approximates the integral 
    //!
    #![doc=include_str!("ellie.svg")]
    //!
    //! of amplitude `phi` and modulus `m`, using the arithmetic - geometric mean algorithm.
    //!
    //! Returns `NAN` if either input is `NAN` or if `m > 1.0`.
    //!
    //! Returns `phi` if `phi` is infinite.
    //!
    //! Returns `-m` if `m` is infinte.
    //!
    //! ## Accuracy
    //!
    //! Tested at random arguments with `phi` in `[-10, 10]` and `m` in `[0, 1]`.
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
    //!     <td>150000</td>
    //!     <td>3.3e-15</td>
    //!     <td>1.4e-16</td>
    //! </tr>
    //!</table>
 
    if phi.is_nan() || m.is_nan() || m > 1.0 {
        return f64::NAN;
    } else if phi.is_infinite() {
        return phi;
    } else if m.is_infinite() {
        return -m;
    } else if m == 0.0 {
        return phi;
    }
    let mut lphi = phi;
    let mut npio2 = (lphi / M_PI_2).floor();
    if npio2.abs() % 2.0 == 1.0 {
        npio2 += 1.0;
    }
    lphi -= npio2 * M_PI_2;
    let sign = if lphi < 0.0 {
        lphi = -lphi;
        -1
    } else {
        1
    };
    let a = 1.0 - m;
    let ellpe_m = ellpe(m);
    let mut temp = if a == 0.0 {
        lphi.sin()
    } else if a > 1.0 {
        ellie_neg_m(lphi, m)
    } else if lphi < 0.135 {
        let m11= (((((-7.0/2816.0)*m + (5.0/1056.0))*m - (7.0/2640.0))*m
                    + (17.0/41580.0))*m - (1.0/155925.0))*m;
        let m9 = ((((-5.0/1152.0)*m + (1.0/144.0))*m - (1.0/360.0))*m
                    + (1.0/5670.0))*m;
        let m7 = ((-m/112.0 + (1.0/84.0))*m - (1.0/315.0))*m;
        let m5 = (-m/40.0 + (1.0 / 30.0))*m;
        let m3 = -m/6.0;
        let p2 = lphi * lphi;

        ((((m11*p2 + m9)*p2 + m7)*p2 + m5)*p2 + m3)*p2*lphi + lphi
    } else {
        let mut t = lphi.tan();
        let mut b = a.sqrt();
        /* Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
        * for pointing out an instability near odd multiples of pi/2.  */
        if t.abs() > 10.0 {
            /* Transform the amplitude */
            let mut e = 1.0 / (b * t);
            /* ... but avoid multiple recursions.  */
            if e.abs() < 10.0 {
                e = e.atan();
                let mut temp2 = ellpe_m + m * lphi.sin() * e.sin() - ellie(e, m);
                if sign < 0 {
                    temp2 = -temp2;
                }
                temp2 += npio2 * ellpe_m;
                return temp2;
            }
        }

        let mut c = m.sqrt();
        let mut a = 1.0;
        let mut d: isize = 1;
        let mut e = 0.0;
        let mut mod_: isize = 0;
    
        while (c / a).abs() > MACHEP {
            let temp2 = b / a;
            lphi = lphi + (t * temp2).atan() + (mod_ as f64) * M_PI;
            let denom = 1. - temp2 * t * t;
            if denom.abs() > 10.0 * MACHEP {
                t = t * (1.0 + temp2) / denom;
                mod_ = ((lphi + M_PI_2) / M_PI) as isize;
            }
            else {
                t = lphi.tan();
                mod_ = ((lphi - t.atan())/M_PI).floor() as isize;
            }
            c = (a - b) / 2.0;
            let temp2 = (a * b).sqrt();
            a = (a + b) / 2.0;
            b = temp2;
            d += d;
            e += c * (lphi).sin();
        }
    
        let mut temp2 = ellpe_m / ellpk(1.0 - m);
        temp2 *= (t.atan() + (mod_ as f64) * M_PI) / ((d as f64) * a);
        temp2 += e;
        temp2
    };

    if sign < 0 {
       temp = -temp;
    }
    temp += npio2 * ellpe_m;
    temp
 }
 
fn max3(a: f64, b: f64, c: f64) -> f64 {
   a.max(b.max(c))
}
 
 /* To calculate legendre's incomplete elliptical integral of the second kind for
  * negative m, we use a power series in phi for small m*phi*phi, an asymptotic
  * series in m for large m*phi*phi* and the relation to Carlson's symmetric
  * integrals, R_F(x,y,z) and R_D(x,y,z).
  * 
  * E(phi, m) = sin(phi) * R_F(cos(phi)^2, 1 - m * sin(phi)^2, 1.0)
  *             - m * sin(phi)^3 * R_D(cos(phi)^2, 1 - m * sin(phi)^2, 1.0) / 3
  *             
  *           = R_F(c-1, c-m, c) - m * R_D(c-1, c-m, c) / 3
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
fn ellie_neg_m(phi: f64, m: f64) -> f64 {
    //double x, y, z, x1, y1, z1, ret, Q;
    //double a0f, af, xf, yf, zf, e2f, e3f, scalef;
    //double a0d, ad, seriesn, seriesd, xd, yd, zd, e2d, e3d, e4d, e5d, scaled;
    let mut n: isize = 0;
    let mpp: f64 = (m * phi) * phi;
    
    if -mpp < 1e-6 && phi < -m {
        return phi + (mpp * phi * phi / 30.0 - mpp * mpp / 40.0 - mpp / 6.0) * phi;
    }
 
    if -mpp > 1e6 {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();
 
        let a = -cosm1(phi);
        let b1 = (4.0 * sp * sm / (1.0 + cp)).ln();
        let b = -(0.5 + b1) / 2.0 / m;
        let c = (0.75 + cp/sp/sp - b1) / 16.0 / m / m;
        return (a + b + c) * sm;
    }
 
    let (scalef, scaled, x, y, z) = if phi > 1e-153 && m > -1e200 {
        let s = phi.sin();
        let csc2 = 1.0 / s / s;
        // scalef = 1.0;
        // scaled = m / 3.0;
        // x = 1.0 / phi.tan() / phi.tan();
        // y = csc2 - m;
        // z = csc2;
        (1.0, m / 3.0, 1.0 / phi.tan() / phi.tan(), csc2 - m, csc2)
    }
    else {
        // scalef = phi;
        // scaled = mpp * phi / 3.0;
        // x = 1.0;
        // y = 1.0 - mpp;
        // z = 1.0;
        (phi, mpp * phi / 3.0, 1.0, 1.0 - mpp, 1.0)
    };
     
    if x == y && x == z {
        return (scalef + scaled / x) / (x.sqrt());
    }
 
    let a0f = (x + y + z) / 3.0;
    let mut af = a0f;
    let a0d = (x + y + 3.0*z) / 5.0;
    let mut ad = a0d;
    let mut x1 = x; 
    let mut y1 = y; 
    let mut z1 = z; 
    let mut seriesd = 0.0; 
    let mut seriesn = 1.0;
     /* Carlson gives 1/pow(3*r, 1.0/6.0) for this constant. if r == eps,
      * it is ~338.38. */
    let mut q = 400.0 * max3((a0f-x).abs(), (a0f-y).abs(), (a0f-z).abs());
     
    while q > af.abs() && q > ad.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx*sy + sx*sz + sy*sz;
        seriesd += seriesn / (sz * (z1 + lam));
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        af = (x1 + y1 + z1) / 3.0;
        ad = (ad + lam) / 4.0;
        n += 1;
        q /= 4.0;
        seriesn /= 4.0;
    }
 
    let xf = (a0f - x) / af / ((1 << (2 * n)) as f64);
    let yf = (a0f - y) / af / ((1 << (2 * n)) as f64);
    let zf = -(xf + yf);
 
    let e2f = xf*yf - zf*zf;
    let e3f = xf*yf*zf;
 
    let mut ret = scalef * (1.0 - e2f/10.0 + e3f/14.0 + e2f*e2f/24.0
                     - 3.0*e2f*e3f/44.0) / af.sqrt();
 
    let xd = (a0d - x) / ad / ((1 << (2 * n)) as f64);
    let yd = (a0d - y) / ad / ((1 << (2 * n)) as f64);
    let zd = -(xd + yd)/3.0;
 
    let e2d = xd*yd - 6.0*zd*zd;
    let e3d = (3.0 * xd * yd - 8.0 * zd * zd) * zd;
    let e4d = 3.0*(xd*yd - zd*zd)*zd*zd;
    let e5d = xd*yd*zd*zd*zd;
 
    ret -= scaled * (1.0 - 3.0*e2d/14.0 + e3d/6.0 + 9.0*e2d*e2d/88.0
                      - 3.0*e4d/22.0 - 9.0*e2d*e3d/52.0 + 3.0*e5d/26.0)
                      /((1 << (2 * n)) as f64) / ad / (ad.sqrt());
    ret -= 3.0 * scaled * seriesd;
    ret
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ellie_trivials() {
        assert_eq!(ellie(f64::NAN, 0.5).is_nan(), true);
        assert_eq!(ellie(0.0, f64::NAN).is_nan(), true);
        assert_eq!(ellie(0.0, 1.0 + 1e-10).is_nan(), true);
        assert_eq!(ellie(0.0, 10.0).is_nan(), true);
        assert_eq!(ellie(0.0, f64::INFINITY).is_nan(), true);
        assert_eq!(ellie(f64::INFINITY, 0.5), f64::INFINITY);
        assert_eq!(ellie(-f64::INFINITY, 0.5), -f64::INFINITY);
        assert_eq!(ellie(0.0, -f64::INFINITY), f64::INFINITY);
        assert_eq!(ellie(0.0, 0.0), 0.0);
        assert_eq!(ellie(10.0, 0.0), 10.0);
        assert_eq!(ellie(-10.0, 0.0), -10.0);
    }

    #[test]
    fn ellie_m_1() {
        assert_eq!(ellie(0.0, 1.0), 0.0);

        assert_eq!(ellie(1e-10, 1.0), 1e-10);
        assert_eq!(ellie(0.1, 1.0), 0.09983341664682815);
        assert_eq!(ellie(1.0, 1.0), 0.8414709848078965);
        assert_eq!(ellie(5.0, 1.0), 3.0410757253368614);
        assert_eq!(ellie(10.0, 1.0), 6.54402111088937);
        assert_eq!(ellie(100.0, 1.0), 63.49363435889025);
        assert_eq!(ellie(1e10, 1.0), 6366197723.512493);
        assert_eq!(ellie(1e20, 1.0), 6.3661977236758135e+19);

        assert_eq!(ellie(-1e-10, 1.0), -1e-10);
        assert_eq!(ellie(-0.1, 1.0), -0.09983341664682815);
        assert_eq!(ellie(-1.0, 1.0), -0.8414709848078965);
        assert_eq!(ellie(-5.0, 1.0), -3.0410757253368614);
        assert_eq!(ellie(-10.0, 1.0), -6.54402111088937);
        assert_eq!(ellie(-100.0, 1.0), -63.49363435889025);
        assert_eq!(ellie(-1e10, 1.0), -6366197723.512493);
        assert_eq!(ellie(-1e20, 1.0), -6.3661977236758135e+19);
    }

    #[test]
    fn ellie_m_neg() {
        assert_eq!(ellie(0.0, -0.1), 0.0);

        assert_eq!(ellie(1e-10, -0.1), 1e-10);
        assert_eq!(ellie(0.1, -0.1), 0.10001663087782428);
        assert_eq!(ellie(1.0, -0.1), 1.0134826965558748);
        assert_eq!(ellie(10.0, -0.1), 10.234331619985799);
        assert_eq!(ellie(1e10, -0.1), 10245497761.144846);
        assert_eq!(ellie(1e20, -0.1), 1.0245497761134428e+20);

        assert_eq!(ellie(-1e-10, -0.1), -1e-10);
        assert_eq!(ellie(-0.1, -0.1), -0.10001663087782428);
        assert_eq!(ellie(-1.0, -0.1), -1.0134826965558748);
        assert_eq!(ellie(-10.0, -0.1), -10.234331619985799);
        assert_eq!(ellie(-1e10, -0.1), -10245497761.144846);
        assert_eq!(ellie(-1e20, -0.1), -1.0245497761134428e+20);

        assert_eq!(ellie(1e-10, -1.0), 1e-10);
        assert_eq!(ellie(0.1, -1.0), 0.10016608571999015);
        assert_eq!(ellie(1.0, -1.0), 1.123887722945525);
        assert_eq!(ellie(10.0, -1.0), 12.064284079205486);
        assert_eq!(ellie(1e10, -1.0), 12160067234.3396);
        assert_eq!(ellie(1e20, -1.0), 1.2160067234249798e+20);

        assert_eq!(ellie(-1e-10, -1.0), -1e-10);
        assert_eq!(ellie(-0.1, -1.0), -0.10016608571999015);
        assert_eq!(ellie(-1.0, -1.0), -1.123887722945525);
        assert_eq!(ellie(-10.0, -1.0), -12.064284079205486);
        assert_eq!(ellie(-1e10, -1.0), -12160067234.3396);
        assert_eq!(ellie(-1e20, -1.0), -1.2160067234249798e+20);
    }

    #[test]
    fn ellie_small_lphi() {
        assert_eq!(ellie(0.1349, 0.1), 0.13485922238528766);
        assert_eq!(ellie(0.1349 + M_PI * 3.0, 0.1), 9.319405043771868);
        assert_eq!(ellie(-0.1349 - M_PI * 3.0, 0.1), -9.319405043771868);
        assert_eq!(ellie(-0.1 + M_PI * 15.0, 0.1), 45.82274574278697);
        assert_eq!(ellie(0.1 - M_PI * 55.0, 0.1), -168.28335669460805);
        assert_eq!(ellie(M_PI * 5.0, 0.1), 15.307576368977633);

        assert_eq!(ellie(-0.1349, 0.5), -0.13469588961866266);
        assert_eq!(ellie(-0.1349 + M_PI * 3.0, 0.5), 7.96916739666739);
        assert_eq!(ellie(0.1349 - M_PI * 3.0, 0.5), -7.96916739666739);
        assert_eq!(ellie(0.1 + M_PI * 55.0, 0.5), 148.6707436861049);
        assert_eq!(ellie(-0.1 - M_PI * 15.0, 0.5), -40.61923320229088);
        assert_eq!(ellie(-M_PI * 5.0, 0.5), -13.506438810476755);

        assert_eq!(ellie(-0.1349, 0.9), -0.13453219907093278);
        assert_eq!(ellie(0.1349 + M_PI * 1.0, 0.9), 2.3440816644790794);
        assert_eq!(ellie(-0.1349 - M_PI * 2.0, 0.9), -4.5536311298872265);
        assert_eq!(ellie(-0.1 + M_PI * 10.0, 0.9), 21.995644556553664);
        assert_eq!(ellie(0.1 - M_PI * 4.0, 0.9), -8.738347764104784);
        assert_eq!(ellie(-M_PI * 13.0, 0.9), -28.724143050305905);
    }

    #[test]
    fn ellie_nominal() {
        assert_eq!(ellie(M_PI, 1e-20), 3.141592653589793);
        assert_eq!(ellie(100.0, 1e-10), 99.99999999748908);
        assert_eq!(ellie(-10.0, 1e-5), -9.99997614113725);
        assert_eq!(ellie(0.1351, 0.1), 0.1350590411576287);
        assert_eq!(ellie(-0.1351, 0.1), -0.1350590411576287);
        assert_eq!(ellie(-M_PI, 0.5), -2.701287762095351);
        assert_eq!(ellie(-100.0, 0.9), -70.19665651542336);
        assert_eq!(ellie(30.0 * M_PI + 0.1, 1.0 - 1e-10), 60.09983345384451);
        assert_eq!(ellie(-30.0 * M_PI + 0.1, 1.0 - 1e-15), -59.90016658335372);

        // Falls into t.abs() > 10 case
        assert_eq!(ellie(-M_PI_2, 0.9), -1.1047747327040733);
        assert_eq!(ellie(-3.0 * M_PI_2 - 0.05, 1e-10), -4.762388980264382);
        assert_eq!(ellie(31.0 * M_PI_2, 1e-20), 48.69468613064179);
        assert_eq!(ellie(15.0 * M_PI_2 + 0.05, 1.0 - 1e-10), 15.001249749389917);
        assert_eq!(ellie(-15.0 * M_PI_2, 1.0 - 1e-20), -15.0);
        assert_eq!(ellie(M_PI_2, 0.1), 1.5307576368977633);
        assert_eq!(ellie(M_PI_2 - 0.01, 0.1), 1.5212707863494297);
        assert_eq!(ellie(M_PI_2 + 0.01, 0.1), 1.540244487446097);
    }
}
