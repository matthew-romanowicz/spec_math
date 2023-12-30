/*
* (C) Copyright John Maddock 2006.
* Use, modification and distribution are subject to the
* Boost Software License, Version 1.0. (See accompanying file
*  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)
*/

#![allow(clippy::excessive_precision)]

use crate::cephes64::consts::EULER;
use crate::cephes64::polevl::polevl;
use crate::cephes64::unity::log1p;
use crate::cephes64::gamma::{gamma, lgam};
use crate::cephes64::igam::{igam, igamc, igam_fac};

const A: [f64; 4] = [0.213623493715853, 4.28342155967104,
        11.6616720288968, 3.31125922108741];
const B: [f64; 5] = [0.3611708101884203e-1, 1.27364489782223,
        6.40691597760039, 6.61053765625462, 1.0];

fn find_inverse_s(p: f64, q: f64) -> f64
{
    /*
    * Computation of the Incomplete Gamma Function Ratios and their Inverse
    * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
    * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
    * December 1986, Pages 377-393.
    *
    * See equation 32.
    */

    let t = if p < 0.5 {
        (-2.0 * p.ln()).sqrt()
    }
    else {
        (-2.0 * q.ln()).sqrt()
    };
    let s = t - polevl(t, &A, 3) / polevl(t, &B, 4);
    if p < 0.5 {
        -s
    } else {
        s
    }
}


fn didonato_sn(a: f64, x: f64, n: usize, tolerance: f64) -> f64
{
    /*
    * Computation of the Incomplete Gamma Function Ratios and their Inverse
    * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
    * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
    * December 1986, Pages 377-393.
    *
    * See equation 34.
    */
    let mut sum: f64 = 1.0;

    if n >= 1 {
        let mut partial: f64 = x / (a + 1.0);

        sum += partial;
        for i in 2..=n {
            partial *= x / (a + i as f64);
            sum += partial;
            if partial < tolerance {
                break;
            }
        }
    }
    
    sum
}


fn find_inverse_gamma(a: f64, p: f64, q: f64) -> f64
{
    /*
    * In order to understand what's going on here, you will
    * need to refer to:
    *
    * Computation of the Incomplete Gamma Function Ratios and their Inverse
    * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
    * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
    * December 1986, Pages 377-393.
    */

    if a == 1.0 {
        if q > 0.9 {
            -log1p(-p)
        }
        else {
            -q.ln()
        }
    } else if a < 1.0 {
        let g: f64 = gamma(a);
        let b: f64 = q * g;

        if (b > 0.6) || ((b >= 0.45) && (a >= 0.3)) {
            /* DiDonato & Morris Eq 21:
            *
            * There is a slight variation from DiDonato and Morris here:
            * the first form given here is unstable when p is close to 1,
            * making it impossible to compute the inverse of Q(a,x) for small
            * q. Fortunately the second form works perfectly well in this case.
            */
            let u: f64 = if (b * q > 1e-8) && (q > 1e-5) {
                (p * g * a).powf(1.0 / a)
            } else {
                ((-q / a) - EULER).exp()
            };
            u / (1.0 - (u / (a + 1.0)))
        } else if (a < 0.3) && (b >= 0.35) {
            /* DiDonato & Morris Eq 22: */
            let t = (-EULER - b).exp();
            let u = t * t.exp();
            t * u.exp()
        } else if (b > 0.15) || (a >= 0.3) {
            /* DiDonato & Morris Eq 23: */
            let y = -b.ln();
            let  u = y - (1.0 - a) * y.ln();
            y - (1.0 - a) * u.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln()
        } else if b > 0.1 { // TODO: This case may be impossible to reach
            /* DiDonato & Morris Eq 24: */
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            y - (1.0 - a) * u.ln() - ((u * u + 2.0 * (3.0 - a) * u + (2.0 - a) * (3.0 - a)) / (u * u + (5.0 - a) * u + 2.0)).ln()
        } else {
            /* DiDonato & Morris Eq 25: */
            let y = -b.ln();
            let c1 = (a - 1.0) * y.ln();
            let c1_2 = c1 * c1;
            let c1_3 = c1_2 * c1;
            let c1_4 = c1_2 * c1_2;
            let a_2 = a * a;
            let a_3 = a_2 * a;

            let c2 = (a - 1.0) * (1.0 + c1);
            let c3 = (a - 1.0) * (-(c1_2 / 2.0)
                                + (a - 2.0) * c1
                                + (3.0 * a - 5.0) / 2.0);
            let c4 = (a - 1.0) * ((c1_3 / 3.0) - (3.0 * a - 5.0) * c1_2 / 2.0
                                + (a_2 - 6.0 * a + 7.0) * c1
                                + (11.0 * a_2 - 46.0 * a + 47.0) / 6.0);
            let c5 = (a - 1.0) * (-(c1_4 / 4.0)
                                + (11.0 * a - 17.0) * c1_3 / 6.0
                                + (-3.0 * a_2 + 13.0 * a -13.0) * c1_2
                                + (2.0 * a_3 - 25.0 * a_2 + 72.0 * a - 61.0) * c1 / 2.0
                                + (25.0 * a_3 - 195.0 * a_2 + 477.0 * a - 379.0) / 12.0);

            let y_2 = y * y;
            let y_3 = y_2 * y;
            let y_4 = y_2 * y_2;
            y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4)
        }
    } else {
        /* DiDonato and Morris Eq 31: */
        let s = find_inverse_s(p, q);

        let s_2 = s * s;
        let s_3 = s_2 * s;
        let s_4 = s_2 * s_2;
        let s_5 = s_4 * s;
        let ra = a.sqrt();

        let mut w = a + s * ra + (s_2 - 1.0) / 3.0;
        w += (s_3 - 7.0 * s) / (36.0 * ra);
        w -= (3.0 * s_4 + 7.0 * s_2 - 16.0) / (810.0 * a);
        w += (9.0 * s_5 + 256.0 * s_3 - 433.0 * s) / (38880.0 * a * ra);

        if (a >= 500.0) && ((1.0 - w / a).abs() < 1e-6) {
            w
        } else if p > 0.5 {
            if w < 3.0 * a {
                w
            } else {
                let d = 2.0_f64.max(a * (a - 1.0));
                let lg = lgam(a);
                let lb = q.ln() + lg;
                if lb < -d * 2.3 {
                    /* DiDonato and Morris Eq 25: */
                    let y = -lb;
                    let c1 = (a - 1.0) * y.ln();
                    let c1_2 = c1 * c1;
                    let c1_3 = c1_2 * c1;
                    let c1_4 = c1_2 * c1_2;
                    let a_2 = a * a;
                    let a_3 = a_2 * a;

                    let c2 = (a - 1.0) * (1.0 + c1);
                    let c3 = (a - 1.0) * (-(c1_2 / 2.0)
                        + (a - 2.0) * c1
                        + (3.0 * a - 5.0) / 2.0);
                    let c4 = (a - 1.0) * ((c1_3 / 3.0)
                        - (3.0 * a - 5.0) * c1_2 / 2.0
                        + (a_2 - 6.0 * a + 7.0) * c1
                        + (11.0 * a_2 - 46.0 * a + 47.0) / 6.0);
                    let c5 = (a - 1.0) * (-(c1_4 / 4.0)
                        + (11.0 * a - 17.0) * c1_3 / 6.0
                        + (-3.0 * a_2 + 13.0 * a -13.0) * c1_2
                        + (2.0 * a_3 - 25.0 * a_2 + 72.0 * a - 61.0) * c1 / 2.0
                        + (25.0 * a_3 - 195.0 * a_2 + 477.0 * a - 379.0) / 12.0);

                    let y_2 = y * y;
                    let y_3 = y_2 * y;
                    let y_4 = y_2 * y_2;
                    y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4)
                } else {
                    /* DiDonato and Morris Eq 33: */
                    let u = -lb + (a - 1.0) * w.ln() - (1.0 + (1.0 - a) / (1.0 + w)).ln();
                    -lb + (a - 1.0) * u.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln()
                }
            }
        } else {
            let mut z = w;
            let ap1 = a + 1.0;
            let ap2 = a + 2.0;
            if w < 0.15 * ap1 {
                /* DiDonato and Morris Eq 35: */
                let v = p.ln() + lgam(ap1);
                z = ((v + w) / a).exp();
                let mut s = log1p(z / ap1 * (1.0 + z / ap2));
                z = ((v + z - s) / a).exp();
                s = log1p(z / ap1 * (1.0 + z / ap2));
                z = ((v + z - s) / a).exp();
                s = log1p(z / ap1 * (1.0 + z / ap2 * (1.0 + z / (a + 3.0))));
                z = ((v + z - s) / a).exp();
            }

            if (z <= 0.01 * ap1) || (z > 0.7 * ap1) {
                z
            } else {
                /* DiDonato and Morris Eq 36: */
                let ls = didonato_sn(a, z, 100, 1e-4).ln();
                let v = p.ln() + lgam(ap1);
                let z = ((v + z - ls) / a).exp();
                z * (1.0 - (a * z.ln() - z - v + ls) / (a - z))
            }
        }
    }
}


pub fn igami(a: f64, p: f64) -> f64 {

    if a.is_nan() || p.is_nan() || (a < 0.0) || !(0.0..=1.0).contains(&p) {
        //sf_error("gammaincinv", SF_ERROR_DOMAIN, NULL);
        return f64::NAN;
    } else if p == 0.0 {
        return 0.0;
    } else if p == 1.0 {
        return f64::INFINITY;
    } else if p > 0.9 {
        return igamci(a, 1.0 - p);
    }

    let mut x = find_inverse_gamma(a, p, 1.0 - p);
    /* Halley's method */
    for _ in 0..3 {
        let fac = igam_fac(a, x);
        if fac == 0.0 {
            break;
        }
        let f_fp = (igam(a, x) - p) * x / fac;
        /* The ratio of the first and second derivatives simplifies */
        let fpp_fp = -1.0 + (a - 1.0) / x;
        if fpp_fp.is_infinite() {
            /* Resort to Newton's method in the case of overflow */
            x -= f_fp;
        } else {
            x -= f_fp / (1.0 - 0.5 * f_fp * fpp_fp);
        }
    }

    x
}


pub fn igamci(a: f64, q: f64) -> f64 {
    // int i;
    // double x, fac, f_fp, fpp_fp;

    if a.is_nan() || q.is_nan() {
        return f64::NAN;
    }
    else if (a < 0.0) || !(0.0..=1.0).contains(&q) {
        //sf_error("gammainccinv", SF_ERROR_DOMAIN, NULL);
    }
    else if q == 0.0 {
        return f64::INFINITY;
    }
    else if q == 1.0 {
        return 0.0;
    }
    else if q > 0.9 {
        return igami(a, 1.0 - q);
    }

    let mut x = find_inverse_gamma(a, 1.0 - q, q);
    for _ in 0..3 {
        let fac = igam_fac(a, x);
        if fac == 0.0 {
            break;
        }
        let f_fp = -(igamc(a, x) - q) * x / fac;
        let fpp_fp = -1.0 + (a - 1.0) / x;
        if fpp_fp.is_infinite() {
            x -= f_fp;
        } else {
            x -= f_fp / (1.0 - 0.5 * f_fp * fpp_fp);
        }
    }

    x
}

#[cfg(test)]
mod igami_tests {
    use super::*;

    #[test]
    fn igami_trivials() {
        assert_eq!(igami(f64::NAN, 0.5).is_nan(), true);
        assert_eq!(igami(0.5, f64::NAN).is_nan(), true);
        assert_eq!(igami(0.5, -0.1).is_nan(), true);
        assert_eq!(igami(-0.1, 0.5).is_nan(), true);
        assert_eq!(igami(0.5, 1.1).is_nan(), true);
        assert_eq!(igami(0.5, 0.0), 0.0);
        assert_eq!(igami(0.5, 1.0), f64::INFINITY);
    }

    #[test]
    fn igami_large_p() {
        assert_eq!(igami(0.1, 0.91), 0.30723850581594636);
        assert_eq!(igami(10.0, 0.91), 14.44370657056981);
        assert_eq!(igami(1e5, 0.91), 100424.24923098792);
        assert_eq!(igami(1e10, 0.91), 10000134075.769241);
        assert_eq!(igami(1e20, 0.91), 1.0000000001340755e+20);

        assert_eq!(igami(0.2, 0.95), 1.0305257996426673);
        assert_eq!(igami(20.0, 0.95), 27.87923963944351);
        assert_eq!(igami(2e5, 0.95), 200736.1689801604);
        assert_eq!(igami(2e10, 0.95), 20000232617.99925);
        assert_eq!(igami(2e20, 0.95), 2.0000000002326174e+20);

        assert_eq!(igami(0.5, 0.99999), 9.755710482333134);
        assert_eq!(igami(50.0, 0.99999), 86.04947101807929);
        assert_eq!(igami(5e5, 0.99999), 503021.4648378374);
        assert_eq!(igami(5e10, 0.99999), 50000953664.302956);
        assert_eq!(igami(5e20, 0.99999), 5.000000000953659e+20);
    }

    #[test]
    fn igami_a_1() {
        assert_eq!(igami(1.0, 1e-10), 1.0000000000499969e-10);
        assert_eq!(igami(1.0, 0.05), 0.05129329438755053);
        assert_eq!(igami(1.0, 0.5), 0.6931471805599455);
        assert_eq!(igami(1.0, 0.95), 2.9957322735539895);
        assert_eq!(igami(1.0, 1.0 - 1e-10), 23.02585084720009);
    }

    #[test]
    fn igami_a_small() {
        assert_eq!(igami(1e-1, 1e-10), 6.073048362408172e-101);
        assert_eq!(igami(1e-2, 0.05), 4.465535018910624e-131);
        assert_eq!(igami(1e-3, 0.5), 5.244206408274974e-302);
        assert_eq!(igami(1e-4, 0.95), 9.669567020847056e-224);
        assert_eq!(igami(1e-5, 1.0 - 1e-10), 9.198940494977181);

        assert_eq!(igami(0.2, 1e-10), 6.525480843457197e-51);
        assert_eq!(igami(0.2, 0.05), 2.0392131101128526e-07);
        assert_eq!(igami(0.2, 0.5), 0.02074633919282486);
        assert_eq!(igami(0.2, 0.95), 1.0305257996426673);
        assert_eq!(igami(0.2, 1.0 - 1e-10), 19.102745525758326);

        assert_eq!(igami(0.5, 1e-10), 7.85398163397448e-21);
        assert_eq!(igami(0.5, 0.05), 0.001966070000009761);
        assert_eq!(igami(0.5, 0.5), 0.227468211559786);
        assert_eq!(igami(0.5, 0.95), 1.920729410347062);
        assert_eq!(igami(0.5, 1.0 - 1e-10), 20.910728101491394);

        assert_eq!(igami(0.8, 1e-10), 2.8934156942517426e-13);
        assert_eq!(igami(0.8, 0.05), 0.02189753286142773);
        assert_eq!(igami(0.8, 0.5), 0.5013512263630726);
        assert_eq!(igami(0.8, 0.95), 2.5951436134799155);
        assert_eq!(igami(0.8, 1.0 - 1e-10), 22.244784919035467);

        assert_eq!(igami(1.0 - 1e-5, 1e-10), 9.997655388750058e-11);
        assert_eq!(igami(1.0 - 1e-5, 0.05), 0.05129150185635944);
        assert_eq!(igami(1.0 - 1e-5, 0.5), 0.6931375001152216);
        assert_eq!(igami(1.0 - 1e-5, 0.95), 2.9957129055994884);
        assert_eq!(igami(1.0 - 1e-5, 1.0 - 1e-10), 23.025813291912808);

        assert_eq!(igami(0.2999, 0.8), 0.4598827808614538);
        assert_eq!(igami(0.8, 0.88), 1.7760548371325577);
    }

    #[test]
    fn igami_a_large() {
        assert_eq!(igami(1.0 + 1e-5, 1e-10), 1.000234511583614e-10);
        assert_eq!(igami(1.0 + 1e-5, 0.05), 0.05129508695002124);
        assert_eq!(igami(1.0 + 1e-5, 0.5), 0.6931568610118308);
        assert_eq!(igami(1.0 + 1e-5, 0.95), 2.995751641450343);
        assert_eq!(igami(1.0 + 1e-5, 1.0 - 1e-10), 23.02588840235435);

        assert_eq!(igami(10.0, 1e-10), 0.47272209260635223);
        assert_eq!(igami(10.0, 0.05), 5.425405697091295);
        assert_eq!(igami(10.0, 0.5), 9.668714614714128);
        assert_eq!(igami(10.0, 0.95), 15.705216422115459);
        assert_eq!(igami(10.0, 1.0 - 1e-10), 44.627857114101495);

        assert_eq!(igami(1e5, 1e-10), 98001.50416544829);
        assert_eq!(igami(1e5, 0.05), 99480.42074678872);
        assert_eq!(igami(1e5, 0.5), 99999.6666668642);
        assert_eq!(igami(1e5, 0.95), 100520.7162815659);
        assert_eq!(igami(1e5, 1.0 - 1e-10), 102024.80680796504);

        assert_eq!(igami(1e15, 1e-10), 999999805579286.9);
        assert_eq!(igami(1e15, 0.05), 999999947985161.8);
        assert_eq!(igami(1e15, 0.5), 999999999999999.6);
        assert_eq!(igami(1e15, 0.95), 1000000052014839.4);
        assert_eq!(igami(1e15, 1.0 - 1e-10), 1000000201163275.0);
    }
}

#[cfg(test)]
mod igamci_tests {
    use super::*;

    #[test]
    fn igamci_trivials() {
        assert_eq!(igamci(f64::NAN, 0.5).is_nan(), true);
        assert_eq!(igamci(0.5, f64::NAN).is_nan(), true);
        assert_eq!(igamci(0.5, -0.1).is_nan(), true);
        assert_eq!(igamci(-0.1, 0.5).is_nan(), true);
        assert_eq!(igamci(0.5, 1.1).is_nan(), true);
        assert_eq!(igamci(0.5, 1.0), 0.0);
        assert_eq!(igamci(0.5, 0.0), f64::INFINITY);
    }

    #[test]
    fn igamci_large_q() {
        assert_eq!(igamci(0.1, 0.91), 2.117541029697004e-11);
        assert_eq!(igamci(10.0, 0.91), 6.086398538483476);
        assert_eq!(igamci(1e5, 0.91), 99576.28251823064);
        assert_eq!(igamci(1e10, 0.91), 9999865924.762508);
        assert_eq!(igamci(1e20, 0.91), 9.999999998659245e+19);

        assert_eq!(igamci(0.2, 0.95), 2.0392131101128603e-07);
        assert_eq!(igamci(20.0, 0.95), 13.254651598346555);
        assert_eq!(igamci(2e5, 0.95), 199264.96804850159);
        assert_eq!(igamci(2e10, 0.95), 19999767383.13778);
        assert_eq!(igamci(2e20, 0.95), 1.9999999997673826e+20);

        assert_eq!(igamci(0.5, 0.99999), 7.853981634314245e-11);
        assert_eq!(igamci(50.0, 0.99999), 25.355468772812532);
        assert_eq!(igamci(5e5, 0.99999), 496989.9946857008);
        assert_eq!(igamci(5e10, 0.99999), 49999046347.15658);
        assert_eq!(igamci(5e20, 0.99999), 4.999999999046341e+20);
    }

    #[test]
    fn igamci_a_1() {
        assert_eq!(igamci(1.0, 1e-10), 23.025850929940457);
        assert_eq!(igamci(1.0, 0.05), 2.9957322735539913);
        assert_eq!(igamci(1.0, 0.5), 0.6931471805599455);
        assert_eq!(igamci(1.0, 0.95), 0.05129329438755058);
        assert_eq!(igamci(1.0, 1.0 - 1e-10), 1.0000000827903717e-10);
    }

    #[test]
    fn igamci_a_small() {
        assert_eq!(igamci(1e-5, 1e-10), 9.198940570224709);
        assert_eq!(igamci(1e-4, 0.05), 9.669567020851293e-224);
        assert_eq!(igamci(1e-3, 0.5), 5.244206408274974e-302);
        assert_eq!(igamci(1e-2, 0.95), 4.4655350189108117e-131);
        assert_eq!(igamci(1e-1, 1.0 - 1e-10), 6.073053387272597e-101);

        assert_eq!(igamci(0.2, 1e-10), 19.102745605319942);
        assert_eq!(igamci(0.2, 0.05), 1.0305257996426687);
        assert_eq!(igamci(0.2, 0.5), 0.02074633919282486);
        assert_eq!(igamci(0.2, 0.95), 2.0392131101128603e-07);
        assert_eq!(igamci(0.2, 1.0 - 1e-10), 6.525483543061197e-51);

        assert_eq!(igamci(0.5, 1e-10), 20.910728182380645);
        assert_eq!(igamci(0.5, 0.05), 1.9207294103470642);
        assert_eq!(igamci(0.5, 0.5), 0.2274682115597862);
        assert_eq!(igamci(0.5, 0.95), 0.0019660700000097655);
        assert_eq!(igamci(0.5, 1.0 - 1e-10), 7.853982933657245e-21);

        assert_eq!(igamci(0.8, 1e-10), 22.244785001068575);
        assert_eq!(igamci(0.8, 0.05), 2.5951436134799146);
        assert_eq!(igamci(0.8, 0.5), 0.5013512263630726);
        assert_eq!(igamci(0.8, 0.95), 0.021897532861427763);
        assert_eq!(igamci(0.8, 1.0 - 1e-10), 2.8934159935046064e-13);

        assert_eq!(igamci(1.0 - 1e-5, 1e-10), 23.025813374653136);
        assert_eq!(igamci(1.0 - 1e-5, 0.05), 2.995712905599489);
        assert_eq!(igamci(1.0 - 1e-5, 0.5), 0.6931375001152216);
        assert_eq!(igamci(1.0 - 1e-5, 0.95), 0.05129150185635947);
        assert_eq!(igamci(1.0 - 1e-5, 1.0 - 1e-10), 9.997656215968052e-11);

        assert_eq!(igamci(0.2999, 0.2), 0.4598827808614532);
        assert_eq!(igamci(0.8, 0.12), 1.7760548371325586);
    }

    #[test]
    fn igamci_a_large() {
        assert_eq!(igamci(1.0 + 1e-5, 1e-10), 23.025888485094754);
        assert_eq!(igamci(1.0 + 1e-5, 0.05), 2.9957516414503433);
        assert_eq!(igamci(1.0 + 1e-5, 0.5), 0.6931568610118308);
        assert_eq!(igamci(1.0 + 1e-5, 0.95), 0.05129508695002131);
        assert_eq!(igamci(1.0 + 1e-5, 1.0 - 1e-10), 1.0002345943425605e-10);

        assert_eq!(igamci(10.0, 1e-10), 44.627857217059066);
        assert_eq!(igamci(10.0, 0.05), 15.70521642211546);
        assert_eq!(igamci(10.0, 0.5), 9.668714614714128);
        assert_eq!(igamci(10.0, 0.95), 5.425405697091293);
        assert_eq!(igamci(10.0, 1.0 - 1e-10), 0.472722096692631);

        assert_eq!(igamci(1e5, 1e-10), 102024.80681203725);
        assert_eq!(igamci(1e5, 0.05), 100520.7162815659);
        assert_eq!(igamci(1e5, 0.5), 99999.6666668642);
        assert_eq!(igamci(1e5, 0.95), 99480.42074678872);
        assert_eq!(igamci(1e5, 1.0 - 1e-10), 98001.50416941273);

        assert_eq!(igamci(1e15, 1e-10), 1000000201163275.4);
        assert_eq!(igamci(1e15, 0.05), 1000000052014839.4);
        assert_eq!(igamci(1e15, 0.5), 999999999999999.6);
        assert_eq!(igamci(1e15, 0.95), 999999947985161.8);
        assert_eq!(igamci(1e15, 1.0 - 1e-10), 999999805579287.2);
    }
}