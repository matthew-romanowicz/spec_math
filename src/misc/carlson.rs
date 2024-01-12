
const MAX_ITER: isize = 10000;

fn sum_carry(large: f64, small: f64, carry: &mut f64) -> f64 {
    let s = large + small;
    *carry += small - (s - large);
    s
}

fn q_sum(a0: f64, g0: f64, p0: f64, q0: f64) -> (f64, f64) {
    let mut a = a0;
    let mut g = g0;
    let mut p2 = p0 * p0;
    let mut p_next = p0;
    let mut q = q0;
    let mut s = q;
    let mut c: f64 = 0.0;
    for _ in 0..MAX_ITER {
        let a_next = 0.5 * (a + g);

        // Store g_next in it's squared form so it can be used in place
        // of a * g in other expressions. g will be assigned to 
        // g_next2.sqrt() at the end of the iteration.
        let g_next2 = a * g;

        let p2_ag = p2 + g_next2;

        p_next = p2_ag / (2.0 * p_next);
        let eps = (p2 - g_next2) / p2_ag;
        let q_next = 0.5 * q * eps;

        a = a_next;
        g = g_next2.sqrt();

        // Check that both s and a - g have converged
        if s + q_next == s && (a - g).abs() < 1e-15 { 
            c += q_next;
            return (s + c, 0.5 * (a + g));
        } else {
            s = sum_carry(s, q_next, &mut c);
            q = q_next;

            // p squared is calculated for the next iteration since
            // it's used multiple times. In the one case that p is
            // needed, p_next can be used because it is before the
            // reassignment of p_next.
            p2 = p_next * p_next;
        }
    }
    (s + c, 0.5 * (a + g))
}

pub fn ellippi(n: f64, m: f64) -> f64 {
    //! Complete elliptic integral of the third kind.
    //!
    //! # Description:
    //!
    //! Approximates the complete elliptic integral of the third kind
    //! using the arithmetic-geometric mean algorithm from
    //! [DLMF 19.8](https://dlmf.nist.gov/19.8#i.p3)
    //!
    //! # Domain:
    //!
    //! Function returns `NAN` if `n == 1.0` or `m >= 1.0`
    //!
    //! # Accuracy:
    //!
    //! Relative error:
    //!
    //!<table>
    //! <tr>
    //!     <th>Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>-100 <= n < 1.0
    //!
    //! -100 <= m < 1.0</td>
    //!     <td>30000</td>
    //!     <td>2.8e-15</td>
    //!     <td>2.5e-16</td>
    //! </tr>
    //!</table>
    if m >= 1.0 {
        return f64::NAN;
    }
    let k_p = (1.0 - m).sqrt();
    if n < 1.0 {
        let (q, agm) = q_sum(1.0, k_p, (1.0 - n).sqrt(), 1.0);
        let coef = std::f64::consts::FRAC_PI_4 / agm;
        coef * (2.0 + n / (1.0 - n) * q)
    } else if n > 1.0 {
        let (q, agm) = q_sum(1.0, k_p, (1.0 - m / n).sqrt(), 1.0);
        let coef = std::f64::consts::FRAC_PI_4 / agm;
        coef * (m / (m - n)) * q
    } else {
        f64::NAN
    }
}

pub fn elliprc(x: f64, y: f64) -> f64 {
    //! Degenerate Carlson symmetric elliptic integral of the first kind
    //!
    //! # Description
    //!
    //! Approximates Carlson's degnerate symmetric elliptric integral
    //! of the first kind (R_C(x, y)) using the technique from 
    //!
    //! B. C. Carlson, “Numerical computation of real or complex elliptic 
    //! integrals,” Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
    //! 
    //! https://arxiv.org/abs/math/9409227 https://doi.org/10.1007/BF02198293
    //!
    //! # Domain:
    //!
    //! Returns `NAN` if `x < 0.0` and infinity if `y == 0.0`.
    //!
    //! # Accuracy:
    //!
    //! Relative error:
    //!
    //!<table>
    //! <tr>
    //!     <th>Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>0 <= x < 100.0
    //!
    //! -100 <= y < 100.0</td>
    //!     <td>1999000</td>
    //!     <td>7.0e-16</td>
    //!     <td>1.4e-16</td>
    //! </tr>
    //!</table>

    // Peak: 7.00474240313977e-16
    // RMS: 1.385189339345323e-16
    // # Samples: 1999000

    if y < 0.0 && 0.0 <= x {
        if x == 0.0{
            // https://dlmf.nist.gov/19.6#E15
            return 0.0;
        } else {
            return (x / (x - y)).sqrt() * elliprc(x - y, -y);
        }
    } else if y == 0.0 {
        return f64::INFINITY;
    } else if x < 0.0 {
        return f64::NAN;
    }

    let mut x = x;
    let mut y = y;
    let a0 = (x + 2.0 * y) / 3.0;
    let mut a = a0;
    let r: f64 = 1e-16;
    let mut q = (3.0 * r).powf(-0.125) * (a - x).abs();
    let mut s = y - a0;
    loop {
        s *= 0.25;
        q *= 0.25;
        let lam = 2.0 * x.sqrt() * y.sqrt() + y;
        let an = (a + lam) * 0.25;
        let xn = (x + lam) * 0.25;
        let yn = (y + lam) * 0.25;
        a = an;
        x = xn;
        y = yn;
        if !(q >= a.abs()) {
            break;
        }
    }
    a = (x + y + y) / 3.0;
    s /= a;
    (1.0 + s * s * (3.0/10.0 + s * (1.0 / 7.0 + s * (3.0/ 8.0 + 
        s * (9.0/22.0 + s * (159.0/208.0 + s * 9.0/8.0)))))) / a.sqrt()
}

pub fn elliprf(x: f64, y: f64, z: f64) -> f64 {
    //! Carlson symmetric elliptic integral of the first kind
    //!
    //! Uses technique from 
    //!
    //! B. C. Carlson, “Numerical computation of real or complex elliptic 
    //! integrals,” Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
    //! 
    //! https://arxiv.org/abs/math/9409227 https://doi.org/10.1007/BF02198293
    //!
    //! # Domain
    //!
    //! Returns `NAN` if any paramters are less than zero and infinity if
    //! two or more of the paramters are equal to zero.

    // Peak: 6.92703689768557e-16
    // RMS: 1.6860763161239592e-16
    // # Samples: 970299

    // https://dlmf.nist.gov/19.20#i
    if x <= 0.0 || y <= 0.0 || z <= 0.0 {
        if (x == 0.0 && (y == 0.0 || z == 0.0)) || (y == 0.0 && z == 0.0) {
            // at least two parameters must be nonzero
            return f64::INFINITY;
        } else if x < 0.0 || y < 0.0 || z < 0.0 {
            // x, y, and z must be greater than or equal to 0.0
            return f64::NAN;
        }
        // Remaining case is one parameter is equal to 0.0
    }

    let r: f64 = 1e-16;
    let mut an = (x + y + z) / 3.0;
    let mut q = (3.0 * r).powf(-1.0 / 6.0) * (an - x).abs().max((an - y).abs()).max((an - z).abs());
    let mut a0_x = an - x;
    let mut a0_y = an - y;

    // https://dlmf.nist.gov/19.26#E18
    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    while q >= an.abs() {
        let xs = xn.sqrt();
        let ys = yn.sqrt();
        let zs = zn.sqrt();

        let lam = xs * ys + ys * zs + zs * xs;

        a0_x *= 0.25;
        a0_y *= 0.25;
        q *= 0.25;

        an = (an + lam) * 0.25;
        xn = (xn + lam) * 0.25;
        yn = (yn + lam) * 0.25;
        zn = (zn + lam) * 0.25;
    }

    let xf = a0_x / an;
    let yf = a0_y / an;
    let zf = -xf - yf;

    let e2 = xf * yf - zf * zf; 
    let e3 = xf * yf * zf;

    (1.0 - e2 * 1.0 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0) / an.sqrt()
}

pub fn elliprj(x: f64, y: f64, z: f64, p: f64) -> f64 {
    //! Carlson symmetric elliptic integral of the third kind
    //!
    //! Uses technique from 
    //!
    //! B. C. Carlson, “Numerical computation of real or complex elliptic 
    //! integrals,” Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
    //! 
    //! https://arxiv.org/abs/math/9409227 https://doi.org/10.1007/BF02198293
    //!
    //! with special processing near asmyptotes.
    //!
    //! # Domain
    //!
    //! Returns `NAN` if `x`, `y`, or `z` are less than zero and infinity if
    //! two or more of `x`, `y`, and `z` are equal to zero or `p` is equal 
    //! to zero.

    // 1 <= x < 25; 1 <= y < 25; 1 <= z < 25, 1 <= p < 25
    // Peak: 8.43612808234365e-16
    // RMS: 1.9681268183030496e-16
    // # Samples: 324888


    // https://dlmf.nist.gov/19.16#E5
    if x <= 0.0 || y <= 0.0 || z <= 0.0 {
        if (x == 0.0 && y == 0.0) || z == 0.0 {
            // at least two parameters must be nonzero and z must be nonzero
            return f64::INFINITY;
        } else if x < 0.0 || y < 0.0 || z < 0.0 {
            // x, y, and z must be greater than or equal to 0.0
            return f64::NAN;
        } 
        // Remaining case is one parameter is equal to 0.0
    } 
    
    if p < 0.0 {
        let x1 = x.min(y).min(z);
        let z1 = x.max(y).max(z);
        let y1 = x.min(y).max(x.max(y).min(z));
        let q = -p;
        let p_y = (z1 - y1) * (y1 - x1) / (y1 + q);
        let p1 = p_y + y1;
        let rhs = p_y * elliprj(x1, y1, z1, p1) - 3.0 * elliprf(x1, y1, z1) 
                    + 3.0 * ((x1 * y1 * z1) / (x1 * z1 + p1 * q)).sqrt() *
                        elliprc(x1 * z1 + p1 * q, p1 * q);
        return rhs / (y1 + q);
    } else if p == 0.0 {
        return f64::INFINITY;
    } 

    // Check for asmyptotes:

    let c = (x + y + z) / 3.0;
    if (c / p).sqrt() < 1e-20 {
        // https://dlmf.nist.gov/19.27#E11
        return 3.0 * elliprf(x, y, z) / p;
    }
    let f = (x * y * z).powf(1.0/3.0);
    if p / f < 1e-20 {
        // https://dlmf.nist.gov/19.27#E12
        return 3.0 * ((4.0 * f / p).ln() - 2.0)/ (2.0 * (x * y * z).sqrt());
    }

    let r: f64 = 1e-16;
    let a0 = (x + y + z + 2.0 * p) / 5.0;
    let mut an = a0;
    let mut delta = (p - x) * (p - y) * (p - z);
    let mut q = (0.25 * r).powf(-1.0 / 6.0) * (an - x).abs().max((an - y).abs()).max((an - z).abs()).max((an - p).abs());

    // https://dlmf.nist.gov/19.26#E18
    let mut xn = x;
    let mut yn = y;
    let mut zn = z;
    let mut pn = p;

    let mut pow4 = 1.0;
    let mut acc = 0.0;

    while q >= an.abs() {
        let xs = xn.sqrt();
        let ys = yn.sqrt();
        let zs = zn.sqrt();
        let ps = pn.sqrt();

        let lam = xs * ys + ys * zs + zs * xs;
        let dm = (ps + xs) * (ps + ys) * (ps + zs);
        let em = delta / (dm * dm);
        acc += pow4 / dm * elliprc(1.0, 1.0 + em);

        q *= 0.25;
        pow4 *= 0.25;
        delta *= 0.015625; // 4.powi(-3)

        an = (an + lam) * 0.25;
        xn = (xn + lam) * 0.25;
        yn = (yn + lam) * 0.25;
        zn = (zn + lam) * 0.25;
        pn = (pn + lam) * 0.25;
    }

    let xf = (a0 - x) * pow4 / an;
    let yf = (a0 - y) * pow4 / an;
    let zf = (a0 - z) * pow4 / an;
    let p = 0.5 * (-xf - yf - zf);
    let p2 = p * p;
    let p3 = p2 * p;

    let e2 = xf * yf + xf * zf + yf * zf -3.0 * p2;
    let e3 = xf * yf * zf + 2.0 * e2 * p + 4.0 * p3;
    let e4 = (2.0 * xf * yf * zf + e2 * p + 3.0 * p3) * p;
    let e5 = xf * yf * zf * p2;

    let poly = 1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0
         - 3.0 * e4 / 22.0 - 9.0 * e2 * e3 / 52.0 + 3.0 * e5 / 26.0;

    pow4 / (an * an.sqrt()) * poly + 6.0 * acc
}

pub fn elliprd(x: f64, y: f64, z: f64) -> f64 {
    //! Carlson elliptic integral symmetric in only two variables
    //!
    //! # Domain
    //!
    //! Returns `NAN` if any paramters are less than zero and infinity if
    //! two or more of the paramters are equal to zero.

    //! Peak: 7.82419490337662e-16
    //! RMS: 1.7278186787171422e-16
    //! # Samples: 970299

    // https://dlmf.nist.gov/19.16#E5
    if x <= 0.0 || y <= 0.0 || z <= 0.0 {
        if (x == 0.0 && y == 0.0) || z == 0.0 {
            // at least two parameters must be nonzero and z must be nonzero
            return f64::INFINITY;
        } else if x < 0.0 || y < 0.0 || z < 0.0 {
            // x, y, and z must be greater than or equal to 0.0
            return f64::NAN;
        }
        // Remaining case is one parameter is equal to 0.0
    }

    let r: f64 = 1e-16;
    let mut an = (x + y + 3.0 * z) / 5.0;
    let mut q = (0.25 * r).powf(-1.0 / 6.0) * (an - x).abs().max((an - y).abs()).max((an - z).abs());
    let mut a0_x = an - x;
    let mut a0_y = an - y;

    // https://dlmf.nist.gov/19.26#E18
    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    let mut scale = 1.0;
    let mut acc = 0.0;

    while q >= an.abs() {
        let xs = xn.sqrt();
        let ys = yn.sqrt();
        let zs = zn.sqrt();

        let lam = xs * ys + ys * zs + zs * xs;
        acc += scale / (zs * (zn + lam));

        a0_x *= 0.25;
        a0_y *= 0.25;
        q *= 0.25;
        scale *= 0.25;

        an = (an + lam) * 0.25;
        xn = (xn + lam) * 0.25;
        yn = (yn + lam) * 0.25;
        zn = (zn + lam) * 0.25;
    }

    let xf = a0_x / an;
    let yf = a0_y / an;
    let zf = -(xf + yf) / 3.0;

    let e2 = xf * yf - 6.0 * zf * zf;
    let e3 = (3.0 * xf * yf - 8.0 * zf * zf) * zf;
    let e4 = 3.0 * (xf * yf - zf * zf) * zf * zf;
    let e5 = xf * yf * zf * zf * zf;

    let poly = 1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0
         - 3.0 * e4 / 22.0 - 9.0 * e2 * e3 / 52.0 + 3.0 * e5 / 26.0;

    scale / (an * an.sqrt()) * poly + 3.0 * acc
}

pub fn elliprg(x: f64, y: f64, z: f64) -> f64 {
    //! Carlson symmetric elliptic integral of the second kind
    //!
    //! Uses technique from 
    //!
    //! B. C. Carlson, “Numerical computation of real or complex elliptic 
    //! integrals,” Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
    //! 
    //! https://arxiv.org/abs/math/9409227 https://doi.org/10.1007/BF02198293
    //!
    //! # Domain
    //!
    //! Returns `NAN` if any paramters are less than zero and infinity if
    //! two or more of the paramters are equal to zero.

    // Peak: 3.47377852281619e-15
    // RMS: 1.563092788497549e-16
    // # Samples: 1000000

    // https://dlmf.nist.gov/19.16#E2_5
    if x <= 0.0 || y <= 0.0 || z <= 0.0 {
        if x == 0.0 {
            if y == 0.0 {
                return 0.5 * z.sqrt();
            } else if z == 0.0 {
                return 0.5 * y.sqrt();
            }
            // Remaining case is that two variables are nonzero
        } else if y == 0.0 && z == 0.0 {
            return 0.5 * x.sqrt();
        } else if x < 0.0 || y < 0.0 || z < 0.0 {
            // x, y, and z must be greater than or equal to 0.0
            return f64::NAN;
        }
    }

    let (y, z) = if z == 0.0 {(z, y)} else {(y, z)};

    0.5 * (z * elliprf(x, y, z) - (x - z) * (y - z) * elliprd(x, y, z) / 3.0 + 
            (x * y / z).sqrt())
}

pub fn ellippi_inc(phi: f64, m: f64, n: f64) -> f64 {
    //! Incomplete elliptic integral of the third kind with m = k^2
    //!
    //! # Domain:
    //! 
    //! m <= csc(phi)^2
    //! if cos(phi) == 0, then m != csc^2(phi)
    //! n != csc^2(phi)
    
    //let c = phi.sin().powi(-2);
    let s = phi.sin();
    let c = phi.cos();
    let x = c * c;
    let y = 1.0 - m * s * s;
    let rf = elliprf(x, y, 1.0);
    let rj = elliprj(x, y, 1.0, 1.0 - n * s * s);
    s * rf + n * s * s * s * rj / 3.0
}

#[cfg(test)]
mod elliprc_tests {
    use std::io::Write;
    use super::*;

    #[test]
    fn elliprc_trivials() {
        // if x < 0, returns NAN
        assert!(elliprc(-1.0, 1.0).is_nan());
        assert!(elliprc(-10.0, -10.0).is_nan());

        // If y == 0, returns infinity
        assert_eq!(elliprc(1.0, 0.0), f64::INFINITY);

        // If x == 0 and y < 0, returns 0.0
        assert_eq!(elliprc(0.0, -1.0), 0.0);
    }

    #[test]
    fn elliprc_values() {
        assert_eq!(elliprc(2.0, 3.0), 0.6154797086703873);
        assert_eq!(elliprc(3.0, 2.0), 0.6584789484624084);
        assert_eq!(elliprc(5.0, 5.0), 0.44721359549995793);
        assert_eq!(elliprc(5.0, -3.0), 0.37934451107024508);
        assert_eq!(elliprc(3.0, -5.0), 0.25198049661585725);
    }

    #[test]
    fn elliprc_relations() {
        for x in (0..100).map(|i| i as f64) {
            for y in (-100..100).map(|i| i as f64) {
                let a = elliprc(x, y);

                // https://dlmf.nist.gov/19.2#E17
                let b = if 0.0 <= x && x < y { // case_1
                    (x / y).sqrt().acos() / (y - x).sqrt()

                } else if 0.0 < y && y < x { // case_2
                    (x / y).sqrt().acosh() / (x - y).sqrt()

                } else if  y < 0.0 && 0.0 <= x { // case_3
                    if x == 0.0{
                        // https://dlmf.nist.gov/19.6#E15
                        0.0
                    } else {
                        (x / (x - y)).sqrt() * elliprc(x - y, -y)
                    }

                } else if x == y {
                    // https://dlmf.nist.gov/19.6#E15
                    1.0 / x.sqrt()
                } else if y == 0.0 {
                    f64::INFINITY
            
                } else {
                    f64::NAN
            
                };

                if b == 0.0 {
                    assert!(a < 1e-13);
                } else if b.is_infinite() {
                    assert_eq!(a, b);
                } else {
                    let err = (a - b) / b;
                    assert!(err.abs() < 1e-14);
                }
            }

        }
    }

    // #[test]
    // fn elliprc_write() {

    //     let mut f = std::fs::File::create("elliprc.bin").unwrap();
    //     for x in (0..1000).map(|i| i as f64 * 0.1) {
    //         for y in (-1000..1000).map(|i| i as f64 * 0.1) {
    //             if y != 0.0 {
    //                 let v = elliprc(x, y);
    //                 f.write_all(&(x).to_ne_bytes());
    //                 f.write_all(&y.to_ne_bytes());
    //                 f.write_all(&v.to_ne_bytes());
    //             }
    //         }
    //     }
    //     f.sync_all();
    // }

}

#[cfg(test)]
mod ellippi_inc_tests {
    use super::*;

    #[test]
    fn ellippi_inc_domain() {
        ellippi_inc(std::f64::consts::FRAC_PI_2, 1.0, 1.0);
        // -pi/2 to pi/2 in 21 steps
        for phi in (-10..=10).map(|i| i as f64 * std::f64::consts::FRAC_PI_2 * 0.1) {
            if phi == 0.0 {
                // m unbounded when csc(phi) is infinite
                assert!(ellippi_inc(phi, 1e20, 1.0).is_finite());
            } else {
                let csc2 = phi.sin().powi(-2);

                // m <= csc(phi)^2
                assert!(ellippi_inc(phi, csc2 - 1e-4, 1.1).is_finite());
                assert!(ellippi_inc(phi, csc2 + 1e-4, 1.1).is_nan());

                // if phi.cos().powi(2) == 0.0 {
                //     assert_eq!(ellippi_inc(phi, csc2, 1.1), f64::INFINITY);
                // } else {
                //     println!("{} {} {}", phi, csc2, ellippi_inc(phi, csc2, 1.1));
                //     assert!(ellippi_inc(phi, csc2, 1.1).is_finite());
                // }
            }
        }
    }

    #[test]
    fn ellippi_inc_values() {
        assert_eq!(ellippi_inc(0.2, 3.0, 0.1), 0.20447159058227227);
        assert_eq!(ellippi_inc(0.3, 3.0, 5.0), 0.3826427688731495);
    }
}

#[cfg(test)]
mod ellippi_tests {
    use super::*;

    #[test]
    fn ellippi_trivials() {
        assert!(ellippi(0.5, 1.0).is_nan());
        assert!(ellippi(0.5, 10.0).is_nan());
        assert!(ellippi(1.0, 0.5).is_nan());

        //assert_eq!(ellippi(0.0, -f64::INFINITY), 0.0);
        assert!(ellippi(0.0, f64::INFINITY).is_nan());
    }

    #[test]
    fn ellippi_values() {
        assert_eq!(ellippi(0.5, 0.5), 2.7012877620953506);
        assert_eq!(ellippi(-100.0, -50.0), 0.11250440171406012);

        assert_eq!(ellippi(-8.0, -80.0), 0.22170030862941437);
        assert_eq!(ellippi(-3.0, -15.0), 0.4376884507289381);

        assert_eq!(ellippi(2.0, 0.5), -0.31354468346518394);
    }
}

#[cfg(test)]
mod elliprf_tests {
    use std::io::Write;
    use super::*;

    #[test]
    fn elliprf_trivials() {
        assert!(elliprf(-1e-20, 1.0, 1.0).is_nan());
        assert!(elliprf(1.0, -1e-20, 1.0).is_nan());
        assert!(elliprf(1.0, 1.0, -1e-20).is_nan());
        assert_eq!(elliprf(0.0, 0.0, 1.0), f64::INFINITY);
        assert_eq!(elliprf(0.0, 1.0, 0.0), f64::INFINITY);
        assert_eq!(elliprf(1.0, 0.0, 0.0), f64::INFINITY);
    }

    #[test]
    fn elliprf_domain() {
        assert_eq!(elliprf(1e-300, 1e-300, 0.0), 1.5707963267948969e+150);
        assert_eq!(elliprf(1e-300, 0.0, 1e-300), 1.5707963267948969e+150);
        assert_eq!(elliprf(0.0, 1e-300, 1e-300), 1.5707963267948969e+150);

        assert_eq!(elliprf(1e300, 1e-300, 0.0), 6.921618222593337e-148);
        assert_eq!(elliprf(1e-300, 0.0, 1e300), 6.921618222593337e-148);
        assert_eq!(elliprf(0.0, 1e300, 1e-300), 6.921618222593336e-148);
    }

    // #[test]
    // fn elliprf_write() {
    //     let mut f = std::fs::File::create("elliprf_test3.bin").unwrap();
    //     for x in (1..100).map(|i| i as f64) {
    //         for y in (1..100).map(|i| i as f64) {
    //             for z in (1..100).map(|i| i as f64) {
    //                 let v = elliprf(x, y, z);
    //                 f.write_all(&x.to_ne_bytes());
    //                 f.write_all(&y.to_ne_bytes());
    //                 f.write_all(&z.to_ne_bytes());
    //                 f.write_all(&v.to_ne_bytes());
    //             }
    //         }
    //     }
    //     f.sync_all();
    // }

    #[test]
    fn elliprf_values() {
        assert_eq!(elliprf(1.0, 3.0, 3.0), 0.67551085885604);
        assert_eq!(elliprf(10.0, 5.0, 20.0), 0.30637969128266374);
        assert_eq!(elliprf(31002.856213369905, 31003.856213369905, 
                31004.856213369905), 0.005679265120850845);


        // First lemniscate constant
        assert_eq!(elliprf(0.0, 1.0, 2.0), 1.31102877714606);
    }

    // #[test]
    // fn elliprf_timing() {
    //     let now = std::time::Instant::now();
    //     let mut s = 0.0;
    //     for x in (1..100).map(|i| i as f64) {
    //         for y in (1..100).map(|i| i as f64) {
    //             for z in (1..100).map(|i| i as f64) {
    //                 s += elliprf(x, y, z);
    //             }
    //         }
    //     }
    //     println!("t={}, s={}", now.elapsed().as_micros(), s);
    // }
}

#[cfg(test)]
mod elliprj_tests {
    use std::io::Write;
    use super::*;

    #[test]
    fn elliprj_values() {
        assert_eq!(elliprj(0.0, 1.0, 2.0, 3.0),  0.7768862377858231);
        assert_eq!(elliprj(2.0, 3.0, 4.0, 5.0),  0.1429757966715675);
        assert_eq!(elliprj(2.0, 3.0, 4.0, -0.5),  0.2472381970305159);
        assert_eq!(elliprj(2.0, 3.0, 4.0, -5.0),  -0.12711230042963909);
    }

    #[test]
    fn elliprj_domain() {
        //assert_eq!(elliprj(0.0, 1e-300, 1.0, 1.0), 1037.3221749306806);
        assert_eq!(elliprj(0.0, 1e-130, 1e-130, 1.0), 4.71238898038469e65);
        assert_eq!(elliprj(0.0, 1e-300, 1e-300, 1.0), 4.712388980384691e150);
        assert_eq!(elliprj(1.0, 1.0, 1.0, 1e-43), 147.5961800397958);

        assert_eq!(elliprj(0.02447174185242323, 0.0, 1.0, -0.07308108396233459),
                -49.037105497507994);
        

    }

    // #[test]
    // fn elliprj_write() {
    //     let mut f = std::fs::File::create("elliprj.bin").unwrap();
    //     for x in (1..25).map(|i| i as f64) {
    //         for y in (1..25).map(|i| i as f64) {
    //             for z in (1..25).map(|i| i as f64) {
    //                 for p in (-25..25).map(|i| i as f64) {
    //                     let v = elliprj(x, y, z, p);
    //                     f.write_all(&x.to_ne_bytes());
    //                     f.write_all(&y.to_ne_bytes());
    //                     f.write_all(&z.to_ne_bytes());
    //                     f.write_all(&p.to_ne_bytes());
    //                     f.write_all(&v.to_ne_bytes());
    //                 }
    //             }
    //         }
    //     }
    //     f.sync_all();
    // }
}

#[cfg(test)]
mod elliprd_tests {
    use std::io::Write;
    use super::*;

    #[test]
    fn elliprd_values() {
        assert_eq!(elliprd(1.0, 2.0, 3.0), 0.29046028102899063);

        // Second lemniscate constant
        assert_eq!(elliprd(0.0, 2.0, 1.0) / 3.0, 0.5990701173677963);
    }

    #[test]
    fn elliprd_trivials() {
        assert!(elliprd(-1e-20, 1.0, 1.0).is_nan());
        assert!(elliprd(1.0, -1e-20, 1.0).is_nan());
        assert!(elliprd(1.0, 1.0, -1e-20).is_nan());
        assert_eq!(elliprd(0.0, 0.0, 1.0), f64::INFINITY);
        assert_eq!(elliprd(0.0, 1.0, 0.0), f64::INFINITY);
        assert_eq!(elliprd(1.0, 0.0, 0.0), f64::INFINITY);
        assert_eq!(elliprd(1.0, 1.0, 0.0), f64::INFINITY);
    }

    #[test]
    fn elliprd_degenerate() {

        // https://dlmf.nist.gov/19.20#E18
        for x in (1..100).map(|i| i as f64 * 0.5) {
            let a1 = elliprd(0.0, x  , x  );
            let a2 = elliprd(x  , 0.0, x  );
            let b = 3.0 * std::f64::consts::FRAC_PI_4 / (x * x.sqrt());
            assert!((b - a1) / b < 1e-15);
            assert!((b - a2) / b < 1e-15);
        }
    }

    // #[test]
    // fn elliprd_write() {
    //     let mut f = std::fs::File::create("elliprd.bin").unwrap();
    //     for x in (1..100).map(|i| i as f64) {
    //         for y in (1..100).map(|i| i as f64) {
    //             for z in (1..100).map(|i| i as f64) {
    //                 let v = elliprd(x, y, z);
    //                 f.write_all(&x.to_ne_bytes());
    //                 f.write_all(&y.to_ne_bytes());
    //                 f.write_all(&z.to_ne_bytes());
    //                 f.write_all(&v.to_ne_bytes());
    //             }
    //         }
    //     }
    //     f.sync_all();
    // }

    // #[test]
    // fn elliprd_timing() {
    //     let now = std::time::Instant::now();
    //     let mut s = 0.0;
    //     for x in (1..100).map(|i| i as f64) {
    //         for y in (1..100).map(|i| i as f64) {
    //             for z in (1..100).map(|i| i as f64) {
    //                 s += elliprd(x, y, z);
    //             }
    //         }
    //     }
    //     println!("t={}, s={}", now.elapsed().as_micros(), s);

    //     let now = std::time::Instant::now();
    //     let mut s = 0.0;
    //     for x in (1..100).map(|i| i as f64) {
    //         for y in (1..100).map(|i| i as f64) {
    //             for z in (1..100).map(|i| i as f64) {
    //                 s += elliprd2(x, y, z);
    //             }
    //         }
    //     }
    //     println!("t={}, s={}", now.elapsed().as_micros(), s);
    //     todo!();
    // }
}

#[cfg(test)]
mod elliprg_tests {
    use std::io::Write;
    use super::*;

    #[test]
    fn elliprg_values() {
        assert_eq!(elliprg(1.0, 2.0, 3.0), 1.4018470999908951);

        // Second lemniscate constant
        //assert_eq!(elliprd(0.0, 2.0, 1.0) / 3.0, 0.599070117367796);
    }

    #[test]
    fn elliprd_trivials() {
        assert!(elliprg(-1e-20, 1.0, 1.0).is_nan());
        assert!(elliprg(1.0, -1e-20, 1.0).is_nan());
        assert!(elliprg(1.0, 1.0, -1e-20).is_nan());
        assert_eq!(elliprg(0.0, 0.0, 0.0), 0.0);
        assert_eq!(elliprg(1.0, 1.0, 1.0), 1.0);
    }

    #[test]
    fn elliprg_degenerate() {

        // https://dlmf.nist.gov/19.20#E18
        for x in (1..100).map(|i| i as f64 * 0.5) {
            let a1 = elliprg(0.0, x  , x  );
            let a2 = elliprg(x  , 0.0, x  );
            let a3 = elliprg(x  , x  , 0.0);
            let b = std::f64::consts::FRAC_PI_4 * x.sqrt();
            assert!((b - a1) / b < 1e-15);
            assert!((b - a2) / b < 1e-15);
            assert!((b - a3) / b < 1e-15);
        }
    }

    // #[test]
    // fn elliprg_write() {
    //     let mut f = std::fs::File::create("elliprg.bin").unwrap();
    //     for x in (0..100).map(|i| i as f64) {
    //         for y in (0..100).map(|i| i as f64) {
    //             for z in (0..100).map(|i| i as f64) {
    //                 let v = elliprg(x, y, z);
    //                 f.write_all(&x.to_ne_bytes());
    //                 f.write_all(&y.to_ne_bytes());
    //                 f.write_all(&z.to_ne_bytes());
    //                 f.write_all(&v.to_ne_bytes());
    //             }
    //         }
    //     }
    //     f.sync_all();
    // }
}