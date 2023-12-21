/*                                                     ellpj.c
 *
 *     Jacobian Elliptic Functions
 *
 *
 *
 * SYNOPSIS:
 *
 * double u, m, sn, cn, dn, phi;
 * int ellpj();
 *
 * ellpj( u, m, _&sn, _&cn, _&dn, _&phi );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 * Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
 * and dn(u|m) of parameter m between 0 and 1, and real
 * argument u.
 *
 * These functions are periodic, with quarter-period on the
 * real axis equal to the complete elliptic integral
 * ellpk(m).
 *
 * Relation to incomplete elliptic integral:
 * If u = ellik(phi,m), then sn(u|m) = sin(phi),
 * and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
 *
 * Computation is by means of the arithmetic-geometric mean
 * algorithm, except when m is within 1e-9 of 0 or 1.  In the
 * latter case with m close to 1, the approximation applies
 * only for phi < pi/2.
 *
 * ACCURACY:
 *
 * Tested at random points with u between 0 and 10, m between
 * 0 and 1.
 *
 *            Absolute error (* = relative error):
 * arithmetic   function   # trials      peak         rms
 *    IEEE      phi         10000       9.2e-16*    1.4e-16*
 *    IEEE      sn          50000       4.1e-15     4.6e-16
 *    IEEE      cn          40000       3.6e-15     4.4e-16
 *    IEEE      dn          10000       1.3e-12     1.8e-14
 *
 *  Peak error observed in consistency check using addition
 * theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
 * the above relation to the incomplete elliptic integral.
 * Accuracy deteriorates when u is large.
 *
 */
 
 /*                                                     ellpj.c         */
 
 
 /*
  * Cephes Math Library Release 2.0:  April, 1987
  * Copyright 1984, 1987 by Stephen L. Moshier
  * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
  */
 
 /* Scipy changes:
  * - 07-18-2016: improve evaluation of dn near quarter periods
  */
 
use super::consts::{M_PI_2, MACHEP};



pub fn ellpj(u: f64, m: f64) -> (f64, f64, f64, f64) // sn, cn, dn, ph
{
//! Jacobian elliptic functions
//!
//! ## Description:
//! 
//! Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m), and dn(u|m) of parameter `m` between `0` and `1`, and real argument `u`.
//!
//! These functions are periodic, with quarter-period on the real axis equal to the complete elliptic integral `ellpk(m)`.
//!
//! Relation to incomplete elliptic integral:
//! If u = ellik(phi,m), then sn(u|m) = sin(phi),
//! and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
//!
//! Computation is by means of the arithmetic-geometric mean algorithm, except when 
//! `m` is within `1e-9` of `0` or `1`.  In the latter case with `m` close to `1`, the approximation 
//! applies only for `phi < pi/2`.
//!
//! ## Accuracy
//!
//! Tested at random points with `u` between `0` and `10`, `m` between `0` and `1`.
//!
//! Relative Error:
//!<table>
//! <tr>
//!     <th>Arithmetic</th>
//!     <th>Function</th>
//!     <th># Trials</th>
//!     <th>Peak</th>
//!     <th>RMS</th>
//! </tr>
//! <tr>
//!     <td>IEEE</td>
//!     <td>phi</td>
//!     <td>10000</td>
//!     <td>9.2e-16*</td>
//!     <td>1.4e-16*</td>
//! </tr>
//! <tr>
//!     <td>IEEE</td>
//!     <td>sn</td>
//!     <td>10000</td>
//!     <td>4.1e-15</td>
//!     <td>4.6e-16</td>
//! </tr>
//! <tr>
//!     <td>IEEE</td>
//!     <td>cn</td>
//!     <td>10000</td>
//!     <td>3.6e-15</td>
//!     <td>4.4e-16</td>
//! </tr>
//! <tr>
//!     <td>IEEE</td>
//!     <td>dn</td>
//!     <td>10000</td>
//!     <td>1.3e-12</td>
//!     <td>1.8e-14</td>
//! </tr>
//!</table>
//!
//! Peak error observed in consistency check using addition theorem for sn(u+v) 
//! was `4e-16` (absolute).  Also tested by the above relation to the incomplete 
//! elliptic integral. Accuracy deteriorates when `u` is large.

    /* Check for special cases */
    if m < 0.0 || m > 1.0 || m.is_nan() {
        //sf_error("ellpj", SF_ERROR_DOMAIN, NULL);
        (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
    } else if m < 1.0e-9 {
        let t = u.sin();
        let b = u.cos();
        let ai = 0.25 * m * (u - t * b);
        (
            t - ai * b,
            b + ai * t,
            1.0 - 0.5 * m * t * t,
            u - ai
        )
    } else if m >= 0.9999999999 {
        // TODO: values become NaN as phi becomes large (> 100)
        let mut ai = 0.25 * (1.0 - m);
        let b = u.cosh();
        let t = u.tanh();
        let phi = 1.0 / b;
        let twon = b * u.sinh();
        let sn = t + ai * (twon - u) / (b * b);
        let ph = 2.0 * u.exp().atan() - M_PI_2 + ai * (twon - u) / b;
        ai *= t * phi;
        let cn = phi - ai * (twon - u);
        let dn = phi + ai * (twon + u);
        (sn, cn, dn, ph)

    } else {

        /* A. G. M. scale. See DLMF 22.20(ii) */
        let mut a: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut b = (1.0 - m).sqrt();
        let mut c: [f64; 9] = [m.sqrt(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut twon = 1.0;
        let mut i = 0;

        while (c[i] / a[i]).abs() > MACHEP {
            if i > 7 {
                //sf_error("ellpj", SF_ERROR_OVERFLOW, NULL);
                break;
            }
            let ai = a[i];
            i += 1;
            c[i] = (ai - b) / 2.0;
            let t = (ai * b).sqrt();
            a[i] = (ai + b) / 2.0;
            b = t;
            twon *= 2.0;
        }

        /* backward recurrence */
        let mut phi = twon * a[i] * u;
        while i >= 1 {
            let t = c[i] * phi.sin() / a[i];
            b = phi;
            phi = (t.asin() + phi) / 2.0;
            i -= 1;
        }

        let sn = phi.sin();
        let t = phi.cos();
        let cn = t;
        let dnfac = (phi - b).cos();
        /* See discussion after DLMF 22.20.5 */
        let dn = if dnfac.abs() < 0.1 {
            (1.0 - m*sn*sn).sqrt()
        }
        else {
            t / dnfac
        };
        let ph = phi;
        (sn, cn, dn, ph)
    }
}

#[cfg(test)]
mod ellpj_tests {
    use super::*;

    fn tuple_is_nan(t: (f64, f64, f64, f64)) -> bool {
        t.0.is_nan() && t.1.is_nan() && t.2.is_nan() && t.3.is_nan()
    }

    #[test]
    fn ellpj_trivials() {
        assert_eq!(tuple_is_nan(ellpj(0.0, -f64::INFINITY)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, -1e20)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, -1.0)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, -1e-20)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, 1.0 + 1e-15)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, 2.0)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, 1e20)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, f64::INFINITY)), true);
        assert_eq!(tuple_is_nan(ellpj(0.0, f64::NAN)), true);
    }

    #[test]
    fn ellpj_m_small() {
        assert_eq!(ellpj(0.0, 1e-9 - 1e-20), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 5e-10), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 1e-10), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 1e-20), (0.0, 1.0, 1.0, 0.0));

        assert_eq!(ellpj(M_PI_2, 1e-9 - 1e-20), (1.0, 3.926991429271371e-10, 0.9999999995, 1.5707963264021976));
        assert_eq!(ellpj(M_PI_2, 5e-10), (1.0, 1.9634960208170204e-10, 0.99999999975, 1.570796326598547));
        assert_eq!(ellpj(M_PI_2, 1e-10), (1.0, 3.926996940221238e-11, 0.99999999995, 1.5707963267556266));
        assert_eq!(ellpj(M_PI_2, 1e-20), (1.0, 6.123626694818465e-17, 1.0, 1.5707963267948966));

        assert_eq!(ellpj(-M_PI_2, 1e-9 - 1e-20), (-1.0, 3.926991429271371e-10, 0.9999999995, -1.5707963264021976));
        assert_eq!(ellpj(-M_PI_2, 5e-10), (-1.0, 1.9634960208170204e-10, 0.99999999975, -1.570796326598547));
        assert_eq!(ellpj(-M_PI_2, 1e-10), (-1.0, 3.926996940221238e-11, 0.99999999995, -1.5707963267556266));
        assert_eq!(ellpj(-M_PI_2, 1e-20), (-1.0, 6.123626694818465e-17, 1.0, -1.5707963267948966));
    
        assert_eq!(ellpj(1000.0, 1e-9 - 1e-20), (0.8268794000026128, 0.5623792829144593, 0.9999999996581351, 999.9999997501162));
        assert_eq!(ellpj(1000.0, 5e-10), (0.8268794702673077, 0.5623791796025811, 0.9999999998290675, 999.9999998750582));
        assert_eq!(ellpj(1000.0, 1e-10), (0.8268795264790636, 0.5623790969530785, 0.9999999999658135, 999.9999999750116));
        assert_eq!(ellpj(1000.0, 1e-20), (0.8268795405320025, 0.5623790762907029, 1.0, 1000.0));
    
        // TODO: Accuracy dwindles with very large u values
        //assert_eq!(ellpj(-1e10, 1e-9 - 1e-20), (2.6703050818551843, -0.3456454400472321, 0.9999999998811689, -9999999997.5));
        //assert_eq!(ellpj(-1e10, 5e-10), (1.578905553476641, 0.2637370913118563, 0.9999999999405844, -9999999998.75));
        //assert_eq!(ellpj(-1e10, 1e-10), (0.7057859307563442, 0.7512431164088771, 0.9999999999881168, -9999999999.75));
        //assert_eq!(ellpj(-1e10, 1e-20), (0.487506025098098, 0.8731196226709446, 1.0, -10000000000.0));
    }

    #[test]
    fn ellpj_m_large() {
        assert_eq!(ellpj(0.0, 0.9999999999), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 0.99999999995), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 0.9999999999999999), (0.0, 1.0, 1.0, 0.0));

        assert_eq!(ellpj(M_PI_2, 0.9999999999), (0.9171523356839659, 0.39853681529997453, 0.3985368154055066, 1.1608753910106864));
        assert_eq!(ellpj(M_PI_2, 0.99999999995), (0.9171523356756202, 0.3985368153191806, 0.3985368153719467, 1.1608753909897453));
        assert_eq!(ellpj(M_PI_2, 0.9999999999999999), (0.9171523356672744, 0.39853681533838664, 0.39853681533838675, 1.1608753909688043));

        assert_eq!(ellpj(-M_PI_2, 0.9999999999), (-0.9171523356839659, 0.39853681529997453, 0.3985368154055066, -1.1608753910106864));
        assert_eq!(ellpj(-M_PI_2, 0.99999999995), (-0.9171523356756202, 0.3985368153191806, 0.3985368153719467, -1.1608753909897453));
        assert_eq!(ellpj(-M_PI_2, 0.9999999999999999), (-0.9171523356672744, 0.39853681533838664, 0.39853681533838675, -1.1608753909688043));

        assert_eq!(ellpj(10.0, 0.9999999999), (0.9999999959026927, 9.052452851700366e-05, 9.107519020403078e-05, 1.5707058022662561));
        assert_eq!(ellpj(10.0, 0.99999999995), (0.9999999958901927, 9.066219392741045e-05, 9.093752477092401e-05, 1.5707056646008453));
        assert_eq!(ellpj(10.0, 0.9999999999999999), (0.9999999958776927, 9.079985903213866e-05, 9.07998596434959e-05, 1.57070552693574));

        assert_eq!(ellpj(-10.0, 0.9999999999), (-0.9999999959026927, 9.052452851700366e-05, 9.107519020403078e-05, -1.570705802266256));
        assert_eq!(ellpj(-10.0, 0.99999999995), (-0.9999999958901927, 9.066219392741045e-05, 9.093752477092401e-05, -1.570705664600845));
        assert_eq!(ellpj(-10.0, 0.9999999999999999), (-0.9999999958776927, 9.079985903213866e-05, 9.07998596434959e-05, -1.5707055269357397));
    }

    #[test]
    fn ellpj_nominal() {
        assert_eq!(ellpj(0.0, 1e-9), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 0.1), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 0.5), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 0.9), (0.0, 1.0, 1.0, 0.0));
        assert_eq!(ellpj(0.0, 0.9999999998), (0.0, 1.0, 1.0, 0.0));

        assert_eq!(ellpj(M_PI_2, 1e-9), (1.0, 3.926992677185533e-10, 0.9999999995, 1.5707963264021974));
        assert_eq!(ellpj(M_PI_2, 0.1), (0.9992196178831643, 0.03949880045549695, 0.9487655218881756, 1.5312872484118731));
        assert_eq!(ellpj(M_PI_2, 0.5), (0.9797395592908522, 0.20027579973767845, 0.7211485269903025, 1.3691569109849622));
        assert_eq!(ellpj(M_PI_2, 0.9), (0.9330319394105544, 0.3597935519708204, 0.4653023318615559, 1.2027497085690817));
        assert_eq!(ellpj(M_PI_2, 0.9999999998), (0.9171523357006574, 0.39853681526156237, 0.39853681547262654, 1.1608753910525684));

        assert_eq!(ellpj(-M_PI_2, 1e-9), (-1.0, 3.926992677185533e-10, 0.9999999995, -1.5707963264021974));
        assert_eq!(ellpj(-M_PI_2, 0.1), (-0.9992196178831643, 0.03949880045549695, 0.9487655218881756, -1.5312872484118731));
        assert_eq!(ellpj(-M_PI_2, 0.5), (-0.9797395592908522, 0.20027579973767845, 0.7211485269903025, -1.3691569109849622));
        assert_eq!(ellpj(-M_PI_2, 0.9), (-0.9330319394105544, 0.3597935519708204, 0.4653023318615559, -1.2027497085690817));
        assert_eq!(ellpj(-M_PI_2, 0.9999999998), (-0.9171523357006574, 0.39853681526156237, 0.39853681547262654, -1.1608753910525684));

        assert_eq!(ellpj(10.0, 1e-9), (-0.5440211088874444, -0.8390715303744225, 0.9999999998520201, 9.999999997614118));
        assert_eq!(ellpj(10.0, 0.1), (-0.3191099739995604, -0.9477176924031755, 0.9948954128195586, 9.74956817373828));
        assert_eq!(ellpj(10.0, 0.5), (0.8588125059527789, -0.5122900346669922, 0.7944938890951594, 8.39183082303414));
        assert_eq!(ellpj(10.0, 0.9), (-0.3030609695105877, 0.9529711689024504, 0.9577779721226483, 5.9752822590846275));
        assert_eq!(ellpj(10.0, 0.9999999998), (0.9999999959275412, 9.024919769607522e-05, 9.1350521069685e-05, 1.570706077597078));

        assert_eq!(ellpj(-10.0, 1e-9), (0.5440211088874444, -0.8390715303744225, 0.9999999998520201, -9.999999997614118));
        assert_eq!(ellpj(-10.0, 0.1), (0.3191099739995604, -0.9477176924031755, 0.9948954128195586, -9.74956817373828));
        assert_eq!(ellpj(-10.0, 0.5), (-0.8588125059527789, -0.5122900346669922, 0.7944938890951594, -8.39183082303414));
        assert_eq!(ellpj(-10.0, 0.9), (0.3030609695105877, 0.9529711689024504, 0.9577779721226483, -5.9752822590846275));
        assert_eq!(ellpj(-10.0, 0.9999999998), (-0.9999999959275412, 9.024919769607522e-05, 9.1350521069685e-05, -1.570706077597078));

        // TODO: Accuracy dwindles with very large u values
        //assert_eq!(ellpj(1e10, 1e-9), (-0.13197543311589455, -0.9912529874123316, 1.0, 9999999997.5));
        //assert_eq!(ellpj(1e10, 0.1), (0.14661588716523885, 0.9891935006007418, 0.9989246831572828, 9741726904.003225));
        //assert_eq!(ellpj(1e10, 0.5), (-0.21298576345957113, 0.9770553027150222, 0.9885942041755555, 8472130847.907698));
        //assert_eq!(ellpj(1e10, 0.9), (0.9409120825609791, -0.33865093073954655, 0.4507947460481793, 6092863472.955908));
        //assert_eq!(ellpj(1e10, 0.9999999998), (-0.9999999604071443, -0.00028139955516449293, 0.00028175469600747507, 1251366686.3300655));

        //assert_eq!(ellpj(-1e10, 1e-9), (0.13197543311589455, -0.9912529874123316, 1.0, -9999999997.5));
        //assert_eq!(ellpj(-1e10, 0.1), (-0.14661588716523885, 0.9891935006007418, 0.9989246831572828, -9741726904.003225));
        //assert_eq!(ellpj(-1e10, 0.5), (0.21298576345957113, 0.9770553027150222, 0.9885942041755555, -8472130847.907698));
        //assert_eq!(ellpj(-1e10, 0.9), (-0.9409120825609791, -0.33865093073954655, 0.4507947460481793, -6092863472.955908));
        //assert_eq!(ellpj(-1e10, 0.9999999998), (0.9999999604071443, -0.00028139955516449293, 0.00028175469600747507, -1251366686.3300655));
    }

    // #[test]
    // fn ellpk_small() {
    //     assert_eq!(ellpk(1.1e-16), 19.759320015170093); // Should actually be 19.759320015170094
    //     assert_eq!(ellpk(1.0e-20), 24.412145291060348); // Should actually round to 24.412145291060347
    // }
}