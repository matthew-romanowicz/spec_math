/*
* Cephes Math Library Release 2.2:  January, 1991
* Copyright 1984, 1991 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#![allow(clippy::excessive_precision)]

use crate::utils::frexp;

const CBRT2: f64 = 1.2599210498948731647672;
const CBRT4: f64 = 1.5874010519681994747517;
const CBRT2I: f64 = 0.79370052598409973737585;
const CBRT4I: f64 = 0.62996052494743658238361;

pub fn scalbn(x: f64, mut n: i32) -> f64 {
    let x1p1023 = f64::from_bits(0x7fe0000000000000); // 0x1p1023 === 2 ^ 1023
    let x1p53 = f64::from_bits(0x4340000000000000); // 0x1p53 === 2 ^ 53
    let x1p_1022 = f64::from_bits(0x0010000000000000); // 0x1p-1022 === 2 ^ (-1022)

    let mut y = x;

    if n > 1023 {
        y *= x1p1023;
        n -= 1023;
        if n > 1023 {
            y *= x1p1023;
            n -= 1023;
            if n > 1023 {
                n = 1023;
            }
        }
    } else if n < -1022 {
        /* make sure final n < -53 to avoid double
        rounding in the subnormal range */
        y *= x1p_1022 * x1p53;
        n += 1022 - 53;
        if n < -1022 {
            y *= x1p_1022 * x1p53;
            n += 1022 - 53;
            if n < -1022 {
                n = -1022;
            }
        }
    }
    y * f64::from_bits(((0x3ff + n) as u64) << 52)
}

pub fn ldexp(x: f64, n: i32) -> f64 {
    scalbn(x, n)
}

pub fn cbrt(x: f64) -> f64
{
    //! Cube root
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the cube root of the argument, which may be negative.
    //!
    //! Range reduction involves determining the power of 2 of
    //! the argument.  A polynomial of degree 2 applied to the
    //! mantissa, and multiplication by the cube root of 1, 2, or 4
    //! approximates the root to within about 0.1%.  Then Newton's
    //! iteration is used three times to converge to an accurate
    //! result.
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
    //!     <td>0, 1e308</td>
    //!     <td>30000</td>
    //!     <td>1.5e-16</td>
    //!     <td>5.0e-17</td>
    //! </tr>
    //!</table>

    if x.is_infinite() {
        return x;
    }
    if x == 0.0 {
        return x;
    }
    let (sign, x) = if x > 0.0 {
        (1, x)
    } else {
        (-1, -x)
    };

    let z = x;
    /* extract power of 2, leaving
    * mantissa between 0.5 and 1
    */
    let (x, mut e) = frexp(x);

    /* Approximate cube root of number between .5 and 1,
    * peak relative error = 9.2e-6
    */
    let mut x = (((-1.3466110473359520655053e-1 * x
    + 5.4664601366395524503440e-1) * x
    - 9.5438224771509446525043e-1) * x
    + 1.1399983354717293273738e0) * x + 4.0238979564544752126924e-1;

    /* exponent divided by 3 */
    if e >= 0 {
        let mut rem = e;
        e /= 3;
        rem -= 3 * e;
        if rem == 1 {
            x *= CBRT2;
        }else if rem == 2 {
            x *= CBRT4;
        }
    }


        /* argument less than 1 */

    else {
        e = -e;
        let mut rem = e;
        e /= 3;
        rem -= 3 * e;
        if rem == 1 {
            x *= CBRT2I;
        } else if rem == 2 {
            x *= CBRT4I;
        }
        e = -e;
    }

    /* multiply by power of 2 */
    x = ldexp(x, e);

    /* Newton iteration */
    x -= (x - (z / (x * x))) * 0.33333333333333333333;
    x -= (x - (z / (x * x))) * 0.33333333333333333333;

    if sign < 0 {
        -x
    } else {
        x
    }
}

#[cfg(test)]
mod cbrt_tests {
    use super::*;

    #[test]
    fn cbrt_trivials() {
        assert_eq!(cbrt(f64::INFINITY), f64::INFINITY);
        assert_eq!(cbrt(-f64::INFINITY), -f64::INFINITY);
        assert_eq!(cbrt(f64::NAN).is_nan(), true);
        assert_eq!(cbrt(1.0), 1.0);
        assert_eq!(cbrt(0.0), 0.0);
    }

    #[test]
    fn cbrt_e_neg() {
        assert_eq!(cbrt(0.5), 0.7937005259840997);
        assert_eq!(cbrt(0.421875), 0.75);
        assert_eq!(cbrt(0.008), 0.2);
        assert_eq!(cbrt(0.001), 0.1);
        assert_eq!(cbrt(1e-18), 1e-6);
        assert_eq!(cbrt(1e-300), 1e-100);

        assert_eq!(cbrt(-0.5), -0.7937005259840997);
        assert_eq!(cbrt(-0.421875), -0.75);
        assert_eq!(cbrt(-0.008), -0.2);
        assert_eq!(cbrt(-0.001), -0.1);
        assert_eq!(cbrt(-1e-18), -1e-6);
        assert_eq!(cbrt(-1e-300), -1e-100);
    }

    #[test]
    fn cbrt_e_pos() {
        assert_eq!(cbrt(8.0), 2.0);
        assert_eq!(cbrt(100.0), 4.641588833612778);
        assert_eq!(cbrt(3375.0), 15.0);
        assert_eq!(cbrt(1e18), 1e6);
        assert_eq!(cbrt(1e300), 1e100);

        assert_eq!(cbrt(-8.0), -2.0);
        assert_eq!(cbrt(-100.0), -4.641588833612778);
        assert_eq!(cbrt(-3375.0), -15.0);
        assert_eq!(cbrt(-1e18), -1e6);
        assert_eq!(cbrt(-1e300), -1e100);
    }
}