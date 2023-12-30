/*
* Implement sin(pi * x) and cos(pi * x) for real x. Since the periods
* of these functions are integral (and thus representable in double
* precision), it's possible to compute them with greater accuracy
* than sin(x) and cos(x).
*/

use crate::cephes64::M_PI;


/* Compute sin(pi * x). */
pub fn sinpi(x: f64) -> f64
{
    let mut s: f64 = 1.0;
    let mut x = x;

    if x < 0.0 {
        x = -x;
        s = -1.0;
    }

    let r = x % 2.0;
    if r < 0.5 {
        s * (M_PI * r).sin()
    } else if r > 1.5 {
        s * (M_PI * (r - 2.0)).sin()
    }
    else {
        -s * (M_PI * (r - 1.0)).sin()
    }
}


/* Compute cos(pi * x) */
pub fn cospi(x: f64) -> f64
{
    let r = x.abs() % 2.0;
    if r == 0.5 {
        // We don't want to return -0.0
        0.0
    } else if r < 1.0 {
        -(M_PI * (r - 0.5)).sin()
    }
    else {
        (M_PI * (r - 1.5)).sin()
    }
}