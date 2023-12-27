/*                                                     beta.c
*
*     Beta function
*
*
*
* SYNOPSIS:
*
* double a, b, y, beta();
*
* y = beta( a, b );
*
*
*
* DESCRIPTION:
*
*                   -     -
*                  | (a) | (b)
* beta( a, b )  =  -----------.
*                     -
*                    | (a+b)
*
* For large arguments the logarithm of the function is
* evaluated using lgam(), then exponentiated.
*
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE       0,30       30000       8.1e-14     1.1e-14
*
* ERROR MESSAGES:
*
*   message         condition          value returned
* beta overflow    log(beta) > MAXLOG       0.0
*                  a or b <0 integer        0.0
*
*/


/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1984, 1987 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

// #include "mconf.h"

const MAXGAM: f64 = 171.624376956302725;

use crate::cephes64::consts::MAXLOG;
use crate::cephes64::gamma::{gamma, lgam_sgn};

const ASYMP_FACTOR: f64 = 1e6;

// static double lbeta_asymp(double a, double b, int *sgn);
// static double lbeta_negint(int a, double b);
// static double beta_negint(int a, double b);

pub fn beta(a: f64, b: f64) -> f64 {
    //! Beta function
    //!
    //! ## DESCRIPTION:
    //!
    #![doc=include_str!("beta.svg")]
    //!
    //! For large arguments the logarithm of the function is
    //! evaluated using [`lgam`](crate::cephes64::lgam), then exponentiated.
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
    //!     <td>30000</td>
    //!     <td>8.1e-14</td>
    //!     <td>1.1e-14</td>
    //! </tr>
    //!</table>

    let mut sign: isize = 1;

    if a <= 0.0 {
        if a == a.floor() {
            if a == (a as isize) as f64 {
                return beta_negint(a as isize, b);
            } else {
                //sf_error("beta", SF_ERROR_OVERFLOW, NULL);
                return sign as f64 * f64::INFINITY;
            }
        }
    }

    if b <= 0.0 {
        if b == b.floor() {
            if b == (b as isize) as f64 {
                return beta_negint(b as isize, a);
            } else {
                //sf_error("beta", SF_ERROR_OVERFLOW, NULL);
                return sign as f64 * f64::INFINITY;
            }
        }
    }

    let mut a = a;
    let mut b = b;


    if a.abs() < b.abs() {
        let y = a; 
        a = b; 
        b = y;
    }

    if a.abs() > ASYMP_FACTOR * b.abs() && a > ASYMP_FACTOR {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        let y = lbeta_asymp(a, b, &mut sign);
        return sign as f64 * y.exp();
    }

    let mut y = a + b;
    if y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM {
        let mut sgngam: isize = 0; // Value doesn't matter
        y = lgam_sgn(y, &mut sgngam);
        sign *= sgngam;		/* keep track of the sign */
        y = lgam_sgn(b, &mut sgngam) - y;
        sign *= sgngam;
        y = lgam_sgn(a, &mut sgngam) + y;
        sign *= sgngam;
        if y > MAXLOG {
            //sf_error("beta", SF_ERROR_OVERFLOW, NULL);
            return sign as f64 * f64::INFINITY;
        }
        return sign as f64 * y.exp();
    }

    y = gamma(y);
    a = gamma(a);
    b = gamma(b);
    if y == 0.0 {
        //sf_error("beta", SF_ERROR_OVERFLOW, NULL);
        return sign as f64 * f64::INFINITY;
    }

    if (a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs() {
        y = b / y;
        y * a
    } else {
        y = a / y;
        y * b
    }

    //overflow:
    //sf_error("beta", SF_ERROR_OVERFLOW, NULL);
    //return (sign as f64 * f64::INFINITY);
}


/* Natural log of |beta|. */

pub fn lbeta(a: f64, b: f64) -> f64
{
    // double y;
    // int sign;

    let mut sign: isize = 1;

    if a <= 0.0 {
        if a == a.floor() {
            if a == (a as isize) as f64 {
                return lbeta_negint(a as isize, b);
            } else {
                // sf_error("lbeta", SF_ERROR_OVERFLOW, NULL);
                return sign as f64 * f64::INFINITY;
            }
        }
    }

    if b <= 0.0 {
        if b == b.floor() {
            if b == (b as isize) as f64 {
                return lbeta_negint(b as isize, a);
            } else {
                // sf_error("lbeta", SF_ERROR_OVERFLOW, NULL);
                return sign as f64 * f64::INFINITY;
            }
        }
    }

    let mut a = a;
    let mut b = b;

    if a.abs() < b.abs() {
        let y = a; 
        a = b; 
        b = y;
    }

    if a.abs() > ASYMP_FACTOR * b.abs() && a > ASYMP_FACTOR {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        return  lbeta_asymp(a, b, &mut sign)
    }

    let mut y = a + b;
    if y.abs() > MAXGAM || a.abs() > MAXGAM || b.abs() > MAXGAM {
        let mut sgngam: isize = 0; // Value doesn't matter
        y = lgam_sgn(y, &mut sgngam);
        sign *= sgngam;		/* keep track of the sign */
        y = lgam_sgn(b, &mut sgngam) - y;
        sign *= sgngam;
        y = lgam_sgn(a, &mut sgngam) + y;
        sign *= sgngam;
        return y;
    }

    y = gamma(y);
    a = gamma(a);
    b = gamma(b);
    if y == 0.0 {
        // over:
        // sf_error("lbeta", SF_ERROR_OVERFLOW, NULL);
        return sign as f64 * f64::INFINITY;
    }

    if (a.abs() - y.abs()).abs() > (b.abs() - y.abs()).abs() {
        y = b / y;
        y *= a;
    }
    else {
        y = a / y;
        y *= b;
    }

    y.abs().ln()
}

/*
* Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
*/
fn lbeta_asymp(a: f64, b: f64, sgn: &mut isize) -> f64
{
    let mut r = lgam_sgn(b, sgn);
    r -= b * a.ln();

    r += b * (1.0 - b) / (2.0 * a);
    r += b * (1.0 - b) * (1.0 - 2.0 * b) / (12.0 * a * a);
    r + -b * b * (1.0 - b) * (1.0 - b) / (12.0 * a * a * a)
}


/*
* Special case for a negative integer argument
*/

fn beta_negint(a: isize, b: f64) -> f64
{
    if b == (b as isize) as f64 && (1 - a) as f64 - b > 0.0 {
        if b as isize % 2 == 0 {
            beta((1 - a) as f64 - b, b)
        } else {
            -beta((1 - a) as f64 - b, b)
        }
    } else {
        //sf_error("lbeta", SF_ERROR_OVERFLOW, NULL);
        f64::INFINITY
    }
}

fn lbeta_negint(a: isize, b: f64) -> f64
{
    if b == (b as isize) as f64 && (1 - a) as f64 - b > 0.0 {
        lbeta((1 - a) as f64 - b, b)
    } else {
        //sf_error("lbeta", SF_ERROR_OVERFLOW, NULL);
        return f64::INFINITY
    }
}

#[cfg(test)]
mod beta_tests {
    use super::*;

    #[test]
    fn beta_trivials() {
        assert_eq!(beta(f64::NAN, 1.5).is_nan(), true);
        assert_eq!(beta(1.5, f64::NAN).is_nan(), true);
        assert_eq!(beta(0.0, 1.0), f64::INFINITY);
        assert_eq!(beta(1.0, 0.0), f64::INFINITY);
        assert_eq!(beta(f64::INFINITY, 3.5), 0.0);
        assert_eq!(beta(-f64::INFINITY, 3.5), f64::INFINITY);
        assert_eq!(beta(3.5, f64::INFINITY), 0.0);
        assert_eq!(beta(3.5, -f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn beta_norm() {
        assert_eq!(beta(-10.1, 13.4), -1092.1135745260203);
        assert_eq!(beta(35.8, 2.0), 0.0007590478503764876);
        assert_eq!(beta(56.0, -100.01), 3.587060392013282e-31);
        assert_eq!(beta(0.01, -1e-5), -99900.0163139966);
        assert_eq!(beta(1e-20, 1e-15), 1.00001e+20);
    }

    #[test]
    fn beta_large() {
        assert_eq!(beta(100.0, 70.0), 3.740912730581706e-51);
        assert_eq!(beta(100.5, 100.2), 1.357578491028074e-61);
        assert_eq!(beta(1000.0, 100.0), 7.730325902035211e-147);
        assert_eq!(beta(-1357.9, 111.1), -4.308898979808891e-168);
        assert_eq!(beta(1357.9, -111.1), 3.708962683087477e+166)
    }

    #[test]
    fn beta_neg_int() {
        // 1 - a > b
        assert_eq!(beta(-1.0, 1.0), -1.0);
        assert_eq!(beta(-10.0, 5.0), -0.0007936507936507937);
        assert_eq!(beta(-10.0, 10.0), 0.1);
        assert_eq!(beta(-1e8, 15.0), -8.717838273725711e-110);
        assert_eq!(beta(-1e9, 25.0), -6.204485878677699e-202);

        assert_eq!(beta(1.0, -1.0), -1.0);
        assert_eq!(beta(5.0, -10.0), -0.0007936507936507937);
        assert_eq!(beta(10.0, -10.0), 0.1);
        assert_eq!(beta(15.0, -1e8), -8.717838273725711e-110);
        assert_eq!(beta(25.0, -1e9), -6.204485878677699e-202);
    }
}

#[cfg(test)]
mod lbeta_tests {
    use super::*;

    #[test]
    fn lbeta_trivials() {
        assert_eq!(lbeta(f64::NAN, 1.5).is_nan(), true);
        assert_eq!(lbeta(1.5, f64::NAN).is_nan(), true);
        assert_eq!(lbeta(0.0, 1.0), f64::INFINITY);
        assert_eq!(lbeta(1.0, 0.0), f64::INFINITY);
        assert_eq!(lbeta(f64::INFINITY, 3.5), -f64::INFINITY);
        assert_eq!(lbeta(-f64::INFINITY, 3.5), f64::INFINITY);
        assert_eq!(lbeta(3.5, f64::INFINITY), -f64::INFINITY);
        assert_eq!(lbeta(3.5, -f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn lbeta_norm() {
        assert_eq!(lbeta(-10.1, 13.4), 6.9958701568728126);
        assert_eq!(lbeta(35.8, 2.0), -7.18344573858154);
        assert_eq!(lbeta(56.0, -100.01), -70.1028048480498);
        assert_eq!(lbeta(0.01, -1e-5), 11.5119251279399);
        assert_eq!(lbeta(1e-20, 1e-15), 46.05171185983092);
    }

    #[test]
    fn lbeta_large() {
        assert_eq!(lbeta(100.0, 70.0), -116.11251011543409);
        assert_eq!(lbeta(100.5, 100.2), -140.15198808120385);
        assert_eq!(lbeta(1000.0, 100.0), -336.43485764773686);
        assert_eq!(lbeta(-1357.9, 111.1), -385.37361320863965);
        assert_eq!(lbeta(1357.9, -111.1), 383.5398776742895)
    }

    #[test]
    fn lbeta_neg_int() {
        // 1 - a > b
        assert_eq!(lbeta(-1.0, 1.0), 0.0); // TODO: Why is this zero?
        assert_eq!(lbeta(-10.0, 5.0), -7.138866999945524);
        assert_eq!(lbeta(-10.0, 10.0), -2.3025850929940455);
        assert_eq!(lbeta(-1e8, 15.0), -251.11898892654676);
        assert_eq!(lbeta(-1e9, 25.0), -463.296916225548);

        assert_eq!(lbeta(1.0, -1.0), 0.0);
        assert_eq!(lbeta(5.0, -10.0), -7.138866999945524);
        assert_eq!(lbeta(10.0, -10.0), -2.3025850929940455);
        assert_eq!(lbeta(15.0, -1e8), -251.11898892654676);
        assert_eq!(lbeta(25.0, -1e9), -463.296916225548);
    }
}