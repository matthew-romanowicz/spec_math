/*                                                     incbi()
*
*      Inverse of incomplete beta integral
*
*
*
* SYNOPSIS:
*
* double a, b, x, y, incbi();
*
* x = incbi( a, b, y );
*
*
*
* DESCRIPTION:
*
* Given y, the function finds x such that
*
*  incbet( a, b, x ) = y .
*
* The routine performs interval halving or Newton iterations to find the
* root of incbet(a,b,x) - y = 0.
*
*
* ACCURACY:
*
*                      Relative error:
*                x     a,b
* arithmetic   domain  domain  # trials    peak       rms
*    IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
*    IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
*    IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
*    VAX       0,1    .5,100     25000    3.5e-14   1.1e-15
* With a and b constrained to half-integer or integer values:
*    IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
*    IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
* With a = .5, b constrained to half-integer or integer values:
*    IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11
*/


/*
* Cephes Math Library Release 2.4:  March,1996
* Copyright 1984, 1996 by Stephen L. Moshier
*/

use crate::cephes64::consts::{MACHEP, MAXLOG, MINLOG};
use crate::cephes64::ndtri::ndtri;
use crate::cephes64::gamma::lgam;
use crate::cephes64::incbet::incbet;

pub fn incbi(aa: f64, bb: f64, yy0: f64) -> f64
{
    // double a, b, y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh, xt;
    // int i, rflg, dir, nflg;


    //let i: isize = 0;
    if yy0 <= 0.0 {
        return 0.0;
    } else if yy0 >= 1.0 {
        return 1.0;
    }
    let mut x0: f64 = 0.0;
    let mut yl: f64 = 0.0;
    let mut x1: f64 = 1.0;
    let mut yh: f64 = 1.0;
    let mut nflg: isize = 0;

    let mut a: f64;
    let mut b: f64;
    let mut y0: f64;
    let mut d: f64;
    let mut y: f64;
    let mut x: f64;
    let mut lgm: f64;
    let mut yp: f64;
    let mut di: f64;
    let mut dithresh: f64;
    let mut xt: f64;

    let mut rflg: isize;
    let mut dir: isize;

    let mut goto_newt = false;

    if aa <= 1.0 || bb <= 1.0 {
        dithresh = 1.0e-6;
        rflg = 0;
        a = aa;
        b = bb;
        y0 = yy0;
        x = a / (a + b);
        y = incbet(a, b, x);
        //goto ihalve;
    } else {
        dithresh = 1.0e-4;
        /* approximation to inverse function */

        yp = -ndtri(yy0);

        if yy0 > 0.5 {
            rflg = 1;
            a = bb;
            b = aa;
            y0 = 1.0 - yy0;
            yp = -yp;
        } else {
            rflg = 0;
            a = aa;
            b = bb;
            y0 = yy0;
        }

        lgm = (yp * yp - 3.0) / 6.0;
        x = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
        d = yp * (x + lgm).sqrt() / x
            - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0))
            * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * x));
        d = 2.0 * d;
        if d < MINLOG {
            // TODO: Find a test for this case

            x = 1.0;
            //goto under;
            //sf_error("incbi", SF_ERROR_UNDERFLOW, NULL);
            x = 0.0;
            if rflg != 0 {
                if x <= MACHEP {
                    return 1.0 - MACHEP;
                } else {
                    return 1.0 - x;
                }
            } else {
                return x;
            }
        }
        x = a / (a + b * d.exp());
        y = incbet(a, b, x);
        yp = (y - y0) / y0;
        if yp.abs() < 0.2 {
            //goto newt;
            goto_newt = true;
        } 
    }

    'ihalve: loop { //ihalve
    /* Resort to interval halving if not close enough. */

        if !goto_newt {
            dir = 0;
            di = 0.5;
            for i in 0..100 {
                if i != 0 {
                    x = x0 + di * (x1 - x0);
                    if x == 1.0 {
                        // TODO: Find a test for this case
                        x = 1.0 - MACHEP;
                    }
                    if x == 0.0 {
                        di = 0.5;
                        x = x0 + di * (x1 - x0);
                        if x == 0.0 {
                            break 'ihalve; //goto under;
                        } 
                    }
                    y = incbet(a, b, x);
                    yp = (x1 - x0) / (x1 + x0);
                    if yp.abs() < dithresh {
                        goto_newt = true;
                        break; //goto newt; 
                    }
                    yp = (y - y0) / y0;
                    if yp.abs() < dithresh {
                        goto_newt = true;
                        break; //goto newt; 
                    }
                }
                if y < y0 {
                    x0 = x;
                    yl = y;
                    if dir < 0 {
                        dir = 0;
                        di = 0.5;
                    } else if dir > 3 {
                        // TODO: Find a test for this case
                        di = 1.0 - (1.0 - di) * (1.0 - di);
                    } else if dir > 1 {
                        di = 0.5 * di + 0.5;
                    } else {
                        di = (y0 - y) / (yh - yl);
                    }
                    dir += 1;
                    if x0 > 0.75 {
                        if rflg == 1 {
                            rflg = 0;
                            a = aa;
                            b = bb;
                            y0 = yy0;
                        } else {
                            rflg = 1;
                            a = bb;
                            b = aa;
                            y0 = 1.0 - yy0;
                        }
                        x = 1.0 - x;
                        y = incbet(a, b, x);
                        x0 = 0.0;
                        yl = 0.0;
                        x1 = 1.0;
                        yh = 1.0;
                        continue 'ihalve; //goto ihalve;
                    } 
                }
                else {
                    x1 = x;
                    if rflg == 1 && x1 < MACHEP {
                        x = 0.0;
                        break 'ihalve; //goto done;
                    }
                    yh = y;
                    if dir > 0 {
                        dir = 0;
                        di = 0.5;
                    } else if dir < -3 {
                        di = di * di;
                    } else if dir < -1 {
                        di = 0.5 * di;
                    } else {
                        di = (y - y0) / (yh - yl);
                    }
                    dir -= 1;
                }
            }
        } 

        if !goto_newt {
            // TODO: Find a test for this case

            //sf_error("incbi", SF_ERROR_LOSS, NULL);
            if x0 >= 1.0 {
                x = 1.0 - MACHEP;
                break; //goto done;
            } else if x <= 0.0 {
                //under:
                //sf_error("incbi", SF_ERROR_UNDERFLOW, NULL);
                x = 0.0;
                break; //goto done;
            } 
        
        } else {
            goto_newt = false;
        }

        //newt:

        if nflg != 0 {
            break; //goto done;
        }
        nflg = 1;
        lgm = lgam(a + b) - lgam(a) - lgam(b);

        for i in 0..8 {
            /* Compute the function at this point. */
            if i != 0 {
                y = incbet(a, b, x);
            }
            if y < yl {
                x = x0;
                y = yl;
            } else if y > yh {
                x = x1;
                y = yh;
            } else if y < y0 {
                x0 = x;
                yl = y;
            } else {
                x1 = x;
                yh = y;
            }
            if x == 1.0 || x == 0.0 {
                // TODO: Find a test for this case
                break;
            }
            /* Compute the derivative of the function at this point. */
            d = (a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() + lgm;
            if d < MINLOG {
                break 'ihalve; //goto done;
            } else if d > MAXLOG {
                // TODO: Find a test for this case
                break;
            }
            d = d.exp();
            /* Compute the step to the next approximation of x. */
            d = (y - y0) / d;
            xt = x - d;
            if xt <= x0 {
                y = (x - x0) / (x1 - x0);
                xt = x0 + 0.5 * y * (x - x0);
                if xt <= 0.0 {
                    // TODO: Find a test for this case
                    break;
                }
            } 
            if xt >= x1 {
                y = (x1 - x) / (x1 - x0);
                xt = x1 - 0.5 * y * (x1 - x);
                if xt >= 1.0 {
                    // TODO: Find a test for this case
                    break;
                }
            }
            x = xt;
            if (d / x).abs() < 128.0 * MACHEP {
                break 'ihalve; //goto done;
            } 
        }
        /* Did not converge.  */
        dithresh = 256.0 * MACHEP;
        //goto ihalve;
    }

    //done:

    if rflg != 0 {
        if x <= MACHEP {
            1.0 - MACHEP
        } else {
            1.0 - x
        }
    } else {
        x
    }
}

#[cfg(test)]
mod incbi_tests {
    use super::*;

    #[test]
    fn incbi_trivials() {
        assert_eq!(incbi(0.0, 0.0, -1e-10), 0.0);
        assert_eq!(incbi(0.0, 0.0, 1.0 + 1e-10), 1.0);
    }

    #[test]
    fn incbi_small_a_b() {
        assert_eq!(incbi(0.1, 0.1, 0.549), 0.7541648230048816);
        assert_eq!(incbi(0.9, 0.1, 0.1), 0.6052324636813339);
        assert_eq!(incbi(0.1, 0.9, 0.1), 1.1794379185012221e-10);
        assert_eq!(incbi(0.1, 0.1, 0.988), 0.9999999999999999);
        assert_eq!(incbi(0.1, 1.7, 0.995), 0.7919648901556371);
        assert_eq!(incbi(0.1, 23.3, 0.999), 0.136843903826126);
        assert_eq!(incbi(0.1, 90.8, 0.976), 0.011043097040208618);
        assert_eq!(incbi(0.1, 97.8, 0.977), 0.010524395107847806);
        assert_eq!(incbi(98.2, 0.1, 0.024), 0.989788604268831);
        assert_eq!(incbi(0.3, 0.3, 0.5), 0.5);
        assert_eq!(incbi(1e-10, 1e-10, 0.995), 0.9999999999999999);
        assert_eq!(incbi(1e-10, 1e-10, 1e-5), 0.0);
    }

    #[test]
    fn incbi_nominal() {
        assert_eq!(incbi(2.9, 1.1, 0.1), 0.43152176790060076);
        assert_eq!(incbi(1.1, 2.9, 0.1), 0.04544357924658877);
        assert_eq!(incbi(10.0, 35.0, 0.5), 0.2180812669098103);
        assert_eq!(incbi(10.0, 1e9, 0.9), 1.4205990161997306e-08);
        assert_eq!(incbi(1e9, 10.0, 0.5), 0.9999999903312855);
        assert_eq!(incbi(1.5, 99.9, 0.995), 0.06208488307711709);
        assert_eq!(incbi(1e9, 1e9, 0.995), 0.5000287986753604);
        assert_eq!(incbi(1e10, 1e10, 1e-5), 0.4999849213567364);
        assert_eq!(incbi(1e50, 1e50, 1e-5), 0.5000025);
        assert_eq!(incbi(1e50, 1e50, 1e-15), 0.5000000000000002);
        assert_eq!(incbi(1e50, 1e-50, 1e-15), 0.9999999999999999);
        assert_eq!(incbi(1e-50, 1e50, 1e-15), 0.0);
        assert_eq!(incbi(0.5, 1e50, 0.5), 2.2746821155978688e-51);
        assert_eq!(incbi(1e50, 1e50, 1.0 - 1e-10), 0.499999999975);
        // for i in 1..1000 {
        //     for j in 1..1000 {
        //         for k in 1..1000 {
        //             println!("i: {}, j: {}, k: {}", i, j, k);
        //             incbi((1000 - j) as f64 / 10.0, (1000 - i) as f64 / 10.0, k as f64 / 1000.0);
        //             println!("case 2");
        //             incbi((1000 - i) as f64 / 10.0, (1000 - j) as f64 / 10.0, k as f64 / 1000.0);
        //         }
        //     }
        // }
    }
}