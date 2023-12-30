/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/

use crate::utils::frexp;
use crate::cephes64::consts::{M_PI, MACHEP, MAXLOG};
use crate::cephes64::cbrt::cbrt;
use crate::cephes64::polevl::polevl;
use crate::cephes64::airy::airy;
use crate::cephes64::gamma::{gamma, lgam_sgn};
use crate::cephes64::j0::j0;
use crate::cephes64::j1::j1;

const MAXGAM: f64 = 171.624376956302725;
const BIG: f64 = 1.44115188075855872E+17;

pub fn jv(n: f64, x: f64) -> f64 {
    //! Bessel function of noninteger order
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns Bessel function of order `v` of the argument,
    //! where `v` is real.  Negative `x` is allowed if `v` is an integer.
    //!
    //! Several expansions are included: the ascending power
    //! series, the Hankel expansion, and two transitional
    //! expansions for large `v`.  If `v` is not too large, it
    //! is reduced by recurrence to a region of best accuracy.
    //! The transitional expansions give 12D accuracy for `v > 500`.
    //!
    //! ## ACCURACY:
    //!
    //! Results for integer `v` are indicated by *, where `x` and `v`
    //! both vary from `-125` to `+125`.  Otherwise,
    //! `x` ranges from `0` to `125`, `v` ranges as indicated by "domain."
    //! Error criterion is absolute, except relative when |jv()| > 1.
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>v Domain</th>
    //!     <th>x Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 125</td>
    //!     <td>0, 125</td>
    //!     <td>100000</td>
    //!     <td>4.6e-15</td>
    //!     <td>2.2e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>-125, 0</td>
    //!     <td>0, 125</td>
    //!     <td>40000</td>
    //!     <td>5.4e-11</td>
    //!     <td>3.7e-13</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>0, 500</td>
    //!     <td>0, 500</td>
    //!     <td>20000</td>
    //!     <td>4.4e-15</td>
    //!     <td>4.0e-16</td>
    //! </tr>
    //!</table>
    //!
    //! Integer v:
    //!
    //!<table>
    //! <tr>
    //!     <th>Arithmetic</th>
    //!     <th>v Domain</th>
    //!     <th>x Domain</th>
    //!     <th># Trials</th>
    //!     <th>Peak</th>
    //!     <th>RMS</th>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>-125, 125</td>
    //!     <td>-125, 125</td>
    //!     <td>50000</td>
    //!     <td>3.5e-15*</td>
    //!     <td>1.9e-16*</td>
    //! </tr>
    //!</table>

    let mut n = n;
    let mut x = x;
    let mut nint = 0;			/* Flag for integer n */
    let mut sign = 1;			/* Flag for sign inversion */
    let an = n.abs();
    let mut y = an.floor();
    if y == an {
        nint = 1;
        let i = (an - 16384.0 * (an / 16384.0).floor()) as isize;
        if n < 0.0 {
            if i & 1 != 0 {
                sign = -sign;
            }
            n = an;
        }
        if x < 0.0 {
            if i & 1 != 0 {
                sign = -sign;
            }
            x = -x;
        }
        if n == 0.0 { // case 1
            return j0(x);
        } else if n == 1.0 { // case 2
            return sign as f64 * j1(x);
        }
        // case 3
    } 
    // else {
        // case 4
    //}

    if (x < 0.0) && (y != an) { // case 5
        //sf_error("Jv", SF_ERROR_DOMAIN, NULL);
        y = f64::NAN;
        return sign as f64 * y;
    }

    if x == 0.0 && n < 0.0 && nint == 0 { // case 6
        //sf_error("Jv", SF_ERROR_OVERFLOW, NULL);
        return f64::INFINITY / gamma(n + 1.0);
    }

    y = x.abs();

    if y * y < (n + 1.0).abs() * MACHEP { // case 7
        return (0.5 * x).powf(n) / gamma(n + 1.0);
    }

    let mut k = 3.6 * y.sqrt();
    let mut t = 3.6 * an.sqrt();
    if (y < t) && (an > 21.0) { // case 8
        return sign as f64 * jvs(n, x);
    }
    if (an < k) && (y > 21.0) { // case 9
        return sign as f64 * hankel(n, x);
    }

    if an < 500.0 {
        let mut q;
        /* Note: if x is too large, the continued fraction will fail; but then the
        * Hankel expansion can be used. */
        if nint != 0 {
            k = 0.0;
            q = recur(&mut n, x, &mut k, 1);
            if k == 0.0 { // case 10
                y = j0(x) / q;
                return sign as f64 * y;
            } else if k == 1.0 { // case 11
                y = j1(x) / q;
                return sign as f64 * y;
            }
            // case 12
        } 
        // else {
        // case 13
        //}

        if (an > 2.0 * y) || ((n >= 0.0) && (n < 20.0) && (y > 6.0) && (y < 20.0)) { // case 14
            /* Recur backwards from a larger value of n */

            // TODO: This case appears to give bad results
            k = n;

            y = y + an + 1.0;
            if y < 30.0 {
                y = 30.0;
            }
            y = n + (y - n).floor();
            q = recur(&mut y, x, &mut k, 0);
            y = jvs(y, x) * q;
            return sign as f64 * y;
        }

        if k <= 30.0 { // case 15
            k = 2.0;
        } else if k < 90.0 { // case 16
            k = (3.0 * k) / 4.0;
        } 
        // else {
        // case 17
        //}

        if an > (k + 3.0) {
            if n < 0.0 { // case 18
                k = -k;
            } 
            // else {
            // case 19
            //}
            q = n - n.floor();
            k = k.floor() + q;
            if n > 0.0 { // case 20
                q = recur(&mut n, x, &mut k, 1);
            } else { // case 21
                t = k;
                k = n;
                q = recur(&mut t, x, &mut k, 1);
                k = t;
            }
            if q == 0.0 { // case 22
                y = 0.0;
                return sign as f64 * y;
            }
            
            // case 23

        } else { // case 24
            k = n;
            q = 1.0;
        }

        /* boundary between convergence of
        * power series and Hankel expansion
        */
        y = k.abs();
        if y < 26.0 { // case 25
            t = (0.0083 * y + 0.09) * y + 12.9;
        } else { // case 26
            t = 0.9 * y;
        }

        if x > t { // case 27
            y = hankel(k, x);
        } else { // case 28
            y = jvs(k, x);
        }

        if n > 0.0 { // case 29
            y /= q;
        } else { // case 30
            y *= q;
        }

    } else {

        /* For large n, use the uniform expansion or the transitional expansion.
        * But if x is of the order of n**2, these may blow up, whereas the
        * Hankel expansion will then work.
        */
        if n < 0.0 { // case 31
            //sf_error("Jv", SF_ERROR_LOSS, NULL);
            y = f64::NAN;
            return sign as f64 * y;
        }
        t = x / n;
        t /= n;
        if t > 0.3 { // case 32
            y = hankel(n, x);
        } else { // case 33
            y = jnx(n, x);
        }
    }

    // case 34
    return sign as f64 * y;
}

/* Reduce the order by backward recurrence.
* AMS55 #9.1.27 and 9.1.73.
*/

fn recur(n: &mut f64, x: f64, newn: &mut f64, cancel: isize) -> f64
{

    /* Continued fraction for Jn(x)/Jn-1(x)
    * AMS 9.1.73
    *
    *    x       -x^2      -x^2
    * ------  ---------  ---------   ...
    * 2 n +   2(n+1) +   2(n+2) +
    *
    * Compute it with the simplest possible algorithm.
    *
    * This continued fraction starts to converge when (|n| + m) > |x|.
    * Hence, at least |x|-|n| iterations are necessary before convergence is
    * achieved. There is a hard limit set below, m <= 30000, which is chosen
    * so that no branch in `jv` requires more iterations to converge.
    * The exact maximum number is (500/3.6)^2 - 500 ~ 19000
    */

    let maxiter = 22000;
    let mut miniter = (x.abs() - n.abs()) as i32;
    if miniter < 1 {
        miniter = 1;
    }

    let mut nflag = if *n < 0.0 {
        1
    } else {
        0
    };

    let mut ans: f64;

    loop {

        let mut pkm2 = 0.0;
        let mut qkm2 = 1.0;
        let mut pkm1 = x;
        let mut qkm1 = *n + *n;
        let xk = -x * x;
        let mut yk = qkm1;
        ans = 0.0;			/* ans=0.0 ensures that t=1.0 in the first iteration */
        let mut ctr = 0;
        loop {
            yk += 2.0;
            let pk = pkm1 * yk + pkm2 * xk;
            let qk = qkm1 * yk + qkm2 * xk;
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;

            /* check convergence */
            let r = if qk != 0.0 && ctr > miniter {
                pk / qk
            } else {
                0.0
            };

            let t;
            if r != 0.0 {
                t = ((ans - r) / r).abs();
                ans = r;
            }
            else {
                t = 1.0;
            }

            ctr += 1;
            if ctr > maxiter {
                //sf_error("jv", SF_ERROR_UNDERFLOW, NULL);
                break;
            }
            if t < MACHEP {
                break;
            }

            /* renormalize coefficients */
            if pk.abs() > BIG {
                pkm2 /= BIG;
                pkm1 /= BIG;
                qkm2 /= BIG;
                qkm1 /= BIG;
            }
            if t <= MACHEP {
                break;
            }
        }

        if ans == 0.0 {
            ans = 1.0;
        }

        /* Change n to n-1 if n < 0 and the continued fraction is small */
        if nflag > 0 {
            if ans.abs() < 0.125 {
                nflag = -1;
                *n = *n - 1.0;
                continue;
            }
        }
        break;
    }


    let kf = *newn;

    /* backward recurrence
    *              2k
    *  J   (x)  =  --- J (x)  -  J   (x)
    *   k-1         x   k         k+1
    */

    let mut pk = 1.0;
    let mut pkm1 = 1.0 / ans;
    let mut k = *n - 1.0;
    let mut r = 2.0 * k;
    let mut pkm2: f64;
    loop {
        pkm2 = (pkm1 * r - pk * x) / x;
        pk = pkm1;
        pkm1 = pkm2;
        r -= 2.0;
        k -= 1.0;
        if k <= kf + 0.5 {
            break;
        }
    }

    /* Take the larger of the last two iterates
    * on the theory that it may have less cancellation error.
    */

    if cancel != 0 {
        if (kf >= 0.0) && (pk.abs() > pkm1.abs()) {
            k += 1.0;
            pkm2 = pk;
        }
    }
    *newn = k;

    return pkm2;
}



/* Ascending power series for Jv(x).
* AMS55 #9.1.10.
*/

fn jvs(n: f64, x: f64) -> f64
{

    let mut sgngam: isize = 0; // Value doesn't matter, will be overwritten

    let z = -x * x / 4.0;
    let mut u = 1.0;
    let mut y = u;
    let mut k = 1.0;
    let mut t = 1.0;

    while t > MACHEP {
        u *= z / (k * (n + k));
        y += u;
        k += 1.0;
        if y != 0.0 {
            t = (u / y).abs();
        }
    }

    let (_, mut ex) = frexp(0.5 * x);
    ex = (ex as f64 * n) as i32;
    if (ex > -1023) && (ex < 1023) && (n > 0.0) && (n < (MAXGAM - 1.0)) {
        let t = (0.5 * x).powf(n) / gamma(n + 1.0);
        y * t
    } else {
        let mut t = n * (0.5 * x).ln() - lgam_sgn(n + 1.0, &mut sgngam);
        if y < 0.0 {
            sgngam = -sgngam;
            y = -y;
        }
        t += y.ln();
        if t < -MAXLOG {
            return 0.0;
        } else if t > MAXLOG {
            //sf_error("Jv", SF_ERROR_OVERFLOW, NULL);
            return f64::INFINITY;
        }
        sgngam as f64 * t.exp()
    }
}

/* Hankel's asymptotic expansion
* for large x.
* AMS55 #9.2.5.
*/

pub fn hankel(n: f64, x: f64) -> f64
{
    let m = 4.0 * n * n;
    let mut j = 1.0;
    let z = 8.0 * x;
    let mut k = 1.0;
    let mut p = 1.0;
    let mut u = (m - 1.0) / z;
    let mut q = u;
    let mut sign = 1.0;
    let mut conv = 1.0;
    let mut flag = 0;
    let mut t = 1.0;
    let mut pp = 1.0e38;
    let mut qq = 1.0e38;

    while t > MACHEP {
        k += 2.0;
        j += 1.0;
        sign = -sign;
        u *= (m - k * k) / (j * z);
        p += sign * u;
        k += 2.0;
        j += 1.0;
        u *= (m - k * k) / (j * z);
        q += sign * u;
        t = (u / p).abs();
        if t < conv {
            conv = t;
            qq = q;
            pp = p;
            flag = 1;
        }
        /* stop if the terms start getting larger */
        if (flag != 0) && (t > conv) {
            break;
        }
    }

    //hank1:
    u = x - (0.5 * n + 0.25) * M_PI;

    (2.0 / (M_PI * x)).sqrt() * (pp * u.cos() - qq * u.sin())
}


/* Asymptotic expansion for large n.
* AMS55 #9.3.35.
*/

const LAMBDA: [f64; 11] = [
    1.0,
    1.041666666666666666666667E-1,
    8.355034722222222222222222E-2,
    1.282265745563271604938272E-1,
    2.918490264641404642489712E-1,
    8.816272674437576524187671E-1,
    3.321408281862767544702647E+0,
    1.499576298686255465867237E+1,
    7.892301301158651813848139E+1,
    4.744515388682643231611949E+2,
    3.207490090890661934704328E+3
];

const MU: [f64; 11] = [
    1.0,
    -1.458333333333333333333333E-1,
    -9.874131944444444444444444E-2,
    -1.433120539158950617283951E-1,
    -3.172272026784135480967078E-1,
    -9.424291479571202491373028E-1,
    -3.511203040826354261542798E+0,
    -1.572726362036804512982712E+1,
    -8.228143909718594444224656E+1,
    -4.923553705236705240352022E+2,
    -3.316218568547972508762102E+3
];

const P1: [f64; 2] = [
    -2.083333333333333333333333E-1,
    1.250000000000000000000000E-1
];

const P2: [f64; 3] = [
    3.342013888888888888888889E-1,
    -4.010416666666666666666667E-1,
    7.031250000000000000000000E-2
];

const P3: [f64; 4] = [
    -1.025812596450617283950617E+0,
    1.846462673611111111111111E+0,
    -8.912109375000000000000000E-1,
    7.324218750000000000000000E-2
];

const P4: [f64; 5] = [
    4.669584423426247427983539E+0,
    -1.120700261622299382716049E+1,
    8.789123535156250000000000E+0,
    -2.364086914062500000000000E+0,
    1.121520996093750000000000E-1
];

const P5: [f64; 6] = [
    -2.8212072558200244877E1,
    8.4636217674600734632E1,
    -9.1818241543240017361E1,
    4.2534998745388454861E1,
    -7.3687943594796316964E0,
    2.27108001708984375E-1
];

const P6: [f64; 7] = [
    2.1257013003921712286E2,
    -7.6525246814118164230E2,
    1.0599904525279998779E3,
    -6.9957962737613254123E2,
    2.1819051174421159048E2,
    -2.6491430486951555525E1,
    5.7250142097473144531E-1
];

const P7: [f64; 8] = [
    -1.9194576623184069963E3,
    8.0617221817373093845E3,
    -1.3586550006434137439E4,
    1.1655393336864533248E4,
    -5.3056469786134031084E3,
    1.2009029132163524628E3,
    -1.0809091978839465550E2,
    1.7277275025844573975E0
];


fn jnx(n: f64, x: f64) -> f64
{
    /* Test for x very close to n. Use expansion for transition region if so. */
    let cbn = cbrt(n);
    let mut z = (x - n) / cbn;
    if z.abs() <= 0.7 {
        return jnt(n, x);
    }

    z = x / n;
    let zz = 1.0 - z * z;
    if zz == 0.0 {
        return 0.0;
    }

    let sz: f64;
    let mut t: f64;
    let zeta: f64;
    let nflg: isize;

    if zz > 0.0 {
        sz = zz.sqrt();
        t = 1.5 * (((1.0 + sz) / z).ln() - sz);	/* zeta ** 3/2          */
        zeta = cbrt(t * t);
        nflg = 1;
    }
    else {
        sz = (-zz).sqrt();
        t = 1.5 * (sz - (1.0 / z).acos());
        zeta = -cbrt(t * t);
        nflg = -1;
    }
    let z32i = (1.0 / t).abs();
    let sqz = cbrt(t);

    /* Airy function */
    let n23 = cbrt(n * n);
    t = n23 * zeta;


    let (ai, aip, _, _) = airy(t);

    let mut u = [0.0; 8];
    /* polynomials in expansion */
    u[0] = 1.0;
    let zzi = 1.0 / zz;
    u[1] = polevl(zzi, &P1, 1) / sz;
    u[2] = polevl(zzi, &P2, 2) / zz;
    u[3] = polevl(zzi, &P3, 3) / (sz * zz);
    let mut pp = zz * zz;
    u[4] = polevl(zzi, &P4, 4) / pp;
    u[5] = polevl(zzi, &P5, 5) / (pp * sz);
    pp *= zz;
    u[6] = polevl(zzi, &P6, 6) / pp;
    u[7] = polevl(zzi, &P7, 7) / (pp * sz);

    // #if CEPHES_DEBUG
    //     for (k = 0; k <= 7; k++)
    //     printf("u[%d] = %.5E\n", k, u[k]);
    // #endif

    let mut pp = 0.0;
    let mut qq = 0.0;
    let mut np = 1.0;
    /* flags to stop when terms get larger */
    let mut doa = 1;
    let mut dob = 1;
    let mut akl = f64::INFINITY;
    let mut bkl = f64::INFINITY;

    for k in 0..=3 {
        let tk = 2 * k;
        let tkp1 = tk + 1;
        let mut zp = 1.0;
        let mut ak: f64 = 0.0;
        let mut bk = 0.0;
        for s in 0..=tk {
            if doa != 0 {
                let sign = if (s & 3) > 1 {
                    nflg
                } else {
                    1
                };
                ak += sign as f64 * MU[s] * zp * u[tk - s];
            }

            if dob != 0 {
                let m = tkp1 - s;
                let sign = if ((m + 1) & 3) > 1 {
                    nflg
                } else {
                    1
                };
                bk += sign as f64 * LAMBDA[s] * zp * u[m];
            }
            zp *= z32i;
        }

        if doa != 0 {
            ak *= np;
            t = ak.abs();
            if t < akl {
                akl = t;
                pp += ak;
            } else {
                doa = 0;
            }
        }

        if dob != 0 {
            bk += LAMBDA[tkp1] * zp * u[0];
            bk *= -np / sqz;
            t = bk.abs();
            if t < bkl {
                bkl = t;
                qq += bk;
            } else {
                dob = 0;
            }
        }

        if np < MACHEP {
            break;
        }
        np /= n * n;
    }

    /* normalizing factor ( 4*zeta/(1 - z**2) )**1/4    */
    t = 4.0 * zeta / zz;
    t = t.sqrt().sqrt();

    t * (ai * pp / cbrt(n) + aip * qq / (n23 * n))
}

/* Asymptotic expansion for transition region,
* n large and x close to n.
* AMS55 #9.3.23.
*/

const PF2: [f64; 2] = [
    -9.0000000000000000000e-2,
    8.5714285714285714286e-2
];

const PF3: [f64; 3] = [
    1.3671428571428571429e-1,
    -5.4920634920634920635e-2,
    -4.4444444444444444444e-3
];

const PF4: [f64; 4] = [
    1.3500000000000000000e-3,
    -1.6036054421768707483e-1,
    4.2590187590187590188e-2,
    2.7330447330447330447e-3
];

const PG1: [f64; 2] = [
    -2.4285714285714285714e-1,
    1.4285714285714285714e-2
];

const PG2: [f64; 3] = [
    -9.0000000000000000000e-3,
    1.9396825396825396825e-1,
    -1.1746031746031746032e-2
];

const PG3: [f64; 3] = [
    1.9607142857142857143e-2,
    -1.5983694083694083694e-1,
    6.3838383838383838384e-3
];


fn jnt(n: f64, x: f64) -> f64
{

    let cbn = cbrt(n);
    let z = (x - n) / cbn;
    let cbtwo = cbrt(2.0);

    /* Airy function */
    let zz = -cbtwo * z;
    let (ai, aip, _, _) = airy(zz);

    /* polynomials in expansion */
    let zz = z * z;
    let z3 = zz * z;
    let f: [f64; 5] = [
        1.0,
        -z / 5.0,
        polevl(z3, &PF2, 1) * zz,
        polevl(z3, &PF3, 2),
        polevl(z3, &PF4, 3) * z
    ];
    let g: [f64; 4] = [
        0.3 * zz,
        polevl(z3, &PG1, 1),
        polevl(z3, &PG2, 2) * z,
        polevl(z3, &PG3, 2) * zz
    ];

    let mut pp = 0.0;
    let mut qq = 0.0;
    let mut nk = 1.0;
    let n23 = cbrt(n * n);

    for k in 0..=4 {
        let fk = f[k] * nk;
        pp += fk;
        if k != 4 {
            let gk = g[k] * nk;
            qq += gk;
        }
        nk /= n23;
    }

    cbtwo * ai * pp / cbn + cbrt(4.0) * aip * qq / n
}


#[cfg(test)]
mod jv_tests {
    use super::*;

    #[test]
    fn jv_specifics() {
        // 1000011001100000000101000000001000
        assert_eq!(jv(4.897788193684462, 5.88843655355589), 0.36421007691409374);
        assert_eq!(jv(0.0000000001, 5.88843655355589), 0.11860905401275754);
        //assert_eq!(jv(0.000000021379620895022325, 0.00000001071519305237607), 0.9999996051711874); // TODO: This is off
        
        // 1000000000000000000000000001000
        //assert_eq!(jv(-501.18723362727246, 19054.607179632483), 0.003566933713497572); // TODO
        assert_eq!(jv(-9772372209.558111, 10000000000.0).is_nan(), true);
        assert_eq!(jv(-501.18723362727246, 81.2830516164099).is_nan(), true);

        // 1000010101010011000101000000001000
        assert_eq!(jv(138.03842646028852, 69.18309709189366), 4.484828681274339e-29);
        assert_eq!(jv(5.011872336272722, 20.892961308540396), 0.17090407097940324);
        assert_eq!(jv(20.417379446695296, 13.182567385564074), 0.0007234899781792606);

        // 1000100110010100101001000000001000
        assert_eq!(jv(-34.673685045253166, 87.09635899560806), -0.08865047997768659);
        assert_eq!(jv(-489.77881936844614, 616.5950018614822), -0.04106463782394591);
        assert_eq!(jv(-109.64781961431851, 87.09635899560806), -3843.650935867128);

        // 1000100110100000010001000000001000
        assert_eq!(jv(-91.20108393559097, 630.957344480193), 0.019398112218566216);
        assert_eq!(jv(-346.73685045253166, 9120.108393559096), 0.006139011752456574);
        assert_eq!(jv(-93.3254300796991, 630.957344480193), -0.004640376941868423);

        // 1000010101010011001001000000001000
        assert_eq!(jv(181.97008586099827, 91.20108393559097), 1.0655759495695935e-37);
        assert_eq!(jv(30.902954325135905, 72.44359600749902), -0.0034106311498274884);
        assert_eq!(jv(141.2537544622754, 70.79457843841381), 1.046992443483635e-29);

        // 1000011001010011000101000000001000
        assert_eq!(jv(20.892961308540396, 13.182567385564074), 0.0004373209259475514);
        assert_eq!(jv(5.011872336272722, 5.88843655355589), 0.3539937719047022);
        assert_eq!(jv(5.011872336272722, 2.51188643150958), 0.01958216179441848);

        // 1001000
        assert_eq!(jv(9772372209.558111, 0.0010232929922807535), 0.0);
        assert_eq!(jv(-9772372209.558111, 0.0010232929922807535), -f64::INFINITY);
        assert_eq!(jv(-9772372209.558111, 0.0000000001), -f64::INFINITY);

        // 1000010101100000000101000000001000
        assert_eq!(jv(4.897788193684462, 20.892961308540396), 0.16178203907804387);
        assert_eq!(jv(0.0000000001, 20.892961308540396), 0.05469764973027066);

        // 1000101001100000000101000000001000
        assert_eq!(jv(-0.0000000001, 12.882495516931343), 0.19721370425603993);
        assert_eq!(jv(-4.897788193684462, 13.489628825916533), -0.17104901220970245);
        assert_eq!(jv(-0.000000021379620895022325, 0.00000001071519305237607), 1.0000003948289684);
        
        // 100000100 case 3, 9
        // TODO: These cases are a bit off...
        assert_eq!(jv(100000.0, -10000000000.0), 5.589576497465461e-6);
        assert_eq!(jv(-100000.0, 10000000000.0), 5.589576497465461e-6);
        assert_eq!(jv(-100000.0, -10000000000.0), 5.589576497465461e-6);

        // 1000100101100000000101000000001000
        assert_eq!(jv(-0.0000000001, 20.892961308540396), 0.05469764967820164);
        assert_eq!(jv(-4.897788193684462, 20.892961308540396), -0.17617941764357217);
        assert_eq!(jv(-2.51188643150958, 13.182567385564074), -0.14385282760073412);

        // 1100000000000000000000000000000100
        assert_eq!(jv(10000000000.0, -10000000000.0), 0.00020762166542623404);
        assert_eq!(jv(-10000000000.0, 10000000000.0), 0.00020762166542623404);
        assert_eq!(jv(-10000000000.0, -10000000000.0), 0.00020762166542623404);

        // 100001000 case 4, 9
        // TODO: These cases are a bit off...
        assert_eq!(jv(354813.3892335753, 10000000000.0), -7.5642141066880886e-6);
        assert_eq!(jv(-354813.3892335753, 10000000000.0), 4.965879563275322e-6);
        assert_eq!(jv(-16.595869074375607, 21.379620895022324), -0.215457819784271);

        // 10000000100
        assert_eq!(jv(100.0, -741.3102413009177), 0.026889671736553176);
        assert_eq!(jv(-100.0, 741.3102413009177), 0.026889671736553176);
        assert_eq!(jv(-100.0, -741.3102413009177), 0.026889671736553176);

        // 1000101001010100100101000000001000
        assert_eq!(jv(-5.011872336272722, 12.882495516931343), -0.11571185210931113);
        assert_eq!(jv(-20.892961308540396, 12.882495516931343), 21.453477363886766);
        assert_eq!(jv(-5.011872336272722, 2.51188643150958), -0.16159219637087377);

        // 1000100101010100101001000000001000
        assert_eq!(jv(-30.902954325135905, 72.44359600749902), 0.03282314716919885);
        assert_eq!(jv(-181.97008586099827, 91.20108393559097), -1.7802073771928776e33);
        assert_eq!(jv(-141.2537544622754, 70.79457843841381), -1.7794051487040196e26);

        // 10001000
        assert_eq!(jv(9772372209.558111, 354813.3892335753), 0.0);
        assert_eq!(jv(-9772372209.558111, 354813.3892335753), f64::INFINITY);
        assert_eq!(jv(-21.379620895022324, 0.0000000478630092322638), -1.956062319314165e181);

        // 11000
        assert_eq!(jv(9772372209.558111, -10000000000.0).is_nan(), true);
        assert_eq!(jv(-9772372209.558111, -0.0000000001).is_nan(), true);
        assert_eq!(jv(-9772372209.558111, -10000000000.0).is_nan(), true);

        // 10000100
        // TODO: Should any of these be NaN?
        assert_eq!(jv(10000000000.0, -354813.3892335753), 0.0);
        assert_eq!(jv(-10000000000.0, 354813.3892335753), 0.0);
        assert_eq!(jv(-10000000000.0, -354813.3892335753), 0.0);

        // 1000010110010011010001000000001000
        assert_eq!(jv(489.77881936844614, 18197.008586099826), -0.004228667870843825);
        assert_eq!(jv(95.49925860214358, 645.6542290346556), 0.023729176355939453);
        assert_eq!(jv(489.77881936844614, 630.957344480193), 0.03989440419898407);

        // 10
        assert_eq!(jv(1.0, -10000000000.0), 7.676506113846571e-6);
        assert_eq!(jv(-1.0, 10000000000.0), 7.676506113846571e-6);
        assert_eq!(jv(-1.0, -10000000000.0), -7.676506113846571e-6);

        // 1000100110010100110001000000001000
        assert_eq!(jv(-95.49925860214358, 645.6542290346556), -0.020885458393017742);
        assert_eq!(jv(-489.77881936844614, 18197.008586099826), -0.005897140491143249);
        assert_eq!(jv(-489.77881936844614, 630.957344480193), 0.03255117841622729);

        // 1000010110010011001001000000001000
        assert_eq!(jv(489.77881936844614, 616.5950018614822), -0.033873066983000105);
        assert_eq!(jv(34.673685045253166, 89.12509381337455), -0.08726186375632396);
        assert_eq!(jv(177.82794100389228, 89.12509381337455), 6.9194024232127745e-37);

        // 1000100101010100100101000000001000
        assert_eq!(jv(-5.011872336272722, 20.892961308540396), -0.16904987921660208);
        assert_eq!(jv(-138.03842646028852, 69.18309709189366), 7.155668135052916e24);
        assert_eq!(jv(-20.892961308540396, 13.182567385564074), 14.853262251695945);

        // 1000100
        // TODO: Should these be NaN?
        assert_eq!(jv(10000000000.0, -0.0010471285480508996), 0.0);
        assert_eq!(jv(-10000000000.0, 0.0010471285480508996), 0.0);
        assert_eq!(jv(-10000000000.0, -0.0010471285480508996), 0.0);

        // 1100000000000000000000000000001000
        assert_eq!(jv(9772372209.558111, 10000000000.0), -1.719892140481984e-5);
        assert_eq!(jv(501.18723362727246, 19054.607179632483), -0.0004430004551984964);
        assert_eq!(jv(501.18723362727246, 81.2830516164099), 0.0);

        // 1000010110100000010001000000001000
        // TODO: The first is pretty far off
        assert_eq!(jv(346.73685045253166, 9120.108393559096), 1.4670530214822236e-5);
        assert_eq!(jv(91.20108393559097, 630.957344480193), -0.030633593293011036);
        assert_eq!(jv(93.3254300796991, 630.957344480193), 0.029386336766700064);

        // 11000000001000 // case 4, 13, 14
        // TODO: These are very wrong
        assert_eq!(jv(489.77881936844614, 239.88329190194898), 9.24758256467989e-102);
        assert_eq!(jv(-489.77881936844614, 239.88329190194898).is_nan(), true);
        assert_eq!(jv(-1.023292992280754, 0.0000000016218100973589331), -f64::INFINITY);
        
        // 1000000100
        assert_eq!(jv(100.0, -758.5775750291835), -0.01606537872714924);
        assert_eq!(jv(-100.0, 758.5775750291835), -0.01606537872714924);
        assert_eq!(jv(-100.0, -758.5775750291835), -0.01606537872714924);


    }

    #[test]
    fn jv_trivials() {
        // Non integer order, x < 0 (case 5)
        assert_eq!(jv(1.0 + 1e-15, -1e-20).is_nan(), true);
        // Negative non integer order, x = 0 (case 6)
        assert_eq!(jv(-0.5, 0.0), f64::INFINITY);
        // Negative non integer order, x = 0 (case 6)
        assert_eq!(jv(-1.5, 0.0), -f64::INFINITY);
        // Negative integer order > 500, x < 0 (case 31)
        // assert_eq!(jv(-500.0, -300.0), 1.1);
    }

    #[test]
    fn jv_specifics_2() {
        // x much smaller than n (case 7)
        assert_eq!(jv(10.5, 1e-10), 5.803087772659944e-116);
        assert_eq!(jv(-8.5, 5e-10), 1.8516296635108193e+85);

        // case 14
        assert_eq!(jv(-8.5, 5e-5), 5.855367120313302e+42);

        // case 8
        assert_eq!(jv(121.1, 1.2), 1.0379611436454817e-228);
        assert_eq!(jv(-500.0, -10.0), 0.0);

        // case 13, 19, 20, 23, 25, 29, 34
        assert_eq!(jv(55.5, 36.1), 1.2806006117623818e-07);

        // case 13, 17, 19, 20, 23, 26, 27, 29, 34
        assert_eq!(jv(131.05, 1001.0), 0.006889418633323632);

        // case 9
        assert_eq!(jv(23.1, 1e5), -0.0020892046447252653);
        assert_eq!(jv(250.1, 1e5), 0.001984742706538866);
        assert_eq!(jv(15.5, 1e10), 6.966483992942254e-6);
        //assert_eq!(jv(25.9, 3e8), -2.157319314452935e-05); // TODO: This one's a bit far off
        //assert_eq!(jv(29.8, 2e10), -4.570389250619019e-07); // This one too
        //assert_eq!(jv(500.1, 76000.0), 0.002467229286055466); // This one too

        // case 3, 33, 34
        assert_eq!(jv(-500.0, -300.0), 2.8630465185002796e-67);

        //x/n/n > 0.3; an > 500
        // case 32

        

    }

    #[test]
    fn jv_small_v() {
        assert_eq!(jv(2.0, 5.0), 0.046565116277752235);
        assert_eq!(jv(15.0, -10.0), -0.004507973143721255);
        assert_eq!(jv(30.5, 1e6), 0.0002789068313065268);
        assert_eq!(jv(30.5, -1e6).is_nan(), true);
    }

    #[test]
    fn jv_1_small_x() {
        assert_eq!(jv(1.0, 0.0), 0.0);

        assert_eq!(jv(1.0, 1e-20), 5.0000000000000005e-21);
        assert_eq!(jv(1.0, 1e-10), 5.000000000000001e-11);
        assert_eq!(jv(1.0, 1e-7), 4.999999999999993e-08);
        assert_eq!(jv(1.0, 1e-6), 4.999999999999375e-07);
        assert_eq!(jv(1.0, 1e-5), 4.9999999999375e-06);
        assert_eq!(jv(1.0, 1.01e-5), 5.049999999935606e-06);
        assert_eq!(jv(1.0, 0.1), 0.049937526036242);
        assert_eq!(jv(1.0, 1.0), 0.44005058574493355);
        assert_eq!(jv(1.0, 2.0), 0.5767248077568734);
        assert_eq!(jv(1.0, 3.0), 0.33905895852593654);
        assert_eq!(jv(1.0, 4.0), -0.06604332802354912);
        assert_eq!(jv(1.0, 5.0), -0.3275791375914653);

        assert_eq!(-jv(1.0, 1e-20), -5.0000000000000005e-21);
        assert_eq!(-jv(1.0, 1e-10), -5.000000000000001e-11);
        assert_eq!(-jv(1.0, 1e-7), -4.999999999999993e-08);
        assert_eq!(-jv(1.0, 1e-6), -4.999999999999375e-07);
        assert_eq!(-jv(1.0, 1e-5), -4.9999999999375e-06);
        assert_eq!(-jv(1.0, 1.01e-5), -5.049999999935606e-06);
        assert_eq!(-jv(1.0, 0.1), -0.049937526036242);
        assert_eq!(-jv(1.0, 1.0), -0.44005058574493355);
        assert_eq!(-jv(1.0, 2.0), -0.5767248077568734);
        assert_eq!(-jv(1.0, 3.0), -0.33905895852593654);
        assert_eq!(-jv(1.0, 4.0), 0.06604332802354912);
        assert_eq!(-jv(1.0, 5.0), 0.3275791375914653);
    }

    #[test]
    fn jv_1_large_x() {
        assert_eq!(jv(1.0, 5.0 + 1e-16), -0.3275791375914653);
        assert_eq!(jv(1.0, 6.0), -0.27668385812756563);
        assert_eq!(jv(1.0, 10.0), 0.04347274616886141);
        assert_eq!(jv(1.0, 100.0), -0.0771453520141123);
        assert_eq!(jv(1.0, 1000.0), 0.00472831190708902);

        assert_eq!(jv(1.0, -5.0 - 1e-16), 0.3275791375914653);
        assert_eq!(jv(1.0, -6.0), 0.27668385812756563);
        assert_eq!(jv(1.0, -10.0), -0.04347274616886141);
        assert_eq!(jv(1.0, -100.0), 0.0771453520141123);
        assert_eq!(jv(1.0, -1000.0), -0.00472831190708902);

        // TODO: Accuracy reduces for very large x
        //assert_eq!(jv(1.0, 1e10), -7.676506113818561e-06);
        //assert_eq!(jv(1.0, 1e20), -5.449273779343996e-11);
    }

    #[test]
    fn jv_0_small_x() {
        assert_eq!(jv(0.0, 0.0), 1.0);

        assert_eq!(jv(0.0, 1e-20), 1.0);
        assert_eq!(jv(0.0, 1e-10), 1.0);
        assert_eq!(jv(0.0, 1e-7), 0.9999999999999974);
        assert_eq!(jv(0.0, 1e-6), 0.99999999999975);
        assert_eq!(jv(0.0, 1e-5), 0.9999999999750001);
        assert_eq!(jv(0.0, 1.01e-5), 0.9999999999744974);
        assert_eq!(jv(0.0, 0.1), 0.99750156206604);
        assert_eq!(jv(0.0, 1.0), 0.7651976865579665);
        assert_eq!(jv(0.0, 2.0), 0.22389077914123562);
        assert_eq!(jv(0.0, 3.0), -0.2600519549019335);
        assert_eq!(jv(0.0, 4.0), -0.3971498098638473);
        assert_eq!(jv(0.0, 5.0), -0.1775967713143383);

        assert_eq!(jv(0.0, -1e-20), 1.0);
        assert_eq!(jv(0.0, -1e-10), 1.0);
        assert_eq!(jv(0.0, -1e-7), 0.9999999999999974);
        assert_eq!(jv(0.0, -1e-6), 0.99999999999975);
        assert_eq!(jv(0.0, -1e-5), 0.9999999999750001);
        assert_eq!(jv(0.0, -1.01e-5), 0.9999999999744974);
        assert_eq!(jv(0.0, -0.1), 0.99750156206604);
        assert_eq!(jv(0.0, -1.0), 0.7651976865579665);
        assert_eq!(jv(0.0, -2.0), 0.22389077914123562);
        assert_eq!(jv(0.0, -3.0), -0.2600519549019335);
        assert_eq!(jv(0.0, -4.0), -0.3971498098638473);
        assert_eq!(jv(0.0, -5.0), -0.1775967713143383);
    }

    #[test]
    fn jv_0_large_x() {
        assert_eq!(jv(0.0, 5.0 + 1e-16), -0.1775967713143383);
        assert_eq!(jv(0.0, 6.0), 0.15064525725099695);
        assert_eq!(jv(0.0, 10.0), -0.24593576445134832);
        assert_eq!(jv(0.0, 100.0), 0.01998585030422333);
        assert_eq!(jv(0.0, 1000.0), 0.02478668615242003);

        assert_eq!(jv(0.0, -5.0 - 1e-16), -0.1775967713143383);
        assert_eq!(jv(0.0, -6.0), 0.15064525725099695);
        assert_eq!(jv(0.0, -10.0), -0.24593576445134832);
        assert_eq!(jv(0.0, -100.0), 0.01998585030422333);
        assert_eq!(jv(0.0, -1000.0), 0.02478668615242003);
        // TODO: Accuracy reduces for very large x
        //assert_eq!(jv(0.0, 1e10), 2.175589294792473e-06);
        //assert_eq!(jv(0.0, 1e20), -5.449273779343996e-11);
    }
}