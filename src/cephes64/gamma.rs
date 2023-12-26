/*
* Cephes Math Library Release 2.2:  July, 1992
* Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

const P: [f64; 7] = [
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1
];

const Q: [f64; 8] = [
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0
];

const MAXGAM: f64 = 171.624376956302725;
const LOGPI: f64 = 1.14472988584940017414;

/* Stirling's formula for the Gamma function */
const STIR: [f64; 5] = [
    7.87311395793093628397E-4,
    -2.29549961613378126380E-4,
    -2.68132617805781232825E-3,
    3.47222221605458667310E-3,
    8.33333333333482257126E-2,
];

const MAXSTIR: f64 = 143.01608;
const SQTPI: f64 = 2.50662827463100050242E0;

use crate::cephes64::consts::{M_PI, MAXLOG};
use crate::cephes64::polevl::{polevl, p1evl};

/* Gamma function computed by Stirling's formula.
* The polynomial STIR is valid for 33 <= x <= 172.
*/
fn stirf(x: f64) -> f64
{

    if x >= MAXGAM {
        return f64::INFINITY;
    }
    let w = 1.0 / x;
    let w = 1.0 + w * polevl(w, &STIR, 4);
    let y = x.exp();
    let y = if x > MAXSTIR {		
        /* Avoid overflow in pow() */
        let v = x.powf(0.5 * x - 0.25);
        v * (v / y)
    }
    else {
        x.powf(x - 0.5) / y
    };
    SQTPI * y * w
}


pub fn gamma(x: f64) -> f64 {
    //! Gamma function
    //! 
    //! ## DESCRIPTION:
    //! 
    //! Returns Gamma function of the argument.  The result is
    //! correctly signed.
    //!
    //! Arguments |x| <= 34 are reduced by recurrence and the function
    //! approximated by a rational function of degree 6/7 in the
    //! interval (2, 3).  Large arguments are handled by Stirling's
    //! formula. Large negative arguments are made positive using
    //! a reflection formula.
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
    //!     <td>-170, -33</td>
    //!     <td>20000</td>
    //!     <td>2.3e-15</td>
    //!     <td>3.3e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>-33,  33</td>
    //!     <td>20000</td>
    //!     <td>9.4e-16</td>
    //!     <td>2.2e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>33, 171.6</td>
    //!     <td>20000</td>
    //!     <td>2.3e-15</td>
    //!     <td>3.2e-16</td>
    //! </tr>
    //!</table>
    //!
    //! Error for arguments outside the test range will be larger
    //! owing to error amplification by the exponential function.

    let mut sgngam = 1;

    if x.is_infinite() {
        return x;
    }
    let q = x.abs();

    if q > 33.0 {
        let z = if x < 0.0 {
            let mut p = q.floor();
            if p == q {
                // gamnan:
                // sf_error("Gamma", SF_ERROR_OVERFLOW, NULL);
                return f64::INFINITY;
            }
            if (p as isize & 1) == 0 {
                sgngam = -1;
            }
            let mut z = q - p;
            if z > 0.5 {
                p += 1.0;
                z = q - p;
            }
            z = q * (M_PI * z).sin();
            if z == 0.0 {
                return sgngam as f64 * f64::INFINITY;
            }
            z = z.abs();
            M_PI / (z * stirf(q))
        }
        else {
            stirf(x)
        };
        return sgngam as f64 * z;
    }

    let mut z = 1.0;
    let mut x = x;
    while x >= 3.0 {
        x -= 1.0;
        z *= x;
    }

    while x < 0.0 {
        if x > -1.0E-9 {
            if x == 0.0 {
                return f64::INFINITY;
            }
            else {
                return z / ((1.0 + 0.5772156649015329 * x) * x);
            }
        }
        z /= x;
        x += 1.0;
    }

    while x < 2.0 {
        if x < 1.0e-9 {
            if x == 0.0 {
                return f64::INFINITY;
            }
            else {
                return z / ((1.0 + 0.5772156649015329 * x) * x);
            }
        }
        z /= x;
        x += 1.0;
    }

    if x == 2.0 {
        return z;
    }

    x -= 2.0;
    let p = polevl(x, &P, 6);
    let q = polevl(x, &Q, 7);
    
    z * p / q
}



/* A[]: Stirling's formula expansion of log Gamma
* B[], C[]: log Gamma function between 2 and 3
*/
const A: [f64; 5] = [
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
];

const B: [f64; 6] = [
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
];

const C: [f64; 6] = [
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
];

/* log( sqrt( 2*pi ) ) */
const LS2PI: f64 = 0.91893853320467274178;

const MAXLGM: f64 = 2.556348e305;


pub fn lgam(x: f64) -> f64 {
    //! Natural logarithm of Gamma function
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns the base e (2.718...) logarithm of the absolute
    //! value of the Gamma function of the argument.
    //!
    //! For arguments greater than 13, the logarithm of the Gamma
    //! function is approximated by the logarithmic version of
    //! Stirling's formula using a polynomial approximation of
    //! degree 4. Arguments between -33 and +33 are reduced by
    //! recurrence to the interval [2,3] of a rational approximation.
    //! The cosecant reflection formula is employed for arguments
    //! less than -33.
    //!
    //! Arguments greater than MAXLGM return INFINITY and an error
    //! message.  MAXLGM = 2.556348e305 for IEEE arithmetic.
    //!
    //! ## ACCURACY:
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
    //!     <td>0, 3</td>
    //!     <td>28000</td>
    //!     <td>5.4e-16</td>
    //!     <td>1.1e-16</td>
    //! </tr>
    //! <tr>
    //!     <td>IEEE</td>
    //!     <td>2.718, 2.556e305</td>
    //!     <td>40000</td>
    //!     <td>3.5e-16</td>
    //!     <td>8.3e-17</td>
    //! </tr>
    //!</table>
    //!
    //! The error criterion was relative when the function magnitude
    //! was greater than one but absolute when it was less than one.
    //!
    //! The following test used the relative error criterion, though
    //! at certain points the relative error could be much higher than
    //! indicated.
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
    //!     <td>-200, -4</td>
    //!     <td>10000</td>
    //!     <td>4.8e-16</td>
    //!     <td>1.3e-16</td>
    //! </tr>
    //!</table>

    let mut sign = 0;
    lgam_sgn(x, &mut sign)
}

pub fn lgam_sgn(x: f64, sign: &mut isize) -> f64
{
    *sign = 1;

    if x.is_infinite() {
        return x;
    }

    if x < -34.0 {
        let q = -x;
        let w: f64;
        w = lgam_sgn(q, sign);
        let mut p = q.floor();
        if p == q {
        //lgsing:
            //sf_error("lgam", SF_ERROR_SINGULAR, NULL);
            return f64::INFINITY;
        }
        let i = p as isize;
        if (i & 1) == 0 {
            *sign = -1;
        } else {
            *sign = 1;
        }
        let mut z = q - p;
        if z > 0.5 {
            p += 1.0;
            z = p - q;
        }
        z = q * (M_PI * z).sin();
        if z == 0.0 {
            return f64::INFINITY;
        }
        /*     z = log(M_PI) - log( z ) - w; */
        z = LOGPI - z.ln() - w;
        return z;
    }
    let mut x = x;
    if x < 13.0 {
        let mut z = 1.0;
        let mut p = 0.0;
        let mut u = x;
        while u >= 3.0 {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while u < 2.0 {
            if u == 0.0 {
                return f64::INFINITY;
            }
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if z < 0.0 {
            *sign = -1;
            z = -z;
        }
        else {
            *sign = 1;
        }
        if u == 2.0 {
            return z.ln();
        }
        p -= 2.0;
        x = x + p;
        p = x * polevl(x, &B, 5) / p1evl(x, &C, 6);
        return z.ln() + p;
    }

    if x > MAXLGM {
        return (*sign as f64) * f64::INFINITY;
    }

    let mut q = (x - 0.5) * x.ln() - x + LS2PI;
    if x > 1.0e8 {
        return q;
    }
    

    let p = 1.0 / (x * x);
    if x >= 1000.0 {
        q += ((7.9365079365079365079365e-4 * p
            - 2.7777777777777777777778e-3) * p
            + 0.0833333333333333333333) / x;
    } else {
        q += polevl(p, &A, 4) / x;
    }

    q
}

#[cfg(test)]
mod gamma_tests {
    use super::*;

    #[test]
    fn gamma_trivials() {
        assert_eq!(gamma(f64::INFINITY), f64::INFINITY);
        //assert_eq!(gamma(-f64::INFINITY).is_nan(), true);
        assert_eq!(gamma(f64::NAN).is_nan(), true);
        assert_eq!(gamma(0.0), f64::INFINITY);
        assert_eq!(gamma(-10.0), f64::INFINITY);
        assert_eq!(gamma(-1e20), f64::INFINITY);
    }

    #[test]
    fn gamma_neg() {
        assert_eq!(gamma(-0.5), -3.5449077018110318);
        assert_eq!(gamma(-0.999999), -1000000.4227569911);
        assert_eq!(gamma(-10.0 + 1e-15), 155133915.73559144);
        assert_eq!(gamma(-100.5), -3.353690819807678e-159);
    }

    #[test]
    fn gamma_pos() {
        assert_eq!(gamma(0.5), 1.7724538509055159);
        assert_eq!(gamma(1.0), 1.0);
        assert_eq!(gamma(5.0), 24.0);
        assert_eq!(gamma(15.0), 87178291200.0);
        assert_eq!(gamma(33.1), 3.7275934243564036e+35);
        assert_eq!(gamma(40.0), 2.0397882081197442e+46);
        // TODO: Accuracy dwindles as x becomes large
        assert_eq!(gamma(50.5), 4.2904629123519605e+63);
        assert_eq!(gamma(100.0), 9.332621544394417e+155);
        assert_eq!(gamma(101.5), 9.36756791960313e+158);
    }
}

#[cfg(test)]
mod lgam_tests {
    use super::*;

    #[test]
    fn lgam_trivials() {
        assert_eq!(lgam(f64::INFINITY), f64::INFINITY);
        //assert_eq!(gamma(-f64::INFINITY).is_nan(), true);
        assert_eq!(lgam(f64::NAN).is_nan(), true);
        assert_eq!(lgam(0.0), f64::INFINITY);
        assert_eq!(lgam(-10.0), f64::INFINITY);
        assert_eq!(lgam(-1e20), f64::INFINITY);
    }

    #[test]
    fn lgam_very_large() { // 1e8 <= x <= MAXLGM
        assert_eq!(lgam(1e8), 1742068066.1038349);
        assert_eq!(lgam(5e30), 3.4843495351127736e+32);
        assert_eq!(lgam(2e100), 4.599033129599291e+102);
        assert_eq!(lgam(6e200), 2.7678526684082232e+203);
        assert_eq!(lgam(3e300), 2.0726224205606455e+303);
        assert_eq!(lgam(MAXLGM), 1.7951366714594412e+308);
    }

    #[test]
    fn lgam_large() { // 13 <= x < 1e8
        assert_eq!(lgam(13.0), 19.987214495661885);
        assert_eq!(lgam(15.5), 26.536914491115613);
        assert_eq!(lgam(20.9), 42.03380831296665);
        assert_eq!(lgam(100.1), 359.5942717888567);
        assert_eq!(lgam(2222.0), 14898.160014335313);
        assert_eq!(lgam(5e6), 72124735.55845618);
        assert_eq!(lgam(1e8 - 1e-8), 1742068066.1038344);
    }

    #[test]
    fn lgam_small() { // -34 <= x < 13
        assert_eq!(lgam(-34.0 + 1e-12), -60.95166994625487);
        assert_eq!(lgam(-15.0 - 1e-12), -0.2683391645460027);
        assert_eq!(lgam(-10.5), -15.147270590717842);
        assert_eq!(lgam(-0.1), 2.368961332728789);
        assert_eq!(lgam(-1e-20), 46.051701859880914);
        assert_eq!(lgam(1e-20), 46.051701859880914);
        assert_eq!(lgam(0.1), 2.252712651734206);
        assert_eq!(lgam(1.0), 0.0);
        assert_eq!(lgam(5.5), 3.9578139676187165);
        assert_eq!(lgam(12.5), 18.73434751193645);
        assert_eq!(lgam(13.0 - 1e-15), 19.98721449566188);
    }

    #[test]
    fn lgam_large_neg() { // x < -34
        assert_eq!(lgam(-34.0 - 1e-12), -60.95166994626196);
        assert_eq!(lgam(-100.5), -364.9009683094273);
        assert_eq!(lgam(-1e10 + 1e-6), -220258509298.6666);
    }
}