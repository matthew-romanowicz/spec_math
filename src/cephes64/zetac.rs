/*                                                     zetac.c
*
*     Riemann zeta function
*
*
*
* SYNOPSIS:
*
* double x, y, zetac();
*
* y = zetac( x );
*
*
*
* DESCRIPTION:
*
*
*
*                inf.
*                 -    -x
*   zetac(x)  =   >   k   ,   x > 1,
*                 -
*                k=2
*
* is related to the Riemann zeta function by
*
*     Riemann zeta(x) = zetac(x) + 1.
*
* Extension of the function definition for x < 1 is implemented.
* Zero is returned for x > log2(INFINITY).
*
* ACCURACY:
*
* Tabulated values have full machine accuracy.
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      1,50        10000       9.8e-16            1.3e-16
*
*
*/

/*
* Cephes Math Library Release 2.1:  January, 1989
* Copyright 1984, 1987, 1989 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Riemann zeta(x) - 1
* for integer arguments between 0 and 30.
*/
const AZETAC: [f64; 31] = [
    -1.50000000000000000000E0,
    0.0,  /* Not used; zetac(1.0) is infinity. */
    6.44934066848226436472E-1,
    2.02056903159594285400E-1,
    8.23232337111381915160E-2,
    3.69277551433699263314E-2,
    1.73430619844491397145E-2,
    8.34927738192282683980E-3,
    4.07735619794433937869E-3,
    2.00839282608221441785E-3,
    9.94575127818085337146E-4,
    4.94188604119464558702E-4,
    2.46086553308048298638E-4,
    1.22713347578489146752E-4,
    6.12481350587048292585E-5,
    3.05882363070204935517E-5,
    1.52822594086518717326E-5,
    7.63719763789976227360E-6,
    3.81729326499983985646E-6,
    1.90821271655393892566E-6,
    9.53962033872796113152E-7,
    4.76932986787806463117E-7,
    2.38450502727732990004E-7,
    1.19219925965311073068E-7,
    5.96081890512594796124E-8,
    2.98035035146522801861E-8,
    1.49015548283650412347E-8,
    7.45071178983542949198E-9,
    3.72533402478845705482E-9,
    1.86265972351304900640E-9,
    9.31327432419668182872E-10
];

/* 2**x (1 - 1/x) (zeta(x) - 1) = P(1/x)/Q(1/x), 1 <= x <= 10 */
const P: [f64; 9] = [
    5.85746514569725319540E11,
    2.57534127756102572888E11,
    4.87781159567948256438E10,
    5.15399538023885770696E9,
    3.41646073514754094281E8,
    1.60837006880656492731E7,
    5.92785467342109522998E5,
    1.51129169964938823117E4,
    2.01822444485997955865E2,
];

const Q: [f64; 8] = [
    /*  1.00000000000000000000E0, */
    3.90497676373371157516E11,
    5.22858235368272161797E10,
    5.64451517271280543351E9,
    3.39006746015350418834E8,
    1.79410371500126453702E7,
    5.66666825131384797029E5,
    1.60382976810944131506E4,
    1.96436237223387314144E2,
];

/* log(zeta(x) - 1 - 2**-x), 10 <= x <= 50 */
const A: [f64; 11] = [
    8.70728567484590192539E6,
    1.76506865670346462757E8,
    2.60889506707483264896E10,
    5.29806374009894791647E11,
    2.26888156119238241487E13,
    3.31884402932705083599E14,
    5.13778997975868230192E15,
    -1.98123688133907171455E15,
    -9.92763810039983572356E16,
    7.82905376180870586444E16,
    9.26786275768927717187E16,
];

const B: [f64; 10] = [
    /* 1.00000000000000000000E0, */
    -7.92625410563741062861E6,
    -1.60529969932920229676E8,
    -2.37669260975543221788E10,
    -4.80319584350455169857E11,
    -2.07820961754173320170E13,
    -2.96075404507272223680E14,
    -4.86299103694609136686E15,
    5.34589509675789930199E15,
    5.71464111092297631292E16,
    -1.79915597658676556828E16,
];

/* (1-x) (zeta(x) - 1), 0 <= x <= 1 */
const R: [f64; 6] = [
    -3.28717474506562731748E-1,
    1.55162528742623950834E1,
    -2.48762831680821954401E2,
    1.01050368053237678329E3,
    1.26726061410235149405E4,
    -1.11578094770515181334E5,
];

const S: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    1.95107674914060531512E1,
    3.17710311750646984099E2,
    3.03835500874445748734E3,
    2.03665876435770579345E4,
    7.43853965136767874343E4,
];

const TAYLOR0: [f64; 10] = [
    -1.0000000009110164892,
    -1.0000000057646759799,
    -9.9999983138417361078e-1,
    -1.0000013011460139596,
    -1.000001940896320456,
    -9.9987929950057116496e-1,
    -1.000785194477042408,
    -1.0031782279542924256,
    -9.1893853320467274178e-1,
    -1.5,
];

const MAXL2: f64 = 127.0;
const SQRT_2_PI: f64 = 0.79788456080286535587989;

use crate::cephes64::consts::{M_E, M_PI, MACHEP};
use crate::cephes64::polevl::{polevl, p1evl};
use crate::cephes64::lanczos::{lanczos_sum_expg_scaled, LANCZOS_G};
use crate::cephes64::zeta::zeta;

pub fn zetac(x: f64) -> f64 {
    //! Riemann zeta function, minus one
    //!
    //! ## DESCRIPTION:
    //!
    #![doc=include_str!("zetac.svg")]
    //!
    //! is related to the Riemann zeta function by
    //! Riemann zeta(x) = zetac(x) + 1.
    //!
    //! Extension of the function definition for x < 1 is implemented.
    //! Zero is returned for x > log2(INFINITY).
    //!
    //! ## ACCURACY:
    //!
    //! Tabulated values have full machine accuracy.
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
    //!     <td>1, 50</td>
    //!     <td>10000</td>
    //!     <td>9.8e-16</td>
    //!     <td>1.3e-16</td>
    //! </tr>
    //!</table>

    if x.is_nan() {
        x
    } else if x == -f64::INFINITY {
        f64::NAN
    } else if x < 0.0 && x > -0.01 {
        zetac_smallneg(x)
    } else if x < 0.0 {
        zeta_reflection(-x) - 1.0
    } else {
        zetac_positive(x)
    }
}

pub fn riemann_zeta(x: f64) -> f64 {
    //! Riemann zeta function
    if x.is_nan() {
        x
    } else if x == -f64::INFINITY {
        f64::NAN
    } else if x < 0.0 && x > -0.01 {
        1.0 + zetac_smallneg(x)
    } else if x < 0.0 {
        zeta_reflection(-x)
    } else {
        1.0 + zetac_positive(x)
    }
}


/*
* Compute zetac for positive arguments
*/
fn zetac_positive(x: f64) -> f64 {

    if x == 1.0 {
        return f64::INFINITY;
    } else if x >= MAXL2 {
        /* because first term is 2**-x */
        return 0.0;
    }

    /* Tabulated values for integer argument */
    let w = x.floor();
    if w == x {
        let i = x as usize;
        if i < 31 {
            return AZETAC[i];
        }
    }

    if x < 1.0 {
        let w = 1.0 - x;
        polevl(x, &R, 5) / (w * p1evl(x, &S, 5))
    } else if x <= 10.0 {
        let b = 2.0_f64.powf(x) * (x - 1.0);
        let w = 1.0 / x;
        (x * polevl(w, &P, 8)) / (b * p1evl(w, &Q, 8))
    } else if x <= 50.0 {
        let b = 2.0_f64.powf(-x);
        let w = polevl(x, &A, 10) / p1evl(x, &B, 10);
        w.exp() + b
    } else {
        /* Basic sum of inverse powers */
        let mut s: f64 = 0.0;
        let mut a: f64 = 1.0;
        loop {
            a += 2.0;
            let b = a.powf(-x);
            s += b;
            if b / s <= MACHEP {
                break
            }
        }

        let b = 2.0_f64.powf(-x);

        (s + b) / (1.0 - b)
    }
}


/*
* Compute zetac for small negative x. We can't use the reflection
* formula because to double precision 1 - x = 1 and zetac(1) = inf.
*/
fn zetac_smallneg(x: f64) -> f64 {
    polevl(x, &TAYLOR0, 9)
}


/*
* Compute zetac using the reflection formula (see DLMF 25.4.2) plus
* the Lanczos approximation for Gamma to avoid overflow.
*/
fn zeta_reflection(x: f64) -> f64 {
    //double base, large_term, small_term, hx, x_shift;

    let hx = x / 2.0;
    if hx == hx.floor() {
        /* Hit a zero of the sine factor */
        return 0.0;
    }

    /* Reduce the argument to sine */
    let x_shift = x % 4.0; 
    let mut small_term = -SQRT_2_PI * (0.5 * M_PI * x_shift).sin();
    small_term *= lanczos_sum_expg_scaled(x + 1.0) * zeta(x + 1.0, 1.0);

    /* Group large terms together to prevent overflow */
    let base = (x + LANCZOS_G + 0.5) / (2.0 * M_PI * M_E);
    let large_term = base.powf(x + 0.5);
    if !large_term.is_infinite() {
        return large_term * small_term;
    }
    /*
    * We overflowed, but we might be able to stave off overflow by
    * factoring in the small term earlier. To do this we compute
    *
    * (sqrt(large_term) * small_term) * sqrt(large_term)
    *
    * Since we only call this method for negative x bounded away from
    * zero, the small term can only be as small sine on that region;
    * i.e. about machine epsilon. This means that if the above still
    * overflows, then there was truly no avoiding it.
    */
    let large_term = base.powf(0.5 * x + 0.25);
    return (large_term * small_term) * large_term;
}

#[cfg(test)]
mod zetac_tests {
    use super::*;

    #[test]
    fn zetac_trivials() {
        assert_eq!(zetac(f64::NAN).is_nan(), true);
        assert_eq!(zetac(-f64::INFINITY).is_nan(), true);
        assert_eq!(zetac(f64::INFINITY), 0.0);
        assert_eq!(zetac(1.0), f64::INFINITY);
    }

    #[test]
    fn zetac_positive() {
        assert_eq!(zetac(0.0), -1.5);
        assert_eq!(zetac(2.0), 0.6449340668482264);
        assert_eq!(zetac(20.0), 9.539620338727962e-07);
        assert_eq!(zetac(30.0), 9.313274324196682e-10);
        assert_eq!(zetac(31.0), 4.656629065033784e-10);

        assert_eq!(zetac(127.1), 0.0);

        assert_eq!(zetac(1e-20), -1.5);
        assert_eq!(zetac(2e-10), -1.5000000001837879);
        assert_eq!(zetac(3e-5), -1.5000275690588838);
        assert_eq!(zetac(0.04), -1.5384293431032379);
        assert_eq!(zetac(0.5), -2.4603545088095866);
        assert_eq!(zetac(1.0 - 1e-10), -9999999173.019142);

        assert_eq!(zetac(1.0 + 1e-10), 9999999172.17357);
        assert_eq!(zetac(1.0 + 1e-5), 99999.57721573791);
        assert_eq!(zetac(2.1), 0.5602165335033622);
        assert_eq!(zetac(5.3), 0.029332205683219356);
        assert_eq!(zetac(7.5), 0.0058267275365228065);
        assert_eq!(zetac(10.0 - 1e-10), 0.0009945751278877884);

        assert_eq!(zetac(10.0 + 1e-10), 0.000994575127748382);
        assert_eq!(zetac(20.1), 8.900673596929275e-07);
        assert_eq!(zetac(30.5), 6.585473125700447e-10);
        assert_eq!(zetac(40.0), 9.094947840263888e-13);
        assert_eq!(zetac(50.0), 8.881784210930816e-16);

        assert_eq!(zetac(50.0 + 1e-10), 8.881784210315165e-16);
        assert_eq!(zetac(60.0), 8.673617380119933e-19);
        assert_eq!(zetac(100.0), 7.888609052210118e-31);
        assert_eq!(zetac(120.1), 7.01936006836698e-37);
    }

    #[test]
    fn zetac_small_neg() {
        assert_eq!(zetac(-1e-20), -1.5);
        assert_eq!(zetac(-2e-10), -1.4999999998162123);
        assert_eq!(zetac(-3e-5), -1.4999724327468373);
        assert_eq!(zetac(-4e-4), -1.49963258503121);
        assert_eq!(zetac(-5e-3), -1.4954302623133413);
        assert_eq!(zetac(-1e-2 + 1e-15), -1.490909941605338);
    }

    #[test]
    fn zetac_neg() {
        assert_eq!(zetac(-0.01), -1.4909099416053366);
        assert_eq!(zetac(-0.1), -1.4172280407673665);
        assert_eq!(zetac(-1.0), -1.0833333333333335);
        assert_eq!(zetac(-2.0), -1.0);
        assert_eq!(zetac(-15.0), -0.5567401960784306);
        assert_eq!(zetac(-1e2 + 2.2), -3.210779249560814e+74);
        assert_eq!(zetac(-2e2 + 1e-13), 1.0376641210003888e+202);
        assert_eq!(zetac(-1e5 + 2.2), -f64::INFINITY);
    }
}

#[cfg(test)]
mod reimann_zeta_tests {
    use super::*;

    #[test]
    fn riemann_zeta_trivials() {
        assert_eq!(riemann_zeta(f64::NAN).is_nan(), true);
        assert_eq!(riemann_zeta(-f64::INFINITY).is_nan(), true);
        assert_eq!(riemann_zeta(f64::INFINITY), 1.0);
        assert_eq!(riemann_zeta(1.0), f64::INFINITY);
    }

    #[test]
    fn riemann_zeta_positive() {
        assert_eq!(riemann_zeta(0.0), -0.5);
        assert_eq!(riemann_zeta(2.0), 1.6449340668482264);
        assert_eq!(riemann_zeta(20.0), 1.0000009539620338);
        assert_eq!(riemann_zeta(30.0), 1.0000000009313275);
        assert_eq!(riemann_zeta(31.0), 1.0000000004656628);

        assert_eq!(riemann_zeta(127.1), 1.0);

        assert_eq!(riemann_zeta(1e-20), -0.5);
        assert_eq!(riemann_zeta(2e-10), -0.5000000001837879);
        assert_eq!(riemann_zeta(3e-5), -0.5000275690588838);
        assert_eq!(riemann_zeta(0.04), -0.5384293431032379);
        assert_eq!(riemann_zeta(0.5), -1.4603545088095866);
        assert_eq!(riemann_zeta(1.0 - 1e-10), -9999999172.019142);

        assert_eq!(riemann_zeta(1.0 + 1e-10), 9999999173.17357);
        assert_eq!(riemann_zeta(1.0 + 1e-5), 100000.57721573791);
        assert_eq!(riemann_zeta(2.1), 1.5602165335033622);
        assert_eq!(riemann_zeta(5.3), 1.029332205683219356);
        assert_eq!(riemann_zeta(7.5), 1.0058267275365228065);
        assert_eq!(riemann_zeta(10.0 - 1e-10), 1.0009945751278877884);

        assert_eq!(riemann_zeta(10.0 + 1e-10), 1.000994575127748382);
        assert_eq!(riemann_zeta(20.1), 1.0000008900673596);
        assert_eq!(riemann_zeta(30.5), 1.0000000006585472);
        assert_eq!(riemann_zeta(40.0), 1.0000000000009095);
        assert_eq!(riemann_zeta(50.0), 1.0000000000000009);

        assert_eq!(riemann_zeta(50.0 + 1e-10), 1.0000000000000009);
        assert_eq!(riemann_zeta(60.0), 1.0);
        assert_eq!(riemann_zeta(100.0), 1.0);
        assert_eq!(riemann_zeta(120.1), 1.0);
    }

    #[test]
    fn riemann_zeta_small_neg() {
        assert_eq!(riemann_zeta(-1e-20), -0.5);
        assert_eq!(riemann_zeta(-2e-10), -0.49999999981621235);
        assert_eq!(riemann_zeta(-3e-5), -0.49997243274683734);
        assert_eq!(riemann_zeta(-4e-4), -0.49963258503121);
        assert_eq!(riemann_zeta(-5e-3), -0.4954302623133413);
        assert_eq!(riemann_zeta(-1e-2 + 1e-15), -0.49090994160533796);
    }

    #[test]
    fn riemann_zeta_neg() {
        assert_eq!(riemann_zeta(-0.01), -0.4909099416053367);
        assert_eq!(riemann_zeta(-0.1), -0.4172280407673665);
        assert_eq!(riemann_zeta(-1.0), -0.08333333333333338);
        assert_eq!(riemann_zeta(-2.0), 0.0);
        assert_eq!(riemann_zeta(-15.0), 0.4432598039215694);
        assert_eq!(riemann_zeta(-1e2 + 2.2), -3.210779249560814e+74);
        assert_eq!(riemann_zeta(-2e2 + 1e-13), 1.0376641210003888e+202);
        assert_eq!(riemann_zeta(-1e5 + 2.2), -f64::INFINITY);
    }
}