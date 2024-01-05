/*                                                     shichi.c
*
*     Hyperbolic sine and cosine integrals
*
*
*
* SYNOPSIS:
*
* double x, Chi, Shi, shichi();
*
* shichi( x, &Chi, &Shi );
*
*
* DESCRIPTION:
*
* Approximates the integrals
*
*                            x
*                            -
*                           | |   cosh t - 1
*   Chi(x) = eul + ln x +   |    -----------  dt,
*                         | |          t
*                          -
*                          0
*
*               x
*               -
*              | |  sinh t
*   Shi(x) =   |    ------  dt
*            | |       t
*             -
*             0
*
* where eul = 0.57721566490153286061 is Euler's constant.
* The integrals are evaluated by power series for x < 8
* and by Chebyshev expansions for x between 8 and 88.
* For large x, both functions approach exp(x)/2x.
* Arguments greater than 88 in magnitude return INFINITY.
*
*
* ACCURACY:
*
* Test interval 0 to 88.
*                      Relative error:
* arithmetic   function  # trials      peak         rms
*    IEEE         Shi      30000       6.9e-16     1.6e-16
*        Absolute error, except relative when |Chi| > 1:
*    IEEE         Chi      30000       8.4e-16     1.4e-16
*/

/*
* Cephes Math Library Release 2.0:  April, 1987
* Copyright 1984, 1987 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#![allow(clippy::excessive_precision)]


/* x exp(-x) shi(x), inverted interval 8 to 18 */
const S1: [f64; 22] = [
    1.83889230173399459482E-17,
    -9.55485532279655569575E-17,
    2.04326105980879882648E-16,
    1.09896949074905343022E-15,
    -1.31313534344092599234E-14,
    5.93976226264314278932E-14,
    -3.47197010497749154755E-14,
    -1.40059764613117131000E-12,
    9.49044626224223543299E-12,
    -1.61596181145435454033E-11,
    -1.77899784436430310321E-10,
    1.35455469767246947469E-9,
    -1.03257121792819495123E-9,
    -3.56699611114982536845E-8,
    1.44818877384267342057E-7,
    7.82018215184051295296E-7,
    -5.39919118403805073710E-6,
    -3.12458202168959833422E-5,
    8.90136741950727517826E-5,
    2.02558474743846862168E-3,
    2.96064440855633256972E-2,
    1.11847751047257036625E0
];

/* x exp(-x) shi(x), inverted interval 18 to 88 */
const S2: [f64; 23] = [
    -1.05311574154850938805E-17,
    2.62446095596355225821E-17,
    8.82090135625368160657E-17,
    -3.38459811878103047136E-16,
    -8.30608026366935789136E-16,
    3.93397875437050071776E-15,
    1.01765565969729044505E-14,
    -4.21128170307640802703E-14,
    -1.60818204519802480035E-13,
    3.34714954175994481761E-13,
    2.72600352129153073807E-12,
    1.66894954752839083608E-12,
    -3.49278141024730899554E-11,
    -1.58580661666482709598E-10,
    -1.79289437183355633342E-10,
    1.76281629144264523277E-9,
    1.69050228879421288846E-8,
    1.25391771228487041649E-7,
    1.16229947068677338732E-6,
    1.61038260117376323993E-5,
    3.49810375601053973070E-4,
    1.28478065259647610779E-2,
    1.03665722588798326712E0
];

/* x exp(-x) chin(x), inverted interval 8 to 18 */
const C1: [f64; 23] = [
    -8.12435385225864036372E-18,
    2.17586413290339214377E-17,
    5.22624394924072204667E-17,
    -9.48812110591690559363E-16,
    5.35546311647465209166E-15,
    -1.21009970113732918701E-14,
    -6.00865178553447437951E-14,
    7.16339649156028587775E-13,
    -2.93496072607599856104E-12,
    -1.40359438136491256904E-12,
    8.76302288609054966081E-11,
    -4.40092476213282340617E-10,
    -1.87992075640569295479E-10,
    1.31458150989474594064E-8,
    -4.75513930924765465590E-8,
    -2.21775018801848880741E-7,
    1.94635531373272490962E-6,
    4.33505889257316408893E-6,
    -6.13387001076494349496E-5,
    -3.13085477492997465138E-4,
    4.97164789823116062801E-4,
    2.64347496031374526641E-2,
    1.11446150876699213025E0
];

/* x exp(-x) chin(x), inverted interval 18 to 88 */
const C2: [f64; 24] = [
    8.06913408255155572081E-18,
    -2.08074168180148170312E-17,
    -5.98111329658272336816E-17,
    2.68533951085945765591E-16,
    4.52313941698904694774E-16,
    -3.10734917335299464535E-15,
    -4.42823207332531972288E-15,
    3.49639695410806959872E-14,
    6.63406731718911586609E-14,
    -3.71902448093119218395E-13,
    -1.27135418132338309016E-12,
    2.74851141935315395333E-12,
    2.33781843985453438400E-11,
    2.71436006377612442764E-11,
    -2.56600180000355990529E-10,
    -1.61021375163803438552E-9,
    -4.72543064876271773512E-9,
    -3.00095178028681682282E-9,
    7.79387474390914922337E-8,
    1.06942765566401507066E-6,
    1.59503164802313196374E-5,
    3.49592575153777996871E-4,
    1.28475387530065247392E-2,
    1.03665693917934275131E0
];

/* Sine and cosine integrals */

use crate::cephes64::consts::{MACHEP, EULER};
use crate::cephes64::chbevl::chbevl;

// $$\mathrm{Chi}(x) = \mathrm{EUL} + \ln(x) + \int_0^x{\frac{\cosh(t)-1}{t}\,dt}$$
// $$\mathrm{Shi}(x) = \int_0^x{\frac{\sinh(t)}{t}\,dt}$$

pub fn shichi(x: f64) -> (f64, f64) { // si, ci
    //! ## Hyperbolic sine and cosine integrals
    //!
    //! ## DESCRIPTION:
    //!
    //! Approximates the integrals
    //!
    #![doc=include_str!("chi.svg")]
    //!
    #![doc=include_str!("shi.svg")]
    //!
    //! where eul = 0.57721566490153286061 is Euler's constant.
    //! The integrals are evaluated by power series for x < 8
    //! and by Chebyshev expansions for x between 8 and 88.
    //! For large x, both functions approach exp(x)/2x.
    //! Arguments greater than 88 in magnitude return INFINITY.
    //!
    //! ## ACCURACY:
    //!
    //! Test interval 0 to 88.
    //!
    //! Relative error:
    //!
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
    //!     <td>Shi(x)</td>
    //!     <td>30000</td>
    //!     <td>6.9e-16</td>
    //!     <td>1.6e-16</td>
    //! </tr>
    //!</table>
    //!
    //! Absolute error, except relative when |Chi| > 1:
    //!
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
    //!     <td>Chi(x)</td>
    //!     <td>30000</td>
    //!     <td>8.4e-16</td>
    //!     <td>1.4e-16</td>
    //! </tr>
    //!</table>

    let mut x = x;

    let sign = if x < 0.0 {
        x = -x;
        -1
    } else {
        0
    };


    if x == 0.0 {
        return (0.0, -f64::INFINITY)
    }

    let mut s: f64;
    let mut c: f64;

    if x >= 8.0 {
        //goto chb;
        /* Chebyshev series expansions */
        if x < 18.0 {
            let a = (576.0 / x - 52.0) / 10.0;
            let k = x.exp() / x;
            s = k * chbevl(a, &S1, 22);
            c = k * chbevl(a, &C1, 23);
        } else if x <= 88.0 {
            let a = (6336.0 / x - 212.0) / 70.0;
            let k = x.exp() / x;
            s = k * chbevl(a, &S2, 23);
            c = k * chbevl(a, &C2, 24);
        } else if x > 1000.0 {
                if sign != 0 {
                    return (-f64::INFINITY, f64::INFINITY);
                } else {
                    return (f64::INFINITY, f64::INFINITY);
                }
        } else {
            /* Asymptotic expansions
            * http://functions.wolfram.com/GammaBetaErf/CoshIntegral/06/02/
            * http://functions.wolfram.com/GammaBetaErf/SinhIntegral/06/02/0001/
            */
            let a = hyp3f0(0.5, 1.0, 1.0, 4.0/(x*x));
            let b = hyp3f0(1.0, 1.0, 1.5, 4.0/(x*x));
            let si = x.cosh()/x * a + x.sinh()/(x*x) * b;
            let ci = x.sinh()/x * a + x.cosh()/(x*x) * b;
            if sign != 0 {
                return (-si, ci)
            } else {
                return (si, ci)
            }
        }
    } else {

        let z = x * x;

        /*     Direct power series expansion   */
        let mut a = 1.0;
        let mut k = 2.0;
        s = 1.0;
        c = 0.0;

        loop {
            a *= z / k;
            c += a / k;
            k += 1.0;
            a /= k;
            s += a / k;
            k += 1.0;
            if !((a / s).abs() > MACHEP) {
                break;
            }
        }

        s *= x;
    }



    //done:
    if sign != 0 {
        s = -s
    }

    (s, EULER + x.ln() + c)
}


/*
* Evaluate 3F0(a1, a2, a3; z)
*
* The series is only asymptotic, so this requires z large enough.
*/
fn hyp3f0(a1: f64, a2: f64, a3: f64, z: f64) -> f64 {

    let m = z.powf(-1.0/3.0);
    let maxiter = std::cmp::min(m as isize, 50);

    let mut term: f64 = 1.0;
    let mut sum = term;
    for n in 0..maxiter {
        term *= (a1 + n as f64) * (a2 + n as f64) * (a3 + n as f64) * z / (n + 1) as f64;
        sum += term;
        if term.abs() < 1e-13 * sum.abs() || term == 0.0 {
            break;
        }
    }

    let err = term.abs();

    if err > 1e-13 * sum.abs() {
        f64::NAN
    } else {
        sum
    }
}

#[cfg(test)]
mod shichi_tests {
    use super::*;

    #[test]
    fn shichi_trivials() {
        assert_eq!(shichi(f64::NAN).0.is_nan(), true);
        assert_eq!(shichi(f64::NAN).1.is_nan(), true);
        assert_eq!(shichi(0.0), (0.0, -f64::INFINITY));
        assert_eq!(shichi(f64::INFINITY), (f64::INFINITY, f64::INFINITY));
        assert_eq!(shichi(-f64::INFINITY), (-f64::INFINITY, f64::INFINITY));
    }

    #[test]
    fn shichi_medium() { // 8.0 <= x < 18.0
        assert_eq!(shichi(8.0), (220.18996860023054, 220.18993093460776));
        assert_eq!(shichi(10.0), (1246.1144901994235, 1246.1144860424545));
        assert_eq!(shichi(14.0), (46596.256817010515, 46596.25681695486));
        assert_eq!(shichi(18.0 - 1e-10), (1938952.1651163367, 1938952.1651163357));
    
        assert_eq!(shichi(-8.0), (-220.18996860023054, 220.18993093460776));
        assert_eq!(shichi(-10.0), (-1246.1144901994235, 1246.1144860424545));
        assert_eq!(shichi(-14.0), (-46596.256817010515, 46596.25681695486));
        assert_eq!(shichi(-18.0 + 1e-10), (-1938952.1651163367, 1938952.1651163357));
    }

    #[test]
    fn shichi_large() { // 18.0 <= x <= 88.0
        assert_eq!(shichi(18.0), (1938952.165298722, 1938952.1652987215));
        assert_eq!(shichi(30.0), (184486604703.63712, 184486604703.63712));
        assert_eq!(shichi(50.0), (5.292818448565846e+19, 5.292818448565845e+19));
        assert_eq!(shichi(88.0), (9.493446879912849e+35, 9.493446879912846e+35));

        assert_eq!(shichi(-18.0), (-1938952.165298722, 1938952.1652987215));
        assert_eq!(shichi(-30.0), (-184486604703.63712, 184486604703.63712));
        assert_eq!(shichi(-50.0), (-5.292818448565846e+19, 5.292818448565845e+19));
        assert_eq!(shichi(-88.0), (-9.493446879912849e+35, 9.493446879912846e+35));

    }

    #[test]
    fn shichi_small() { // x < 8.0
        assert_eq!(shichi(1e-300), (1e-300, -690.1983122333121));
        assert_eq!(shichi(1e-150), (1e-150, -344.8105482842053));
        assert_eq!(shichi(1e-20), (1e-20, -45.474486194979384));
        assert_eq!(shichi(1e-10), (1e-10, -22.448635265038924));
        assert_eq!(shichi(1e-5), (1.0000000000055556e-05, -10.935709800043695));
        assert_eq!(shichi(0.1), (0.10005557222505701, -1.7228683861943335));
        assert_eq!(shichi(1.0), (1.0572508753757286, 0.8378669409802082));
        assert_eq!(shichi(2.0), (2.501567433354976, 2.4526669226469147));
        assert_eq!(shichi(4.0), (9.817326911233032, 9.813547558823187));
        assert_eq!(shichi(6.0), (42.99506111244568, 42.99470102999352));
        assert_eq!(shichi(8.0 - 1e-20), (220.18996860023054, 220.18993093460776));

        assert_eq!(shichi(-1e-300), (-1e-300, -690.1983122333121));
        assert_eq!(shichi(-1e-150), (-1e-150, -344.8105482842053));
        assert_eq!(shichi(-1e-20), (-1e-20, -45.474486194979384));
        assert_eq!(shichi(-1e-10), (-1e-10, -22.448635265038924));
        assert_eq!(shichi(-1e-5), (-1.0000000000055556e-05, -10.935709800043695));
        assert_eq!(shichi(-0.1), (-0.10005557222505701, -1.7228683861943335));
        assert_eq!(shichi(-1.0), (-1.0572508753757286, 0.8378669409802082));
        assert_eq!(shichi(-2.0), (-2.501567433354976, 2.4526669226469147));
        assert_eq!(shichi(-4.0), (-9.817326911233032, 9.813547558823187));
        assert_eq!(shichi(-6.0), (-42.99506111244568, 42.99470102999352));
        assert_eq!(shichi(-8.0 + 1e-20), (-220.18996860023054, 220.18993093460776));
    }

    #[test]
    fn shichi_asy() { // 88.0 < x <= 1000.0
        assert_eq!(shichi(88.0 + 1e-10), (9.493446880851297e+35, 9.493446880851297e+35));
        assert_eq!(shichi(90.0), (6.857084347536261e+36, 6.857084347536261e+36));
        assert_eq!(shichi(100.0), (1.3577763724269395e+41, 1.3577763724269395e+41));
        assert_eq!(shichi(200.0), (1.8156176165796781e+84, 1.8156176165796781e+84));
        assert_eq!(shichi(400.0), (6.54323640853714e+170, 6.54323640853714e+170));
        assert_eq!(shichi(700.0), (7.254893680262804e+300, 7.254893680262804e+300));
        assert_eq!(shichi(800.0), (f64::INFINITY, f64::INFINITY));
        assert_eq!(shichi(1000.0), (f64::INFINITY, f64::INFINITY));

        assert_eq!(shichi(-88.0 - 1e-10), (-9.493446880851297e+35, 9.493446880851297e+35));
        assert_eq!(shichi(-90.0), (-6.857084347536261e+36, 6.857084347536261e+36));
        assert_eq!(shichi(-100.0), (-1.3577763724269395e+41, 1.3577763724269395e+41));
        assert_eq!(shichi(-200.0), (-1.8156176165796781e+84, 1.8156176165796781e+84));
        assert_eq!(shichi(-400.0), (-6.54323640853714e+170, 6.54323640853714e+170));
        assert_eq!(shichi(-700.0), (-7.254893680262804e+300, 7.254893680262804e+300));
        assert_eq!(shichi(-800.0), (-f64::INFINITY, f64::INFINITY));
        assert_eq!(shichi(-1000.0), (-f64::INFINITY, f64::INFINITY));
    }
}