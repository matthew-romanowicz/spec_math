/*                                                     dawsn.c
*
*     Dawson's Integral
*
*
*
* SYNOPSIS:
*
* double x, y, dawsn();
*
* y = dawsn( x );
*
*
*
* DESCRIPTION:
*
* Approximates the integral
*
*                             x
*                             -
*                      2     | |        2
*  dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
*                          | |
*                           -
*                           0
*
* Three different rational approximations are employed, for
* the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0,10        10000       6.9e-16     1.0e-16
*
*
*/

/*                                                     dawsn.c */


/*
* Cephes Math Library Release 2.1:  January, 1989
* Copyright 1984, 1987, 1989 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

use crate::cephes64::polevl::{polevl, p1evl};

/* Dawson's integral, interval 0 to 3.25 */

const AN: [f64; 10] = [
    1.13681498971755972054E-11,
    8.49262267667473811108E-10,
    1.94434204175553054283E-8,
    9.53151741254484363489E-7,
    3.07828309874913200438E-6,
    3.52513368520288738649E-4,
    -8.50149846724410912031E-4,
    4.22618223005546594270E-2,
    -9.17480371773452345351E-2,
    9.99999999999999994612E-1,
];

const AD: [f64; 11] = [
    2.40372073066762605484E-11,
    1.48864681368493396752E-9,
    5.21265281010541664570E-8,
    1.27258478273186970203E-6,
    2.32490249820789513991E-5,
    3.25524741826057911661E-4,
    3.48805814657162590916E-3,
    2.79448531198828973716E-2,
    1.58874241960120565368E-1,
    5.74918629489320327824E-1,
    1.00000000000000000539E0,
];

/* interval 3.25 to 6.25 */
const BN: [f64; 11] = [
    5.08955156417900903354E-1,
    -2.44754418142697847934E-1,
    9.41512335303534411857E-2,
    -2.18711255142039025206E-2,
    3.66207612329569181322E-3,
    -4.23209114460388756528E-4,
    3.59641304793896631888E-5,
    -2.14640351719968974225E-6,
    9.10010780076391431042E-8,
    -2.40274520828250956942E-9,
    3.59233385440928410398E-11,
];

const BD: [f64; 10] = [
    /*  1.00000000000000000000E0, */
    -6.31839869873368190192E-1,
    2.36706788228248691528E-1,
    -5.31806367003223277662E-2,
    8.48041718586295374409E-3,
    -9.47996768486665330168E-4,
    7.81025592944552338085E-5,
    -4.55875153252442634831E-6,
    1.89100358111421846170E-7,
    -4.91324691331920606875E-9,
    7.18466403235734541950E-11,
];

/* 6.25 to infinity */
const CN: [f64; 5] = [
    -5.90592860534773254987E-1,
    6.29235242724368800674E-1,
    -1.72858975380388136411E-1,
    1.64837047825189632310E-2,
    -4.86827613020462700845E-4,
];

const CD: [f64; 5] = [
    /* 1.00000000000000000000E0, */
    -2.69820057197544900361E0,
    1.73270799045947845857E0,
    -3.93708582281939493482E-1,
    3.44278924041233391079E-2,
    -9.73655226040941223894E-4,
];

// $$\mathrm{dawsn}(x) = \exp(-x^2)\int_{0}^{x}{\exp(t^2)\,dt}$$

pub fn dawsn(xx: f64) -> f64
{
    //! Dawson's Integral
    //!
    //! ## DESCRIPTION:
    //!
    //! Approximates the integral
    //!
    #![doc=include_str!("dawsn.svg")]
    //!
    //! Three different rational approximations are employed, for
    //! the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.
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
    //!     <td>0, 10</td>
    //!     <td>10000</td>
    //!     <td>6.9e-16</td>
    //!     <td>1.0e-16</td>
    //! </tr>
    //!</table>


    let mut sign: isize = 1;
    let mut xx = xx;
    if xx < 0.0 {
        sign = -1;
        xx = -xx;
    }

    if xx < 3.25 {
        let x = xx * xx;
        let y = xx * polevl(x, &AN, 9) / polevl(x, &AD, 10);
        sign as f64 * y
    } else if xx < 6.25 {
        let x = 1.0 / (xx * xx);
        let y = 1.0 / xx + x * polevl(x, &BN, 10) / (p1evl(x, &BD, 10) * xx);
        sign as f64 * 0.5 * y
    } else if xx > 1.0e9 {
        (sign as f64 * 0.5) / xx
    } else {
        /* 6.25 to infinity */
        let x = 1.0 / (xx * xx);
        let y = 1.0 / xx + x * polevl(x, &CN, 4) / (p1evl(x, &CD, 5) * xx);
        sign as f64 * 0.5 * y
    }
}

#[cfg(test)]
mod dawsn_tests {
    use super::*;

    #[test]
    fn dawsn_trivials() {
        assert_eq!(dawsn(f64::NAN).is_nan(), true);
        assert_eq!(dawsn(f64::INFINITY), 0.0);
        assert_eq!(dawsn(-f64::INFINITY), 0.0);
    }

    #[test]
    fn dawsn_small() {
        assert_eq!(dawsn(0.0), 0.0);
        assert_eq!(dawsn(1.0), 0.5380795069127684);
        assert_eq!(dawsn(2.0), 0.30134038892379195);
        assert_eq!(dawsn(3.0), 0.17827103061055835);
        assert_eq!(dawsn(3.25 - 1e-16), 0.162570914560687);

        assert_eq!(dawsn(-1.0), -0.5380795069127684);
        assert_eq!(dawsn(-2.0), -0.30134038892379195);
        assert_eq!(dawsn(-3.0), -0.17827103061055835);
        assert_eq!(dawsn(-3.25 + 1e-16), -0.162570914560687);
    }

    #[test]
    fn dawsn_medium() {
        assert_eq!(dawsn(3.25), 0.162570914560687);
        assert_eq!(dawsn(4.0), 0.1293480012360051);
        assert_eq!(dawsn(5.0), 0.10213407442427684);
        assert_eq!(dawsn(6.0), 0.08454268897454385);
        assert_eq!(dawsn(6.25 - 1e-16), 0.08106609406101173);

        assert_eq!(dawsn(-3.25), -0.162570914560687);
        assert_eq!(dawsn(-4.0), -0.1293480012360051);
        assert_eq!(dawsn(-5.0), -0.10213407442427684);
        assert_eq!(dawsn(-6.0), -0.08454268897454385);
        assert_eq!(dawsn(-6.25 + 1e-16), -0.08106609406101173);
    }

    #[test]
    fn dawsn_nominal() {
        assert_eq!(dawsn(6.25), 0.08106609406101173);
        assert_eq!(dawsn(10.0), 0.05025384718759853);
        assert_eq!(dawsn(100.0), 0.005000250037509379);
        assert_eq!(dawsn(1e5), 5.0000000002500005e-06);
        assert_eq!(dawsn(1e9 - 1e-3), 5.000000000005001e-10);

        assert_eq!(dawsn(-6.25), -0.08106609406101173);
        assert_eq!(dawsn(-10.0), -0.05025384718759853);
        assert_eq!(dawsn(-100.0), -0.005000250037509379);
        assert_eq!(dawsn(-1e5), -5.0000000002500005e-06);
        assert_eq!(dawsn(-1e9 + 1e-3), -5.000000000005001e-10);
    }

    #[test]
    fn dawsn_large() {
        assert_eq!(dawsn(1e9), 5e-10);
        assert_eq!(dawsn(2e15), 2.5e-16);
        assert_eq!(dawsn(3e20), 1.6666666666666666e-21);
        assert_eq!(dawsn(4e50), 1.2499999999999999e-51);
        assert_eq!(dawsn(5e250), 1e-251);

        assert_eq!(dawsn(-1e9), -5e-10);
        assert_eq!(dawsn(-2e15), -2.5e-16);
        assert_eq!(dawsn(-3e20), -1.6666666666666666e-21);
        assert_eq!(dawsn(-4e50), -1.2499999999999999e-51);
        assert_eq!(dawsn(-5e250), -1e-251);
    }
}