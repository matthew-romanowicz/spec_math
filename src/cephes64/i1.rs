/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1985, 1987, 2000 by Stephen L. Moshier
*/

#![allow(clippy::excessive_precision)]

use crate::cephes64::chbevl;

/* Chebyshev coefficients for exp(-x) I1(x) / x
* in the interval [0,8].
*
* lim(x->0){ exp(-x) I1(x) / x } = 1/2.
*/

const A: [f64; 29] = [
    2.77791411276104639959E-18,
    -2.11142121435816608115E-17,
    1.55363195773620046921E-16,
    -1.10559694773538630805E-15,
    7.60068429473540693410E-15,
    -5.04218550472791168711E-14,
    3.22379336594557470981E-13,
    -1.98397439776494371520E-12,
    1.17361862988909016308E-11,
    -6.66348972350202774223E-11,
    3.62559028155211703701E-10,
    -1.88724975172282928790E-9,
    9.38153738649577178388E-9,
    -4.44505912879632808065E-8,
    2.00329475355213526229E-7,
    -8.56872026469545474066E-7,
    3.47025130813767847674E-6,
    -1.32731636560394358279E-5,
    4.78156510755005422638E-5,
    -1.61760815825896745588E-4,
    5.12285956168575772895E-4,
    -1.51357245063125314899E-3,
    4.15642294431288815669E-3,
    -1.05640848946261981558E-2,
    2.47264490306265168283E-2,
    -5.29459812080949914269E-2,
    1.02643658689847095384E-1,
    -1.76416518357834055153E-1,
    2.52587186443633654823E-1
];

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
* in the inverted interval [8,infinity].
*
* lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
*/
const B: [f64; 25] = [
    7.51729631084210481353E-18,
    4.41434832307170791151E-18,
    -4.65030536848935832153E-17,
    -3.20952592199342395980E-17,
    2.96262899764595013876E-16,
    3.30820231092092828324E-16,
    -1.88035477551078244854E-15,
    -3.81440307243700780478E-15,
    1.04202769841288027642E-14,
    4.27244001671195135429E-14,
    -2.10154184277266431302E-14,
    -4.08355111109219731823E-13,
    -7.19855177624590851209E-13,
    2.03562854414708950722E-12,
    1.41258074366137813316E-11,
    3.25260358301548823856E-11,
    -1.89749581235054123450E-11,
    -5.58974346219658380687E-10,
    -3.83538038596423702205E-9,
    -2.63146884688951950684E-8,
    -2.51223623787020892529E-7,
    -3.88256480887769039346E-6,
    -1.10588938762623716291E-4,
    -9.76109749136146840777E-3,
    7.78576235018280120474E-1
];

pub fn i1(x: f64) -> f64 {
    //! Modified Bessel function of order one
    //! 
    //! ## DESCRIPTION:
    //!
    //! Returns modified Bessel function of order one of the
    //! argument.
    //!
    //! The function is defined as i1(x) = -i j1( ix ).
    //!
    //! The range is partitioned into the two intervals [0,8] and
    //! (8, infinity).  Chebyshev polynomial expansions are employed
    //! in each interval.
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
    //!     <td>1.9e-15</td>
    //!     <td>2.1e-16</td>
    //! </tr>
    //!</table>
    let z = x.abs();
    let z = if z <= 8.0 {
        chbevl((z / 2.0) - 2.0, &A, 29) * z * z.exp()
    } else {
        z.exp() * chbevl(32.0 / z - 2.0, &B, 25) / z.sqrt()
    };
    if x < 0.0 {
        -z
    } else {
        z
    }
}

pub fn i1e(x: f64) -> f64 {
    //! Modified Bessel function of order one, exponentially scaled
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns exponentially scaled modified Bessel function
    //! of order one of the argument.
    //!
    //! The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
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
    //!     <td>2.0e-15</td>
    //!     <td>2.0e-16</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::i1`](crate::cephes64::i1)

    let z = x.abs();
    let z = if z <= 8.0 {
        chbevl((z / 2.0) - 2.0, &A, 29) * z
    } else {
        chbevl(32.0 / z - 2.0, &B, 25) / z.sqrt()
    };
    if x < 0.0 {
        -z
    } else {
        z
    }
}

#[cfg(test)]
mod i1_tests {
    use super::*;

    #[test]
    fn i1_trivials() {
        assert_eq!(i1(f64::INFINITY).is_nan(), true);
        assert_eq!(i1(-f64::INFINITY).is_nan(), true);
        assert_eq!(i1(f64::NAN).is_nan(), true);
    }

    #[test]
    fn i1_small_x() {
        assert_eq!(i1(0.0), 0.0);

        assert_eq!(i1(0.1), 0.050062526047092694);
        assert_eq!(i1(1.0), 0.5651591039924851);
        assert_eq!(i1(4.0), 9.759465153704449);
        assert_eq!(i1(7.0), 156.03909286995534);
        assert_eq!(i1(8.0), 399.8731367825599);

        assert_eq!(i1(-0.1), -0.050062526047092694);
        assert_eq!(i1(-1.0), -0.5651591039924851);
        assert_eq!(i1(-4.0), -9.759465153704449);
        assert_eq!(i1(-7.0), -156.03909286995534);
        assert_eq!(i1(-8.0), -399.8731367825599);
    }

    #[test]
    fn i1_large_x() {
        assert_eq!(i1(8.0 + 1e-16), 399.8731367825599);
        assert_eq!(i1(10.0), 2670.988303701255);
        assert_eq!(i1(20.0), 42454973.385127775);
        assert_eq!(i1(100.0), 1.0683693903381625e+42);

        assert_eq!(i1(-8.0 - 1e-16), -399.8731367825599);
        assert_eq!(i1(-10.0), -2670.988303701255);
        assert_eq!(i1(-20.0), -42454973.385127775);
        assert_eq!(i1(-100.0), -1.0683693903381625e+42);
    }
}

#[cfg(test)]
mod i1e_tests {
    use super::*;

    #[test]
    fn i1e_trivials() {
        assert_eq!(i1e(f64::INFINITY), 0.0);
        assert_eq!(i1e(-f64::INFINITY), 0.0);
        assert_eq!(i1e(f64::NAN).is_nan(), true);
    }

    #[test]
    fn i0e_small_x() {
        assert_eq!(i1e(0.0), 0.0);

        assert_eq!(i1e(0.1), 0.045298446808809324);
        assert_eq!(i1e(1.0), 0.2079104153497085);
        assert_eq!(i1e(4.0), 0.1787508395024353);
        assert_eq!(i1e(7.0), 0.1422892347095986);
        assert_eq!(i1e(8.0), 0.13414249329269812);

        assert_eq!(i1e(-0.1), -0.045298446808809324);
        assert_eq!(i1e(-1.0), -0.2079104153497085);
        assert_eq!(i1e(-4.0), -0.1787508395024353);
        assert_eq!(i1e(-7.0), -0.1422892347095986);
        assert_eq!(i1e(-8.0), -0.13414249329269812);

    }

    #[test]
    fn i0e_large_x() {
        assert_eq!(i1e(8.0 + 1e-16), 0.13414249329269812);
        assert_eq!(i1e(10.0), 0.1212626813844555);
        assert_eq!(i1e(20.0), 0.08750622218328867);
        assert_eq!(i1e(100.0), 0.03974415302513025);
        assert_eq!(i1e(1000.0), 0.01261093025692863);
        assert_eq!(i1e(1e10), 3.989422803864723e-06);
        assert_eq!(i1e(1e20), 3.989422804014327e-11);

        assert_eq!(i1e(-8.0 - 1e-16), -0.13414249329269812);
        assert_eq!(i1e(-10.0), -0.1212626813844555);
        assert_eq!(i1e(-20.0), -0.08750622218328867);
        assert_eq!(i1e(-100.0), -0.03974415302513025);
        assert_eq!(i1e(-1000.0), -0.01261093025692863);
        assert_eq!(i1e(-1e10), -3.989422803864723e-06);
        assert_eq!(i1e(-1e20), -3.989422804014327e-11);
    }
}