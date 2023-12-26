/*
* Cephes Math Library Release 2.8:  June, 2000
* Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#![allow(clippy::excessive_precision)]

use crate::cephes64::chbevl;

/* Chebyshev coefficients for exp(-x) I0(x)
* in the interval [0,8].
*
* lim(x->0){ exp(-x) I0(x) } = 1.
*/
const A: [f64; 30] = [
    -4.41534164647933937950E-18,
    3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
    1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
    7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
    2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
    9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
    2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
    6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
    1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
    1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
    1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
    1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
    1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
    4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
    1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
    6.76795274409476084995E-1
];

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
* in the inverted interval [8,infinity].
*
* lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
*/
const B: [f64; 25] = [
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
    4.46562142029675999901E-17,
    3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
    1.77256013305652638360E-15,
    3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
    1.54008621752140982691E-14,
    3.85277838274214270114E-13,
    7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
    1.18891471078464383424E-11,
    4.94060238822496958910E-10,
    3.39623202570838634515E-9,
    2.26666899049817806459E-8,
    2.04891858946906374183E-7,
    2.89137052083475648297E-6,
    6.88975834691682398426E-5,
    3.36911647825569408990E-3,
    8.04490411014108831608E-1
];

pub fn i0(x: f64) -> f64 {
    //! Modified Bessel function of order zero
    //!
    //! ## DESCRIPTION:
    //! 
    //! Returns modified Bessel function of order zero of the
    //! argument.
    //!
    //! The function is defined as i0(x) = j0( ix ).
    //!
    //! The range is partitioned into the two intervals `[0,8]` and
    //! `(8, infinity)`.  Chebyshev polynomial expansions are employed
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
    //!     <td>5.8e-16</td>
    //!     <td>1.4e-16</td>
    //! </tr>
    //!</table>

    let x = x.abs();
    if x <= 8.0 {
        x.exp() * chbevl(x * 0.5 - 2.0, &A, 30)
    } else {
        x.exp() * chbevl(32.0 / x - 2.0, &B, 25) / x.sqrt()
    }

}




pub fn i0e(x: f64) -> f64 {
    //!Modified Bessel function of order zero, exponentially scaled
    //!
    //! ## DESCRIPTION:
    //!
    //! Returns exponentially scaled modified Bessel function
    //! of order zero of the argument.
    //!
    //! The function is defined as i0e(x) = exp(-|x|) j0( ix ).
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
    //!     <td>5.4e-16</td>
    //!     <td>1.2e-16</td>
    //! </tr>
    //!</table>
    //!
    //! See [`cephes64::i0`](crate::cephes64::i0)

    let x = x.abs();
    if x <= 8.0 {
       chbevl(x * 0.5 - 2.0, &A, 30)
    } else {
        chbevl(32.0 / x - 2.0, &B, 25) / x.sqrt()
    }

}

#[cfg(test)]
mod i0_tests {
    use super::*;

    #[test]
    fn i0_trivials() {
        assert_eq!(i0(f64::INFINITY).is_nan(), true);
        assert_eq!(i0(-f64::INFINITY).is_nan(), true);
        assert_eq!(i0(f64::NAN).is_nan(), true);
    }

    #[test]
    fn i0_small_x() {
        assert_eq!(i0(0.0), 1.0);

        assert_eq!(i0(0.1), 1.0025015629340956);
        assert_eq!(i0(1.0), 1.2660658777520082);
        assert_eq!(i0(4.0), 11.30192195213633);
        assert_eq!(i0(7.0), 168.59390851028968);
        assert_eq!(i0(8.0), 427.56411572180474);

        assert_eq!(i0(-0.1), 1.0025015629340956);
        assert_eq!(i0(-1.0), 1.2660658777520082);
        assert_eq!(i0(-4.0), 11.30192195213633);
        assert_eq!(i0(-7.0), 168.59390851028968);
        assert_eq!(i0(-8.0), 427.56411572180474);
    }

    #[test]
    fn i0_large_x() {
        assert_eq!(i0(8.0 + 1e-16), 427.56411572180474);
        assert_eq!(i0(10.0), 2815.716628466254);
        assert_eq!(i0(20.0), 43558282.559553534);
        assert_eq!(i0(100.0), 1.0737517071310738e+42);

        assert_eq!(i0(-8.0 - 1e-16), 427.56411572180474);
        assert_eq!(i0(-10.0), 2815.716628466254);
        assert_eq!(i0(-20.0), 43558282.559553534);
        assert_eq!(i0(-100.0), 1.0737517071310738e+42);
    }
}

#[cfg(test)]
mod i0e_tests {
    use super::*;

    #[test]
    fn i0e_trivials() {
        assert_eq!(i0e(f64::INFINITY), 0.0);
        assert_eq!(i0e(-f64::INFINITY), 0.0);
        assert_eq!(i0e(f64::NAN).is_nan(), true);
    }

    #[test]
    fn i0e_small_x() {
        assert_eq!(i0e(0.0), 1.0);

        assert_eq!(i0e(0.1), 0.9071009257823011);
        assert_eq!(i0e(1.0), 0.46575960759364043);
        assert_eq!(i0e(4.0), 0.20700192122398672);
        assert_eq!(i0e(7.0), 0.15373774467288123);
        assert_eq!(i0e(8.0), 0.1434317818568503);

        assert_eq!(i0e(-0.1), 0.9071009257823011);
        assert_eq!(i0e(-1.0), 0.46575960759364043);
        assert_eq!(i0e(-4.0), 0.20700192122398672);
        assert_eq!(i0e(-7.0), 0.15373774467288123);
        assert_eq!(i0e(-8.0), 0.1434317818568503);
    }

    #[test]
    fn i0e_large_x() {
        assert_eq!(i0e(8.0 + 1e-16), 0.1434317818568503);
        assert_eq!(i0e(10.0), 0.1278333371634286);
        assert_eq!(i0e(20.0), 0.089780311884826);
        assert_eq!(i0e(100.0), 0.03994437929909668);

        assert_eq!(i0e(-8.0 - 1e-16), 0.1434317818568503);
        assert_eq!(i0e(-10.0), 0.1278333371634286);
        assert_eq!(i0e(-20.0), 0.089780311884826);
        assert_eq!(i0e(-100.0), 0.03994437929909668);

        // TODO: Accuracy dwindles for large x
        //assert_eq!(i0e(1000.0), 0.03994437929909668);
        //assert_eq!(i0e(1e10), 3.9894228040641945e-06);
    }
}