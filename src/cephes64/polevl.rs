/*                                                     polevl.c
 *                                                     p1evl.c
 *
 *     Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 * The function p1evl() assumes that c_N = 1.0 so that coefficent
 * is omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */

/*
 * Cephes Math Library Release 2.1:  December, 1988
 * Copyright 1984, 1987, 1988 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Sources:
 * [1] Holin et. al., "Polynomial and Rational Function Evaluation",
 *     https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots/rational.html
 */

/* Scipy changes:
 * - 06-23-2016: add code for evaluating rational functions
 */

// from https://github.com/scipy/scipy/blob/c4ce0c4560bc635867512c4d2ea6db6f666d3eeb/scipy/special/cephes/polevl.h#L67
pub fn polevl(x: f64, coef: &[f64], n: usize) -> f64 {
    let mut ans = coef[0];
    for i in 1..(n + 1) {
        ans = ans * x + coef[i];
    }
    ans
}

/*                                                     p1evl() */
/*                                          N
* Evaluate polynomial when coefficient of x  is 1.0.
* That is, C_{N} is assumed to be 1, and that coefficient
* is not included in the input array coef.
* coef must have length N and contain the polynomial coefficients
* stored as
*     coef[0] = C_{N-1}
*     coef[1] = C_{N-2}
*          ...
*     coef[N-2] = C_1
*     coef[N-1] = C_0
* Otherwise same as polevl.
*/

pub fn p1evl(x: f64, coef: &[f64], n: usize) -> f64
{
    let mut ans = x + coef[0];

    for i in 1..n {
        ans = ans * x + coef[i];
    }
    ans
}

/* Evaluate a rational function. See [1]. */

pub fn ratevl(x: f64, num: &[f64], m: isize, denom: &[f64], n: isize) -> f64 {

    let absx = x.abs();

    let (dir, mut p, y) = if absx > 1.0 {
        /* Evaluate as a polynomial in 1/x. */
        (-1_isize, m, 1.0 / x)
    } else {
        (1_isize, 0_isize, x)
    };

    /* Evaluate the numerator */
    let mut num_ans = num[p as usize];
    p += dir;
    for _ in 1..=m {
        num_ans = num_ans * y + num[p as usize];
        p += dir;
    }

    /* Evaluate the denominator */
    if absx > 1.0 {
        p = n;
    } else {
        p = 0;
    }

    let mut denom_ans = denom[p as usize];
    p += dir;
    for _ in 1..=n {
        denom_ans = denom_ans * y + denom[p as usize];
        p += dir;
    }

    if absx > 1.0 {
        let i = n - m;
        x.powi(i as i32) * num_ans / denom_ans
    } else {
        num_ans / denom_ans
    }
}

#[cfg(test)]
mod polevl_tests {
    use super::*;

    #[test]
    fn polevl_test() {
        let coefs = [2.0, 3.0, 4.0];
        assert_eq!(polevl(1.0, &coefs, 1), 5.0);
        assert_eq!(polevl(1.0, &coefs, 2), 9.0);
        assert_eq!(polevl(2.0, &coefs, 2), 18.0);
    }
}