use crate::cephes64::{beta, lbeta};

const SQRT_2PI_RECIP: f64 = 0.3989422804014327;

pub fn norm_pdf(x: f64) -> f64 {
    //! Normal (Gaussian) probability density function
    SQRT_2PI_RECIP * (-0.5 * x * x).exp()
}

pub fn binom_pmf(k: f64, n: i32, p: f64) -> f64 {
    //! Binomial probability mass function
    
    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    let kf = k.floor();
    let c = (n as f64) / (kf * (n as f64 - kf) * beta(kf, n as f64 - kf));

    if !c.is_infinite() {
        c * p.powi(kf as i32) * (1.0 - p).powi(n - (kf as i32))
    } else {
        // Use lbeta to avoid overflow
        let c = (n as f64).ln() - kf.ln() - (n as f64 - kf).ln() - lbeta(kf, n as f64 - kf);
        (c + kf * p.ln() + (n as f64 - kf) * (1.0 - p).ln()).exp()
    }
}

#[cfg(test)]
mod norm_pdf_tests {
    use super::*;

    #[test]
    fn norm_pdf_trivials() {
        assert_eq!(norm_pdf(f64::NAN).is_nan(), true);
        assert_eq!(norm_pdf(f64::INFINITY), 0.0);
        assert_eq!(norm_pdf(-f64::INFINITY), 0.0);
    }

    #[test]
    fn norm_pdf_values() {
        assert_eq!(norm_pdf(0.0), 0.3989422804014327);

        assert_eq!(norm_pdf(1.0), 0.24197072451914337);
        assert_eq!(norm_pdf(2.0), 0.05399096651318806);
        assert_eq!(norm_pdf(3.0), 0.0044318484119380075);
        assert_eq!(norm_pdf(5.0), 1.4867195147342977e-06);
        assert_eq!(norm_pdf(10.0), 7.69459862670642e-23);
        assert_eq!(norm_pdf(20.0), 5.520948362159764e-88);
        assert_eq!(norm_pdf(35.0), 3.940396277136024e-267);
        assert_eq!(norm_pdf(40.0), 0.0);

        assert_eq!(norm_pdf(-1.0), 0.24197072451914337);
        assert_eq!(norm_pdf(-2.0), 0.05399096651318806);
        assert_eq!(norm_pdf(-3.0), 0.0044318484119380075);
        assert_eq!(norm_pdf(-5.0), 1.4867195147342977e-06);
        assert_eq!(norm_pdf(-10.0), 7.69459862670642e-23);
        assert_eq!(norm_pdf(-20.0), 5.520948362159764e-88);
        assert_eq!(norm_pdf(-35.0), 3.940396277136024e-267);
        assert_eq!(norm_pdf(-40.0), 0.0);
    }
}

#[cfg(test)]
mod binom_pmf_tests {
    use super::*;

    #[test]
    fn bindom_pmf_trivials() {
        assert_eq!(binom_pmf(10.0, 15, 0.0), 0.0);
        assert_eq!(binom_pmf(10.0, 15, 1.0), 0.0);
        assert_eq!(binom_pmf(10.0, 15, -1e-10).is_nan(), true);
        assert_eq!(binom_pmf(10.0, 15, 1.0 + 1e-10).is_nan(), true);
    }

    #[test]
    fn bindom_pmf_values() {
        assert_eq!(binom_pmf(3.0, 5, 0.5), 0.3125);
        assert_eq!(binom_pmf(3.5, 5, 0.5), 0.3125);
        assert_eq!(binom_pmf(150.0, 155, 0.5), 1.5294448135427591e-38);
        assert_eq!(binom_pmf(1500.0, 10000, 0.1), 8.160246663190426e-56);
        assert_eq!(binom_pmf(1500.5, 10000, 0.1), 8.160246663190426e-56);
    }
}