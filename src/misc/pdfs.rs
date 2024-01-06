use crate::cephes64::{beta, lbeta, rgamma, lgam};

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

pub fn nbinom_pmf(k: i32, n: i32, p: f64) -> f64 {
    //! Negative binomial probability mass function

    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }

    // (k + (n - 1)) * (1 - p)^k * p^n / (k*(n - 1)*beta(k, n - 1))

    let log_pmf = ((k + n - 1) as f64).ln() + k as f64 * (1.0 - p).ln() + n as f64 * p.ln() 
        - (k as f64 * (n as f64 - 1.0)).ln() - lbeta(k as f64, n as f64 - 1.0);

    log_pmf.exp()
}

pub fn pois_pmf(k: i32, m: f64) -> f64 {
    //! Poisson probability mass function

    if m < 0.0 || k < 0 {
        return f64::NAN;
    }

    // (k + (n - 1)) * (1 - p)^k * p^n / (k*(n - 1)*beta(k, n - 1))

    let log_pmf = (k as f64) * m.ln() - m - lgam(k as f64 + 1.0);

    log_pmf.exp()
}

pub fn beta_pdf(a: f64, b: f64, x: f64) -> f64 {
    //! Beta probability density function

    if !(0.0..=1.0).contains(&x) {
        return f64::NAN;
    }

    let den = beta(a, b);
    if !den.is_infinite() && den != 0.0 {
        x.powf(a - 1.0) * (1.0 - x).powf(b - 1.0) / den
    } else {
        // Use lbeta to avoid overflow
        let den = lbeta(a, b);
        ((a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() - den).exp()
    }
}

pub fn chi2_pdf(df: f64, x: f64) -> f64 {
    //! Chi-squared probability density function

    if x <= 0.0 {
        if x < 0.0 || (x == 0.0 && df == 1.0) {
            f64::NAN
        } else if df == 2.0 {
            0.5
        } else {
            0.0
        }
    } else if df <= 300.0 { // rgamma will start returning 0 after this value
        (0.5 * x).powf(0.5 * df) * rgamma(0.5 * df) * (-0.5 * x).exp() / x
    } else {
        (0.5 * df * (0.5 * x).ln() - 0.5 * x - x.ln() - lgam(0.5 * df)).exp()
    }
}

pub fn f_pdf(d1: f64, d2: f64, x: f64) -> f64 {
    //! F distribution probability density function

    if x <= 0.0 {
        if x < 0.0 || (x == 0.0 && d1 == 1.0) {
            f64::NAN
        } else {
            0.0
        }
    } else {
        let log_pdf = d1 * 0.5 * (d1).ln() + d2 * 0.5 * d2.ln() + (d1 * 0.5 - 1.0)*x.ln() - 
        ((d1 + d2) * 0.5 * (d1 * x + d2).ln() + lbeta(0.5 * d1, 0.5 * d2));
        log_pdf.exp()
    }
}

pub fn gamma_pdf(a: f64, b: f64, x: f64) -> f64 {
    //! Gamma distribution probability density function

    // TODO: Add tests for this

    if x <= 0.0 {
        f64::NAN
    } else {
        let log_pdf = b * a.ln() + (b - 1.0) * x.ln() - a * x - lgam(b);
        log_pdf.exp()
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
    fn binom_pmf_trivials() {
        assert_eq!(binom_pmf(10.0, 15, 0.0), 0.0);
        assert_eq!(binom_pmf(10.0, 15, 1.0), 0.0);
        assert_eq!(binom_pmf(10.0, 15, -1e-10).is_nan(), true);
        assert_eq!(binom_pmf(10.0, 15, 1.0 + 1e-10).is_nan(), true);
    }

    #[test]
    fn binom_pmf_values() {
        assert_eq!(binom_pmf(3.0, 5, 0.5), 0.3125);
        assert_eq!(binom_pmf(3.5, 5, 0.5), 0.3125);
        assert_eq!(binom_pmf(150.0, 155, 0.5), 1.5294448135427591e-38);
        assert_eq!(binom_pmf(1500.0, 10000, 0.1), 8.160246663190426e-56);
        assert_eq!(binom_pmf(1500.5, 10000, 0.1), 8.160246663190426e-56);
    }
}

#[cfg(test)]
mod nbinom_pmf_tests {
    use super::*;

    #[test]
    fn nbinom_pmf_trivials() {
        assert_eq!(nbinom_pmf(10, 15, 0.0), 0.0);
        assert_eq!(nbinom_pmf(10, 15, 1.0), 0.0);
        assert_eq!(nbinom_pmf(10, 15, -1e-10).is_nan(), true);
        assert_eq!(nbinom_pmf(10, 15, 1.0 + 1e-10).is_nan(), true);
    }

    #[test]
    fn nbinom_pmf_values() {
        assert_eq!(nbinom_pmf(3, 5, 0.5), 0.13671874999999994);
        assert_eq!(nbinom_pmf(150, 155, 0.5), 0.022270220044768158);
        assert_eq!(nbinom_pmf(150, 100, 0.01), 5.370670528582223e-130);
    }
}

#[cfg(test)]
mod pois_pmf_tests {
    use super::*;

    #[test]
    fn pois_pmf_trivials() {
        assert_eq!(pois_pmf(5, 0.0), 0.0);
        assert!(pois_pmf(5, -1e-10).is_nan());
        assert!(pois_pmf(-1, 0.5).is_nan());
    }

    #[test]
    fn pois_pmf_values() {
        assert_eq!(pois_pmf(0, 0.5), 0.6065306597126334);
        assert_eq!(pois_pmf(5, 0.5), 0.0001579506926334984);
        assert_eq!(pois_pmf(150, 0.5), 7.43806550111441e-309);
        assert_eq!(pois_pmf(50, 0.01), 3.255233773313457e-165);
        assert_eq!(pois_pmf(50, 1e2), 1.2231421635189012e-08);
    }
}

#[cfg(test)]
mod beta_pdf_tests {
    use super::*;

    #[test]
    fn beta_pdf_trivials() {
        assert_eq!(beta_pdf(10.0, 15.0, 0.0), 0.0);
        assert_eq!(beta_pdf(10.0, 15.0, 1.0), 0.0);
        assert_eq!(beta_pdf(10.0, 15.0, -1e-10).is_nan(), true);
        assert_eq!(beta_pdf(10.0, 15.0, 1.0 + 1e-10).is_nan(), true);
    }

    #[test]
    fn beta_pdf_values() {
        assert_eq!(beta_pdf(3.0, 5.0, 0.5), 1.6406249999999998);
        assert_eq!(beta_pdf(3.5, 5.0, 0.5), 1.9440689061613787);
        assert_eq!(beta_pdf(150.0, 155.0, 0.5), 13.362132026865142);
        assert_eq!(beta_pdf(1500.0, 10000.0, 0.1), 3.293893418465068e-22);
        assert_eq!(beta_pdf(1500.5, 10000.0, 0.1), 2.884325727129547e-22);
    }
}

#[cfg(test)]
mod chi2_pdf_tests {
    use super::*;

    #[test]
    fn chi2_pdf_trivials() {
        assert_eq!(chi2_pdf(5.0, -1e-20).is_nan(), true);
        assert_eq!(chi2_pdf(1.0, 0.0).is_nan(), true);
        assert_eq!(chi2_pdf(2.0, 0.0), 0.5);
        assert_eq!(chi2_pdf(3.0, 0.0), 0.0);
    }

    #[test]
    fn chi2_pdf_values() {
        assert_eq!(chi2_pdf(0.1, 2.0), 0.009447299158950628);
        assert_eq!(chi2_pdf(1.0, 2.0), 0.10377687435514868);
        assert_eq!(chi2_pdf(100.0, 2.0), 3.023922484977498e-64);

        // TODO: Check to see how accurate these are
        assert_eq!(chi2_pdf(1e5, 1e5), 0.0008920605712562789);
        assert_eq!(chi2_pdf(1e10, 1e10), 2.82093439313804e-6);
    }
}

#[cfg(test)]
mod f_pdf_tests {
    use super::*;

    #[test]
    fn f_pdf_trivials() {
        assert_eq!(f_pdf(1.0, 1.0, 0.0).is_nan(), true);
        assert_eq!(f_pdf(2.0, 1.0, -1e-10).is_nan(), true);
        assert_eq!(f_pdf(2.0, 1.0, 0.0), 0.0);
    }

    #[test]
    fn f_pdf_values() {
        assert_eq!(f_pdf(2.0, 2.0, 2.0), 0.11111111111111115);
        assert_eq!(f_pdf(1e10, 1e-5, 1e10), 4.999187035248361e-16);
    }
}

// #[cfg(test)]
// mod gamma_pdf_tests {
//     use super::*;

//     #[test]
//     fn gamma_pdf_trivials() {
//         assert_eq!(gamma_pdf(1.0, 1.0, 0.0).is_nan(), true);
//         assert_eq!(gamma_pdf(2.0, 1.0, -1e-10).is_nan(), true);
//         assert_eq!(gamma_pdf(2.0, 1.0, 0.0), 0.0);
//     }

//     #[test]
//     fn gamma_pdf_values() {
//         assert_eq!(gamma_pdf(2.0, 2.0, 2.0), 0.29305022221974686);
//         //assert_eq!(gamma_pdf(1e-5, 1e10, 1e10), 4.999187035248361e-16);
//     }
// }