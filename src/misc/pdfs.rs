const SQRT_2PI_RECIP: f64 = 0.3989422804014327;

pub fn norm_pdf(x: f64) -> f64 {
    //! Normal (Gaussian) probability density function
    SQRT_2PI_RECIP * (-0.5 * x * x).exp()
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