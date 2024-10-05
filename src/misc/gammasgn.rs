pub fn gammasgn(x: f64) -> f64
{

    if x.is_nan() {
      x
    } else if x > 0.0 {
        1.0
    } else {
        let fx = x.floor();
        if x - fx == 0.0 {
            0.0
        } else if fx as isize % 2 != 0 {
            -1.0
        } else {
            1.0
        }
    }
}