// from https://github.com/scipy/scipy/blob/c4ce0c4560bc635867512c4d2ea6db6f666d3eeb/scipy/special/cephes/polevl.h#L67
pub fn polevl(x: f64, coef: &[f64], n: usize) -> f64 {
    let mut ans = coef[0];
    for i in 1..(n + 1) {
        ans = ans * x + coef[i];
    }
    ans
}