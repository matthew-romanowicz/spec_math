#![allow(clippy::excessive_precision)]

pub const M_2_PI: f64 = 0.63661977236758134308;  /* 2/pi */
pub const M_PI: f64 = 3.14159265358979323846;
pub const M_PI_2: f64 = 1.57079632679489661923;
pub const M_PI_4: f64 = 0.785398163397448309616;

pub const MACHEP: f64 = 1.11022302462515654042E-16;

// #ifdef DENORMAL
pub const MAXLOG: f64 = 7.09782712893383996732E2;	/* log(DBL_MAX) */
pub const MINLOG: f64 = -7.451332191019412076235E2;	/* log(2**-1075) */
// #else
// pub const MAXLOG: f64 = 7.08396418532264106224E2;	/* log 2**1022 */
// pub const MINLOG: f64 = -7.08396418532264106224E2;	/* log 2**-1022 */
// #endif

pub const SQ2OPI: f64 = 7.9788456080286535587989E-1;	/* sqrt( 2/pi ) */
pub const LOGSQ2: f64 = 3.46573590279972654709E-1;	/* log(2)/2 */
pub const THPIO4: f64 = 2.35619449019234492885;	/* 3*pi/4 */