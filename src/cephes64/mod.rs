//! Functions from the CEPHES library re-implemented in Rust for f64.
//!
//! # Description:
//! This module contains approximations of special functions exactly as seen in the CEPHES library. Algorithms are re-implemented in Rust but are unchanged.

mod consts;
use crate::cephes64::consts::*;
mod polevl;
use crate::cephes64::polevl::polevl;
mod unity;
use crate::cephes64::unity::cosm1;
mod ellpj;
pub use crate::cephes64::ellpj::ellpj;
mod ellpk;
pub use crate::cephes64::ellpk::ellpk;
mod ellpe;
pub use crate::cephes64::ellpe::ellpe;
mod ellik;
pub use crate::cephes64::ellik::ellik;
mod ellie;
pub use crate::cephes64::ellie::ellie;
mod j0;
pub use crate::cephes64::j0::{j0, y0};
mod j1;
pub use crate::cephes64::j1::{j1, y1};