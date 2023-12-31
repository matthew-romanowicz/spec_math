//! Functions from the CEPHES library re-implemented in Rust for f64.
//!
//! # Description:
//! This module contains approximations of special functions exactly as seen in the CEPHES library. Algorithms are re-implemented in Rust but are unchanged.

mod consts;
pub use crate::cephes64::consts::*;
mod cbrt;
pub use crate::cephes64::cbrt::cbrt;
mod sinpi;
pub use crate::cephes64::sinpi::{sinpi, cospi};
mod polevl;
pub use crate::cephes64::polevl::{polevl, p1evl};
mod dawsn;
pub use crate::cephes64::dawsn::dawsn;
mod fresnl;
pub use crate::cephes64::fresnl::fresnl;
mod sici;
pub use crate::cephes64::sici::sici;
mod ndtr;
pub use crate::cephes64::ndtr::{ndtr, erfc, erf};
mod ndtri;
pub use crate::cephes64::ndtri::ndtri;
mod gamma;
//use crate::cephes64::gamma::{lgam_sgn};
pub use crate::cephes64::gamma::{gamma, lgam};
mod rgamma;
pub use crate::cephes64::rgamma::rgamma;
mod psi;
pub use crate::cephes64::psi::psi;
mod beta;
pub use crate::cephes64::beta::{beta, lbeta};
mod incbet;
pub use crate::cephes64::incbet::incbet;
mod incbi;
pub use crate::cephes64::incbi::incbi;
mod btdtr;
pub use crate::cephes64::btdtr::btdtr;
mod bdtr;
pub use crate::cephes64::bdtr::{bdtrc, bdtr, bdtri};
mod zeta;
pub use crate::cephes64::zeta::zeta;
mod zetac;
pub use crate::cephes64::zetac::{zetac, riemann_zeta};
mod chbevl;
use crate::cephes64::chbevl::chbevl;
mod unity;
pub use crate::cephes64::unity::{log1p, log1pmx, expm1, cosm1, lgam1p};
mod lanczos;
mod igam;
pub use crate::cephes64::igam::{igam, igamc};
mod igami;
pub use crate::cephes64::igami::{igami, igamci};
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
mod airy;
pub use crate::cephes64::airy::airy;
mod j0;
pub use crate::cephes64::j0::{j0, y0};
mod j1;
pub use crate::cephes64::j1::{j1, y1};
mod jv;
pub use crate::cephes64::jv::jv;
mod yn;
pub use crate::cephes64::yn::yn;
mod yv;
pub use crate::cephes64::yv::yv;
mod i0;
pub use crate::cephes64::i0::{i0, i0e};
mod i1;
pub use crate::cephes64::i1::{i1, i1e};
mod k0;
pub use crate::cephes64::k0::{k0, k0e};
mod k1;
pub use crate::cephes64::k1::{k1, k1e};
mod kn;
pub use crate::cephes64::kn::kn;