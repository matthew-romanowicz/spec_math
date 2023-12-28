# spec_math

This crate contains approximations for a set of mathematical functions 
commonly referred to as "special functions". The goal of this crate is
to eventually contain a full Rust re-implementation of the CEPHES 
library with minimal changes in modules such as the cephes64
module. The special functions will also be available as traits for easier
use in Rust, potentially with modifications to correct bugs or improve performance
from the original CEPHES implementation.

Currently, the gamma functions, error functions, elliptic integrals, and 
bessel functions are implemented.