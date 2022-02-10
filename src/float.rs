//! Wrapper module to allow switching the float type globally.
//!
//! All of the functions from `f32`/`f64` are glob imported into this module, so that they can be
//! used elsewhere.

/// Type alias for the selected global float type
pub type Float = f64;

// import everything to do with the float
pub use std::f64::consts::*;
pub use std::f64::*;
