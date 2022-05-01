//! A handful of helper functions for working with vectors

use crate::Float;

/// Returns the squared norm of the vector (i.e. || vec ||²)
pub fn norm_squared(vec: &[Float]) -> Float {
    vec.iter().map(|x| x * x).sum()
}

/// Subtracts `other` from `this`, assigning all of the values to the entries in `this`
///
/// ## Panics
///
/// This method panics if the lengths of `this` and `other` are not equal.
pub fn sub_assign(this: &mut [Float], other: &[Float]) {
    assert_eq!(this.len(), other.len());

    this.iter_mut()
        .zip(other.iter())
        .for_each(|(t, o)| *t -= *o);
}
