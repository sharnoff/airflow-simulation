//! Basic, proof-that-it-works generator

use super::*;
use crate::float;

/// A branch generator that always produces mirrored terminal children
#[allow(unused)]
pub struct MirroredTerminalChildGenerator {
    pub compliance: Float,
}

impl BranchGenerator for MirroredTerminalChildGenerator {
    fn make_children(&self, parent: ParentInfo, _depth: usize) -> (ChildInfo, ChildInfo) {
        // The total cross-sectional area must be the same, so if the child radii are equal, they
        // must each be `1/sqrt(2)` times the parent's radius.
        let tube_radius = float::FRAC_1_SQRT_2 * parent.tube_radius;

        // Arbitrary. Chosen to be similar to above so that scaling looks even.
        let length = float::FRAC_1_SQRT_2 * parent.length;
        // each child is 45° off from the parent
        let base_angle = float::FRAC_PI_4;

        let left = ChildInfo {
            tube_radius,
            length,
            angle_from_parent: -base_angle,
            compliance: Some(self.compliance),
        };

        // The right child is the same except for the angle
        let right = ChildInfo {
            angle_from_parent: base_angle,
            ..left
        };

        (left, right)
    }
}
