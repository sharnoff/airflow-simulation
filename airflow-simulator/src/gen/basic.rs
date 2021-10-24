//! Basic, proof-that-it-works generator

use super::*;
use crate::float;

/// A branch generator that always produces mirrored terminal children
#[allow(unused)]
pub struct MirroredTerminalChildGenerator;

impl BranchGenerator for MirroredTerminalChildGenerator {
    fn make_children(&self, parent: ParentInfo, _depth: usize) -> (ChildInfo, ChildInfo) {
        // The total cross-sectional area must be the same, so if the child radii are equal, they
        // must each be `1/sqrt(2)` times the parent's radius.
        let stem_radius = float::FRAC_1_SQRT_2 * parent.stem_radius;

        // Arbitrary. Chosen to be similar to above so that scaling looks even.
        let length = float::FRAC_1_SQRT_2 * parent.length;
        // each child is 45° off from the parent
        let base_angle = float::FRAC_PI_4;

        // Also arbitrary.
        let sack_radius = 2.5 * stem_radius;

        let left = ChildInfo {
            stem_radius,
            length,
            angle_from_parent: -base_angle,
            sack_radius: Some(sack_radius),
        };

        // The right child is the same except for the angle
        let right = ChildInfo {
            angle_from_parent: base_angle,
            ..left
        };

        (left, right)
    }
}
