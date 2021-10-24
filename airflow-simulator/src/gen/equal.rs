//! Generator for trees that are perfectly equal

use super::*;
use crate::float;

/// A branch generator that always produces equal branching children, up to a certain maximum depth
#[allow(unused)]
pub struct EqualChildGenerator {
    pub max_depth: usize,
}

impl BranchGenerator for EqualChildGenerator {
    fn make_children(&self, parent: ParentInfo, depth: usize) -> (ChildInfo, ChildInfo) {
        let stem_radius = float::FRAC_1_SQRT_2 * parent.stem_radius;
        // Arbitrary. A bit less than 1/sqrt(2)
        let length = 0.66 * parent.length;
        // Arbitrary. Results in 1/3 of the circle between the children
        let base_angle = float::FRAC_PI_3;
        // Arbitrary. Looks good
        let sack_radius = Some(2.5 * stem_radius);

        let can_continue = depth + 1 < self.max_depth;

        let left = ChildInfo {
            stem_radius,
            length,
            angle_from_parent: -base_angle,
            sack_radius: if can_continue { None } else { sack_radius },
        };

        let right = ChildInfo {
            angle_from_parent: base_angle,
            ..left
        };

        (left, right)
    }
}
