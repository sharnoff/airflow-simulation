//! Tools for generating `BranchTree`s

use crate::{Float, Point};

pub mod basic;
pub mod equal;

/// Helper trait to abstract away some of the parameterization around generating branches
pub trait BranchGenerator {
    /// Makes a pair of child branches for the given parent branch, using whatever information is
    /// required from the generator itself
    ///
    /// The returned tuple should be: `(left child, right child)` in that order.
    fn make_children(&self, parent: ParentInfo, depth: usize) -> (ChildInfo, ChildInfo);
}

/// Constructs the standard initial parent to use for generating [`BranchTree`]s
///
/// The "standard" initial parent always points downwards from the bottom of the usual 1x1 square
/// from (0,0) to (1,1). The length and stem radius can be provided, because they are usually
/// variable -- unlike the starting position & angle.
///
/// [`BranchTree`]: crate::BranchTree
pub fn standard_init_parent(length: Float, stem_radius: Float) -> ParentInfo {
    ParentInfo {
        // The initial point is at the middle-top of the square from (0,0) to (1,1)
        pos: Point { x: 0.5, y: 1.0 },
        // Angle points downwards -- minus pi/2
        total_angle: -crate::float::FRAC_PI_2,
        length,
        stem_radius,
    }
}

/// The necessary information about a parent branch required in order to generate its children
#[derive(Copy, Clone, Debug)]
pub struct ParentInfo {
    /// The position of the end of the parent that the child will attach to
    pub pos: Point,
    /// The angle of the parent stem, as radians anti-clockwise from the positive X direction
    pub total_angle: Float,
    pub stem_radius: Float,
    /// The full length of the parent stem
    pub length: Float,
}

/// Information about a child
#[derive(Copy, Clone, Debug)]
pub struct ChildInfo {
    /// The change in angle from the parent to this child
    pub angle_from_parent: Float,
    pub length: Float,
    pub stem_radius: Float,
    /// The radius of the sack at the end of the stem, iff this is a terminal branch
    pub sack_radius: Option<Float>,
}

impl ChildInfo {
    /// Given this child's `ParentInfo`, create the `ParentInfo` to continue generating from this
    /// child
    ///
    /// Returns `None` if the child is a terminal branch
    pub fn as_parent(&self, childs_parent: ParentInfo) -> Option<ParentInfo> {
        if self.sack_radius.is_some() {
            return None;
        }

        let total_angle = childs_parent.total_angle + self.angle_from_parent;

        Some(ParentInfo {
            pos: Point {
                // Total angles relative to the coordinate plane are standard here, so the trig is
                // recognizable.
                //
                // x + L*cos(θ)
                x: childs_parent.pos.x + self.length * total_angle.cos(),
                // y + L*sin(θ)
                y: childs_parent.pos.y + self.length * total_angle.sin(),
            },
            total_angle,
            stem_radius: self.stem_radius,
            length: self.length,
        })
    }
}
