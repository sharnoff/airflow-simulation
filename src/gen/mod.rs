//! Tools for generating `BranchTree`s

use crate::{Float, Point};

mod equal;
mod from_json;

pub use equal::EqualChildGenerator;
pub use from_json::{FromJsonGenerator, ParsedConfig as JsonConfig};

/// Helper trait to abstract away some of the parameterization around generating branches
pub trait BranchGenerator {
    /// Makes a pair of child branches for the given parent branch, using whatever information is
    /// required from the generator itself
    ///
    /// The returned tuple should be: `(left child, right child)` in that order.
    fn make_children(&self, parent: ParentInfo, depth: usize) -> (ChildInfo, ChildInfo);
}

/// The necessary information about a parent branch required in order to generate its children
#[derive(Copy, Clone, Debug)]
pub struct ParentInfo {
    /// The "ID" of this parent branch
    ///
    /// ID numbers are assigned in a depth-first order, exploring left children first. For example:
    ///
    /// ```js
    /// {
    ///     id: 0,
    ///     left: {
    ///         id: 1,
    ///         left: { id: 2 },
    ///     },
    ///     right: {
    ///         id: 3,
    ///         right: { id: 4 },
    ///     }
    /// }
    /// ```
    pub id: usize,

    /// The position of the end of the parent that the child will attach to
    pub pos: Point,
    /// The angle of the parent stem, as radians anti-clockwise from the positive X direction
    pub total_angle: Float,
    pub tube_radius: Float,
    /// The full length of the parent stem
    pub length: Float,
}

/// Information about a child
#[derive(Copy, Clone, Debug)]
pub struct ChildInfo {
    /// The change in angle from the parent to this child
    pub angle_from_parent: Float,
    pub length: Float,
    pub tube_radius: Float,
    /// The value of `AcinarRegion.compliance` iff this is a terminal branch
    pub compliance: Option<Float>,
}

impl ChildInfo {
    /// Given this child's `ParentInfo`, create the `ParentInfo` to continue generating from this
    /// child
    ///
    /// Returns `None` if the child is a terminal branch
    pub fn as_parent(
        &self,
        childs_parent: ParentInfo,
        parent_id: &mut usize,
    ) -> Option<ParentInfo> {
        if self.compliance.is_some() {
            return None;
        }

        let total_angle = childs_parent.total_angle + self.angle_from_parent;

        let id = *parent_id + 1;
        *parent_id = id;

        Some(ParentInfo {
            id,
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
            tube_radius: self.tube_radius,
            length: self.length,
        })
    }
}
