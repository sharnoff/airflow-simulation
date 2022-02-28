//! Wrapper module for [`BranchTree`]

use std::ops::{Index, IndexMut};

use crate::gen::{self, BranchGenerator};
use crate::{
    AcinarRegion, Bifurcation, Branch, BranchId, BranchKind, Float, Point, Tube,
    UNSET_END_PRESSURE, UNSET_FLOW_RATE,
};

/// Wrapper around a `Vec<Branch>` to provide nicer access to the internals
#[derive(Debug, Clone, PartialEq)]
pub struct BranchTree {
    bifurcations: Vec<Branch>,
    acinars: Vec<Branch>,
    root_id: BranchId,
    root_start_pos: Point,
}

impl Index<BranchId> for BranchTree {
    type Output = Branch;

    fn index(&self, id: BranchId) -> &Branch {
        match id.deconstruct() {
            (BranchKind::Bifurcation, idx) => &self.bifurcations[idx],
            (BranchKind::Acinar, idx) => &self.acinars[idx],
        }
    }
}

impl IndexMut<BranchId> for BranchTree {
    fn index_mut(&mut self, id: BranchId) -> &mut Branch {
        match id.deconstruct() {
            (BranchKind::Bifurcation, idx) => &mut self.bifurcations[idx],
            (BranchKind::Acinar, idx) => &mut self.acinars[idx],
        }
    }
}

impl BranchTree {
    /// Returns the [`BranchId`] of the root branch
    pub fn root_id(&self) -> BranchId {
        self.root_id
    }

    /// Returns the point at which the root branch starts (i.e. the side of the branch closest to
    /// the trachea)
    pub fn root_start_pos(&self) -> Point {
        self.root_start_pos
    }

    /// Returns the total number of branches stored within the tree
    pub fn count_branches(&self) -> usize {
        self.bifurcations.len() + self.acinars.len()
    }

    /// Returns the number of bifurcation branches in the tree
    pub fn count_bifurcations(&self) -> usize {
        self.bifurcations.len()
    }

    /// Returns the number of acinar branches in the tree
    pub fn count_acinar_regions(&self) -> usize {
        self.acinars.len()
    }

    /// Adds a bifurcation to the tree, returning the `BranchId` that now refers to it
    fn push_bifurcation(&mut self, branch: Bifurcation) -> BranchId {
        let idx = self.bifurcations.len();
        self.bifurcations.push(Branch::Bifurcation(branch));

        BranchId::new(BranchKind::Bifurcation, idx)
    }

    /// Adds an acinar region to the tree, returning the `BranchId` that now refers to it
    fn push_acinar(&mut self, branch: AcinarRegion) -> BranchId {
        let idx = self.acinars.len();
        self.acinars.push(Branch::Acinar(branch));

        BranchId::new(BranchKind::Acinar, idx)
    }

    /// Creates a `BranchTree` from a generator
    ///
    /// Whether the generator is seeded (and/or uses the same seed) is up to its own implementation
    /// details.
    ///
    /// ## Panics
    ///
    /// The value of `start.id` MUST be zero, and will panic if this is not the case.
    pub fn from_generator(
        start: gen::ParentInfo,
        gen: &impl BranchGenerator,
        degraded: bool,
    ) -> Self {
        use gen::{ChildInfo, ParentInfo};

        assert_eq!(start.id, 0);

        // Recursive helper function for generating the tree.
        //
        // Returns (left child, right child)
        fn full_gen(
            tree: &mut BranchTree,
            depth: usize,
            parent: ParentInfo,
            parent_id: &mut usize,
            gen: &impl BranchGenerator,
            degraded: bool,
        ) -> (BranchId, BranchId) {
            let (left, right) = gen.make_children(parent, depth, degraded);

            // Helper closure to fill out a child. Makes a recursive call to `full_gen` and places
            // the child branch into `tree.items`
            let mut make_child = |info: ChildInfo, parent_id: &mut usize| -> BranchId {
                let tube = Tube {
                    angle_from_parent: info.angle_from_parent,
                    radius: info.tube_radius,
                    length: info.length,
                    flow_in: 0.0,
                    end_pressure: 0.0,
                };

                if let Some(as_parent) = info.as_parent(parent, parent_id) {
                    let (child_left, child_right) =
                        full_gen(tree, depth + 1, as_parent, parent_id, gen, degraded);
                    let branch = Bifurcation {
                        tube,
                        left_child: child_left,
                        right_child: child_right,
                    };

                    tree.push_bifurcation(branch)
                } else {
                    let branch = AcinarRegion {
                        tube,
                        volume: 0.0,
                        compliance: info.compliance.unwrap(),
                    };

                    tree.push_acinar(branch)
                }
            };

            (make_child(left, parent_id), make_child(right, parent_id))
        }

        let mut tree = BranchTree {
            bifurcations: Vec::new(),
            acinars: Vec::new(),
            // Set an initially invalid `BranchId`; we'll replace this later. we need to use
            // isize::MAX because `BranchId`s aren't allowed to be greater than it
            root_id: BranchId::new(BranchKind::Bifurcation, isize::MAX as usize),
            root_start_pos: start.pos,
        };
        let mut parent_id = start.id;
        let (left, right) = full_gen(&mut tree, 1, start, &mut parent_id, gen, degraded);

        let root = Bifurcation {
            tube: Tube {
                // Because it's the first branch, its full angle needs to be represented by the
                // angle from its "parent" -- even though that doesn't exist.
                angle_from_parent: start.total_angle,
                radius: start.tube_radius,
                length: start.length,
                flow_in: UNSET_FLOW_RATE,
                end_pressure: UNSET_END_PRESSURE,
            },
            left_child: left,
            right_child: right,
        };

        tree.root_id = tree.push_bifurcation(root);
        tree
    }

    fn structure_state_len(&self) -> usize {
        self.bifurcations.len() + 2 * self.acinars.len()
    }

    /// Returns a vector corresponding to the state of the tree that is subject to degration
    ///
    /// Currently, this is just branch radii and compliance for each acinar region.
    pub fn structure_state(&self) -> Vec<Float> {
        let mut state = vec![0.0; self.structure_state_len()];

        let mut i = 0;
        for b in self.bifurcations.iter().chain(self.acinars.iter()) {
            match b {
                Branch::Bifurcation(Bifurcation { tube, .. }) => {
                    state[i] = tube.radius;
                    i += 1;
                }
                Branch::Acinar(AcinarRegion {
                    tube, compliance, ..
                }) => {
                    state[i] = tube.radius;
                    i += 1;
                    state[i] = *compliance;
                    i += 1;
                }
            }
        }

        assert_eq!(i, state.len());

        state
    }

    /// Sets the "structure" of the tree according, using a vector that matches the return of
    /// [`structure_state`]
    ///
    /// ## Panics
    ///
    /// This method will panic if `state` does not have the same length as what [`structure_state`]
    /// returns.
    ///
    /// [`structure_state`]: Self::structure_state
    pub fn set_structure(&mut self, state: &[Float]) {
        assert_eq!(state.len(), self.structure_state_len());

        let mut i = 0;
        for b in self.bifurcations.iter_mut().chain(self.acinars.iter_mut()) {
            match b {
                Branch::Bifurcation(Bifurcation { tube, .. }) => {
                    tube.radius = state[i];
                    i += 1;
                }
                Branch::Acinar(AcinarRegion {
                    tube, compliance, ..
                }) => {
                    tube.radius = state[i];
                    i += 1;
                    *compliance = state[i];
                    i += 1;
                }
            }
        }

        assert_eq!(i, state.len());
    }

    #[cfg(test)]
    pub fn branches_mut(&mut self) -> impl '_ + Iterator<Item = &'_ mut Branch> {
        self.bifurcations.iter_mut().chain(self.acinars.iter_mut())
    }
}
