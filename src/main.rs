//! Airflow simulator for my 3rd-year computer science project

use std::ops::{Index, IndexMut};

mod branch_id;
mod float;
mod gen;
mod img;
mod sim;

use branch_id::BranchId;
use float::Float;
use gen::equal::EqualChildGenerator;
use gen::BranchGenerator;
use img::{rgb, rgba, ImageConfig};

/// A physical point in 2D space, used to represent one end of a branch
///
/// The values are in units of millimeters
///
/// We treat positive X as to the right and positive Y as up, with both values greater than zero.
/// This distinction matters when we're rendering an image of the tree.
///
/// It's also worth noting that x and y values are typically (though not *necessarily*) greater
/// than zero.
#[derive(Copy, Clone, Debug)]
pub struct Point {
    x: Float,
    y: Float,
}

/// Marker for the type of a branch. Primarily used with [`BranchId`] methods.
#[derive(Copy, Clone, Debug)]
enum BranchKind {
    Bifurcation = 0,
    Acinar = 1,
}

/// An individual node in the tree of branches
///
/// Branches are tyically referred to by their [`BranchId`]. Both `Bifurcation` and `Acinar` share
/// a `tube` field representing their bronchiole, which can be retrieved with the [`tube`] method.
///
/// [`stem`]: Self::stem
#[derive(Copy, Clone, Debug)]
enum Branch {
    Bifurcation(Bifurcation),
    Acinar(AcinarRegion),
}

/// The default value for `Tube.flow_rate` as we're generating a [`BranchTree`]
const UNSET_FLOW_RATE: Float = 0.0;

/// The default value for `Tube.end_pressure` as we're generating a [`BranchTree`]
const UNSET_END_PRESSURE: Float = 0.0;

#[derive(Copy, Clone, Debug)]
struct Tube {
    /// The ange of the tube, relative to its parent. Measured in radians, always within 0..2π
    angle_from_parent: Float,

    /// The radius of the tube, in mm
    radius: Float,
    /// The length of the tube, in mm
    length: Float,
    /// The velocity at which air is flowing out of the tube, towards the trachea (negative if
    /// breathing in)
    ///
    /// Units of mm/s
    flow_rate: Float,
    /// The pressure at the distal end of the tube (i.e. away from the trachea), in MegaPascals
    ///
    /// We don't need to record the pressure at the other end because that's given by this one's
    /// parent (or the root pressure, if this has no parent)
    end_pressure: Float,
}

/// A `Node` that represents a branch with two children
///
/// The directions for "left" and "right" are as if the bifurcation is branching downwards -- the
/// left child corresponds to a negative rotation from where it is attached, and the right child a
/// positive one.
#[derive(Copy, Clone, Debug)]
struct Bifurcation {
    /// Information about the bronchus/bronchiole going into this bifurcation
    tube: Tube,

    left_child: BranchId,
    right_child: BranchId,
}

#[derive(Copy, Clone, Debug)]
struct AcinarRegion {
    /// Information about the final bonchiole leading into this acinar region
    tube: Tube,

    /// The "compliance" of the region -- essentially how much it is willing to stretch
    ///
    /// Essentially volume / pressure; units of mm³/MPa
    compliance: Float,
    /// Total volume of the region; units of mm³/MPa
    volume: Float,
}

impl Branch {
    /// Returns a reference to the shared `Tube` field of the branch
    fn tube(&self) -> &Tube {
        match self {
            Branch::Bifurcation(b) => &b.tube,
            Branch::Acinar(b) => &b.tube,
        }
    }

    /// Returns a mutable reference to the shared `Tube` field of the branch
    fn tube_mut(&mut self) -> &mut Tube {
        match self {
            Branch::Bifurcation(b) => &mut b.tube,
            Branch::Acinar(b) => &mut b.tube,
        }
    }
}

/// Wrapper around a `Vec<Branch>` to provide nicer access to the internals
#[derive(Debug, Clone)]
pub struct BranchTree {
    bifurcations: Vec<Branch>,
    acinars: Vec<Branch>,
    root_id: BranchId,
    root_start_pos: Point,
    // The total number of terminal (Acinar) branches in the tree
    //
    // This is not necessarily items.len() / 2, because we might have unbalanced trees.
    num_terminal: usize,
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

    /// Creates a `BranchTree` from a random generator
    ///
    /// Whether the generator is seeded (and/or uses the same seed) is up to its own implementation
    /// details.
    fn from_generator(start: gen::ParentInfo, gen: &impl BranchGenerator) -> Self {
        use gen::{ChildInfo, ParentInfo};

        // Recursive helper function for generating the tree.
        //
        // Returns (left child, right child)
        fn full_gen(
            tree: &mut BranchTree,
            depth: usize,
            parent: ParentInfo,
            gen: &impl BranchGenerator,
        ) -> (BranchId, BranchId) {
            let (left, right) = gen.make_children(parent, depth);

            // Helper closure to fill out a child. Makes a recursive call to `full_gen` and places
            // the child branch into `tree.items`
            let mut make_child = |info: ChildInfo| -> BranchId {
                let tube = Tube {
                    angle_from_parent: info.angle_from_parent,
                    radius: info.tube_radius,
                    length: info.length,
                    flow_rate: 0.0,
                    end_pressure: 0.0,
                };

                if let Some(as_parent) = info.as_parent(parent) {
                    let (child_left, child_right) = full_gen(tree, depth + 1, as_parent, gen);
                    let branch = Bifurcation {
                        tube,
                        left_child: child_left,
                        right_child: child_right,
                    };

                    tree.push_bifurcation(branch)
                } else {
                    tree.num_terminal += 1;
                    let branch = AcinarRegion {
                        tube,
                        volume: 0.0,
                        compliance: info.compliance.unwrap(),
                    };

                    tree.push_acinar(branch)
                }
            };

            (make_child(left), make_child(right))
        }

        let mut tree = BranchTree {
            bifurcations: Vec::new(),
            acinars: Vec::new(),
            // Set an initially invalid `BranchId`; we'll replace this later.
            root_id: BranchId::new(BranchKind::Bifurcation, usize::MAX),
            root_start_pos: start.pos,
            num_terminal: 0,
        };
        let (left, right) = full_gen(&mut tree, 1, start, gen);

        let root = Bifurcation {
            tube: Tube {
                // Becuase it's the first branch, its full angle needs to be represented by the
                // angle from its "parent" -- even though that doesn't exist.
                angle_from_parent: start.total_angle,
                radius: start.tube_radius,
                length: start.length,
                flow_rate: UNSET_FLOW_RATE,
                end_pressure: UNSET_END_PRESSURE,
            },
            left_child: left,
            right_child: right,
        };

        tree.root_id = tree.push_bifurcation(root);
        tree
    }
}

/// Atmospheric pressure at sea level, in MegaPascals
const ATMOSPHERIC_PRESSURE: Float = 0.101_325;

fn main() {
    let start = gen::ParentInfo {
        pos: Point { x: 0.5, y: 1.0 },
        // Angle points downwards -- minus pi/2
        total_angle: -float::FRAC_PI_2,
        length: 0.3,
        tube_radius: 0.05,
    };
    let img_config = ImageConfig {
        centered_at: Point { x: 0.5, y: 0.5 },
        width: 550,
        height: 550,
        scale: 500.0,
        grid_lines: vec![
            // X & Y axes:
            (rgb(0x000000), img::Axis::X, 0.0, 5),
            (rgb(0x000000), img::Axis::Y, 0.0, 5),
        ],
        background: rgba(0x00000000),
        stem_color: rgb(0xFFFFFF),
        sack_color: rgb(0xFF0000),
    };

    let generator = EqualChildGenerator {
        max_depth: 6,
        // 0.2L/cmH₂O is approximately 0.002 mm³/MPa
        compliance: 0.002,
    };
    let mut tree = BranchTree::from_generator(start, &generator);
    let sim_env = sim::SimulationEnvironment {
        pleural_pressure: 0.0,
        tracheal_pressure: ATMOSPHERIC_PRESSURE,
    };

    sim_env.reset(&mut tree);

    img_config
        .make_image(&tree)
        .save("branch-tree.png")
        .expect("failed to write the image");
}
