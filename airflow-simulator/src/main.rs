//! Airflow simulator for my 3rd-year computer science project

mod float;
mod gen;
mod img;

use float::Float;
use gen::basic::MirroredTerminalChildGenerator;
use gen::BranchGenerator;
use img::{rgb, rgba, ImageConfig};

/// Unique identifier for a branch
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct BranchId(usize);

/// A point in 2D space, used to represent one end of a branch
///
/// We treat positive X as to the right and positive Y as up, generally with both values in the
/// range (0, 1).
///
/// Branch trees typically have their entrypoint at `(0.5, 1)`, with an initial angle directly
/// downwards.
#[derive(Copy, Clone, Debug)]
pub struct Point {
    x: Float,
    y: Float,
}

#[derive(Clone, Debug)]
enum Branch {
    Stem(StemBranch),
    Tail(TailBranch),
}

impl BranchId {
    /// Returns the `BranchId` as a `usize` so that it can be used as an array index
    fn index(&self) -> usize {
        self.0 as usize
    }
}

/// A branch with children
#[derive(Clone, Debug)]
struct StemBranch {
    /// The change in angle from the parent branch
    angle_from_parent: Float,
    length: Float,
    stem_radius: Float,
    left_child: BranchId,
    right_child: BranchId,
}

/// A branch that ends in a sack
#[derive(Copy, Clone, Debug)]
struct TailBranch {
    angle_from_parent: Float,
    length: Float,
    stem_radius: Float,
    sack_radius: Float,
}

impl Branch {
    /// Returns the value of the `angle_from_parent` field of the branch
    fn angle_from_parent(&self) -> Float {
        match self {
            Branch::Stem(s) => s.angle_from_parent,
            Branch::Tail(t) => t.angle_from_parent,
        }
    }

    /// Returns the value of the `length` field of the branch
    fn length(&self) -> Float {
        match self {
            Branch::Stem(s) => s.length,
            Branch::Tail(t) => t.length,
        }
    }

    /// Returns the value of the `stem_radius` field of the branch
    fn stem_radius(&self) -> Float {
        match self {
            Branch::Stem(s) => s.stem_radius,
            Branch::Tail(t) => t.stem_radius,
        }
    }
}

/// Wrapper around a `Vec<Branch>` to provide nicer access to the internals
#[derive(Debug, Clone)]
pub struct BranchTree {
    items: Vec<Branch>,
    root_id: BranchId,
    root_start_pos: Point,
}

impl BranchTree {
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
                let child_branch = if let Some(as_parent) = info.as_parent(parent) {
                    let (child_left, child_right) = full_gen(tree, depth + 1, as_parent, gen);
                    Branch::Stem(StemBranch {
                        angle_from_parent: info.angle_from_parent,
                        length: info.length,
                        stem_radius: info.stem_radius,
                        left_child: child_left,
                        right_child: child_right,
                    })
                } else {
                    Branch::Tail(TailBranch {
                        angle_from_parent: info.angle_from_parent,
                        length: info.length,
                        stem_radius: info.stem_radius,
                        sack_radius: info.sack_radius.unwrap(),
                    })
                };

                let idx = tree.items.len();
                tree.items.push(child_branch);
                BranchId(idx)
            };

            (make_child(left), make_child(right))
        }

        let mut tree = BranchTree {
            items: Vec::new(),
            // Set an initially invalid `BranchId`; we'll replace this later.
            root_id: BranchId(usize::MAX),
            root_start_pos: start.pos,
        };
        let (left, right) = full_gen(&mut tree, 1, start, gen);

        let root = Branch::Stem(StemBranch {
            // Becuase it's the first branch, its full angle needs to be represented by the angle
            // from its "parent" -- even though that doesn't exist.
            angle_from_parent: start.total_angle,
            length: start.length,
            stem_radius: start.stem_radius,
            left_child: left,
            right_child: right,
        });

        tree.root_id = BranchId(tree.items.len());
        tree.items.push(root);

        tree
    }
}

fn main() {
    let start = gen::standard_init_parent(0.3, 0.05);
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

    // let generator = EqualChildGenerator { max_depth: 6 };
    let generator = MirroredTerminalChildGenerator;
    let tree = BranchTree::from_generator(start, &generator);
    img_config
        .make_image(&tree)
        .save("branch-tree.png")
        .expect("failed to write the image");
}
