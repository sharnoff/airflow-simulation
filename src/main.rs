//! Airflow simulator for my 3rd-year computer science project

use std::io::{self, Write};

mod branch_id;
mod float;
mod gen;
mod img;
mod sim;
mod tree;

pub use branch_id::BranchId;

use float::Float;
use gen::basic::MirroredTerminalChildGenerator;
use img::{rgb, rgba, ImageConfig};
use tree::BranchTree;

/// A physical point in 2D space, used to represent one end of a branch
///
/// The values are in units of meters.
///
/// We treat positive X as to the right and positive Y as up, with both values greater than zero.
/// This distinction matters when we're rendering an image of the tree.
///
/// It's also worth noting that x and y values are typically (though not *necessarily*) greater
/// than zero.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Point {
    x: Float,
    y: Float,
}

/// Marker for the type of a branch. Primarily used with [`BranchId`] methods.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum BranchKind {
    Bifurcation = 0,
    Acinar = 1,
}

/// An individual node in the tree of branches
///
/// Branches are tyically referred to by their [`BranchId`]. Both `Bifurcation` and `Acinar` share
/// a `tube` field representing their bronchiole, which can be retrieved with the [`tube`] method.
///
/// [`stem`]: Self::stem
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Branch {
    Bifurcation(Bifurcation),
    Acinar(AcinarRegion),
}

/// The default value for `Tube.flow_rate` as we're generating a [`BranchTree`]
const UNSET_FLOW_RATE: Float = 0.0;

/// The default value for `Tube.end_pressure` as we're generating a [`BranchTree`]
const UNSET_END_PRESSURE: Float = 0.0;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Tube {
    /// The ange of the tube, relative to its parent. Measured in radians, always within 0..2π
    angle_from_parent: Float,

    /// The radius of the tube, in mm
    radius: Float,
    /// The length of the tube, in mm
    length: Float,
    /// The velocity at which air is flowing into the tube, away from the trachea (positive if
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
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Bifurcation {
    /// Information about the bronchus/bronchiole going into this bifurcation
    tube: Tube,

    left_child: BranchId,
    right_child: BranchId,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AcinarRegion {
    /// Information about the final bronchiole leading into this acinar region
    tube: Tube,

    /// The "compliance" of the region -- essentially how much it is willing to stretch
    ///
    /// Essentially volume / pressure; units of m³/Pa = m⁴s²/kg
    compliance: Float,
    /// Total volume of the region; units of m³
    volume: Float,
}

impl Branch {
    /// Returns the kind of branch represented
    fn kind(&self) -> BranchKind {
        match self {
            Branch::Bifurcation(_) => BranchKind::Bifurcation,
            Branch::Acinar(_) => BranchKind::Acinar,
        }
    }

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

/// Atmospheric pressure at sea level, in Pascals
const ATMOSPHERIC_PRESSURE: Float = 101_325.0;

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

    let generator = MirroredTerminalChildGenerator {
        // 0.2L/cmH₂O gives the value below, in m³/Pa:
        compliance: 2.0394e-6,
    };
    let mut tree = BranchTree::from_generator(start, &generator);
    let mut sim_env = sim::SimulationEnvironment {
        pleural_pressure: 500.0,
        tracheal_pressure: ATMOSPHERIC_PRESSURE,
        air_viscosity: 1.8e-5,
        air_density: 1.225,
    };

    sim_env.reset(&mut tree);

    // Make the image for the initial state
    img_config
        .make_image(&tree)
        .save("branch-tree-00.png")
        .expect("failed to write the image");

    const N_IMGS: usize = 100;
    const TIMESTEP: Float = 0.01;

    // Then, make all of the images that we need. We'll display the current "time" that we're
    // evaluating in seconds below.
    for i in 1..N_IMGS {
        print!("\rtime: {:.2}s...", (i - 1) as Float * TIMESTEP);
        io::stdout().flush().expect("failed to flush stdout");

        // increase the pleural pressure; i.e. breathe out
        //
        // TODO: This should actually follow some kind of sinusoidal function. We'll do that later.
        sim_env.pleural_pressure += 5.0;

        if let Err(msg) = sim_env.do_tick(&mut tree, TIMESTEP) {
            println!("\nfailed to do simulation tick: {}", msg);
            return;
        }

        img_config
            .make_image(&tree)
            .save(format!("branch-tree-{:02}.png", i))
            .expect("failed to write the image");
    }
}
