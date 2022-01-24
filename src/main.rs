//! Airflow simulator for my 3rd-year computer science project
//!
//! The main entrypoint is actually in [`cli::run`] ('src/cli.rs'), which handles the various
//! things that we may want to do -- that in turn calls the `run` method on [`AppSettings`]

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::process::exit;

mod branch_id;
mod cli;
mod float;
mod gen;
mod img;
mod point;
mod sim;
mod tree;

pub use branch_id::BranchId;

use float::Float;
use gen::{EqualChildGenerator, FromJsonGenerator};
use img::{rgb, rgba, ImageConfig, PixelCount};
use point::Point;
use sim::SimulationEnvironment;
use tree::BranchTree;

struct AppSettings<'cli> {
    total_time: Float,
    timestep: Float,
    display_method: cli::DisplayMethod<'cli>,
    model: cli::Model,
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

/// We'll say that the total volume of the lungs is about 2 liters. It's not incredibly accurate,
/// but close enough so that we're talking about measurements that probably correspond to a real
/// person
const TOTAL_LUNG_VOLUME: Float = 0.002;

/// Length of the trachea, in meters. Average length is around 12cm
const TRACHEA_LENGTH: Float = 0.12;

/// Radius of the trachea, in meters. Seems like average diameter is around 1.5-2cm. Picked 2cm
/// because it's a simpler number & this doesn't need to be perfectly accurate yet.
const TRACHEA_RADIUS: Float = 0.01;

fn main() {
    // Internally calls `AppSettings::run`
    cli::run()
}

type DisplayCallback = Box<dyn FnMut(&BranchTree, &SimulationEnvironment)>;

impl AppSettings<'_> {
    /// Runs the app until completion, using the settings filled by the `cli` module
    fn run(&self) {
        let (mut tree, bounds) = self.make_tree_and_bounds().unwrap_or_else(|e| {
            eprintln!("{:?}", e.wrap_err("failed to construct model"));
            exit(1)
        });

        let mut sim_env = SimulationEnvironment {
            pleural_pressure: Self::pleural_pressure_at_time(0.0),
            tracheal_pressure: ATMOSPHERIC_PRESSURE,
            air_viscosity: 1.8e-5,
            air_density: 1.225,
        };

        sim_env.reset(&mut tree);

        let mut callback = match &self.display_method {
            cli::DisplayMethod::Csv { file } => self.csv_callback(*file),
            cli::DisplayMethod::Png { file_pattern } => self.png_callback(bounds, file_pattern),
        };

        callback(&tree, &sim_env);

        for i in 1.. {
            // Better to recompute than += timestep because the latter is less precise.
            let current_time = i as Float * self.timestep;

            if self.total_time != 0.0 && i as Float * self.timestep > self.total_time {
                break;
            }

            // The pleural pressure needs to update at each timestep
            sim_env.pleural_pressure = Self::pleural_pressure_at_time(current_time);

            if let Err(msg) = sim_env.do_tick(&mut tree, self.timestep) {
                println!("\nfailed to do simulation tick: {}", msg);
                return;
            }

            callback(&tree, &sim_env);
        }
    }

    /// Creates the `BranchTree` and returns the upper-right corner of the rectangle of the minimum
    /// view into it
    ///
    /// It is assumed that the tree shouldn't display further down or left than (0,0).
    fn make_tree_and_bounds(&self) -> eyre::Result<(BranchTree, Point)> {
        match &self.model {
            cli::Model::FromJson { file } => Self::make_json_tree(&file),
            cli::Model::Symmetric { depth } => Ok(Self::make_symmetric_tree(*depth)),
        }
    }

    fn make_json_tree(file: &Path) -> eyre::Result<(BranchTree, Point)> {
        let gen = FromJsonGenerator::from_file(file)?;
        let tree = BranchTree::from_generator(gen.start_parent_info(), &gen);
        Ok((tree, gen.upper_right()))
    }

    fn make_symmetric_tree(depth: usize) -> (BranchTree, Point) {
        const TOTAL_HEIGHT: Float = TRACHEA_LENGTH * 2.7;
        const TOTAL_WIDTH: Float = TOTAL_HEIGHT * 1.3;

        let start = gen::ParentInfo {
            id: 0,
            pos: Point {
                x: TOTAL_WIDTH / 2.0,
                y: TOTAL_HEIGHT,
            },
            // Angle points downwards -- minus pi/2
            total_angle: -float::FRAC_PI_2,
            length: TRACHEA_LENGTH,
            tube_radius: TRACHEA_RADIUS,
        };

        // Calculate the compliance so that it produces the "correct" total lung volume at
        // atmospheric pressure
        let compliance = {
            let n_acinar = 1 << (depth + 1);
            let volume_per = TOTAL_LUNG_VOLUME / n_acinar as Float;

            // volume = compliance * (pressure - pleural pressure), therefore:
            //
            //   compliance = volume / (pressure - pleural pressure)
            //
            // because we're assuming pleural pressure is zero, we just get:
            volume_per / ATMOSPHERIC_PRESSURE
        };

        println!("info: calculated compliance as {}", compliance);

        let generator = EqualChildGenerator {
            // 0.2L/cmH₂O gives 2.0394e-6, in m³/Pa. For some reason, that *really* doesn't play
            // nice here, so we're using a value that does. <- TODO/FIXME
            compliance,
            // +1 here because EqualChildGenerator counts depth from 1, and we want to count from
            // zero.
            max_depth: depth + 1,
        };
        let tree = BranchTree::from_generator(start, &generator);

        // Realistically, this doesn't *quite* need to be a square, but the actual bounding width
        // is a bit more complicated and I don't really want to deal with that right now.
        let upper_right = Point {
            x: TOTAL_WIDTH,
            y: TOTAL_HEIGHT,
        };

        (tree, upper_right)
    }

    fn csv_callback(&self, file: Option<&str>) -> DisplayCallback {
        // Get the file or handle to stdout:
        let mut writer: Box<dyn io::Write> = match file {
            Some(f) => match File::create(f) {
                Ok(w) => Box::new(w),
                Err(e) => {
                    eprintln!("Failed to open file {:?} for writing: {}", f, e);
                    exit(1)
                }
            },
            None => Box::new(io::stdout()),
        };

        // Off the bat, print out the csv header information we need:
        writeln!(writer, "time,pleural pressure,flow out,total volume")
            .expect("failed to write CSV header");

        let mut i = 0_usize;
        let timestep = self.timestep; // borrow checker gets mad if this isn't out here :(
        Box::new(move |tree: &BranchTree, sim_env: &SimulationEnvironment| {
            let time = i as Float * timestep;

            // Flow out, in mm³/s:
            let flow_out = -tree.root_flow_in() * 1e9;

            // Total volume, in m³/s
            let mut total_volume = 0.0;
            for i in 0..tree.count_acinar_regions() {
                match &tree[BranchId::new(BranchKind::Acinar, i)] {
                    Branch::Acinar(a) => total_volume += a.volume,
                    _ => unreachable!(),
                }
            }

            // And now in mm³/s
            total_volume *= 1e9;

            writeln!(
                writer,
                "{},{},{},{}",
                time, sim_env.pleural_pressure, flow_out, total_volume
            )
            .expect("failed to write CSV line");
            writer.flush().expect("Failed to flush CSV line");

            i += 1;
        })
    }

    fn png_callback(&self, bounds: Point, file_pattern: &str) -> DisplayCallback {
        const MIN_DIM: PixelCount = 500;
        const PAD: PixelCount = 50;

        let scale;

        if bounds.x > bounds.y {
            scale = MIN_DIM as Float / bounds.x;
        } else {
            scale = MIN_DIM as Float / bounds.y;
        }

        let width = (scale * bounds.x) as PixelCount + PAD;
        let height = (scale * bounds.y) as PixelCount + PAD;

        let img_config = ImageConfig {
            centered_at: Point {
                x: bounds.x / 2.0,
                y: bounds.y / 2.0,
            },
            width,
            height,
            scale,
            grid_lines: vec![
                // X & Y axes:
                (rgb(0x000000), img::Axis::X, 0.0, 5),
                (rgb(0x000000), img::Axis::Y, 0.0, 5),
            ],
            background: rgba(0x00000000),
            stem_color: rgb(0xFFFFFF),
            sack_color: rgb(0xFF0000),
        };

        // Image number:
        let mut n = 0;
        let max_imgs = if self.total_time == 0.0 {
            // Using '1' here means that we won't get any padding on the image names. That's ok
            1
        } else {
            // +1 for the initial image at time 0
            (self.total_time / self.timestep).ceil() as u64 + 1
        };

        // Need to clone `file_pattern` out here so that the closure can be 'static.
        let pat = file_pattern.to_owned();

        Box::new(move |tree: &BranchTree, _: &SimulationEnvironment| {
            img_config
                .make_image(tree)
                .save(cli::substitute_png_file_pattern(&pat, n, max_imgs))
                .expect("failed to write an image");

            n += 1;
        })
    }

    /// Returns the pleural pressure (i.e. the "pressure" exerted on the lungs by the diaphragm) at
    /// the given point in time
    ///
    /// The pressure is determined by a sinusoidal function from -500 to 1500 with a period of 4
    /// seconds. It's not the *most* accurate representation in the world, but it works to
    /// demonstrate that the model is performing sensibly.
    ///
    /// We set the pleural pressure equal to zero at time 0s, then increase from there (i.e.
    /// breathe out).
    fn pleural_pressure_at_time(time: Float) -> Float {
        const PERIOD: Float = 4.0;
        const START_PRESSURE: Float = 0.0;
        const LO_PRESSURE: Float = -500.0;
        const HI_PRESSURE: Float = 1500.0;

        debug_assert!((LO_PRESSURE..HI_PRESSURE).contains(&START_PRESSURE));

        // These should all be optimized into constants as well; it's just more annoying to write
        // them as consts.
        let amplitude = (HI_PRESSURE - LO_PRESSURE) / 2.0;
        let center = (HI_PRESSURE + LO_PRESSURE) / 2.0;

        let start_phase_height = (center - START_PRESSURE) / amplitude;
        let shift = start_phase_height.acos();

        let stretch = float::PI * 2.0 / PERIOD;

        amplitude + LO_PRESSURE - amplitude * (time * stretch + shift).cos()
    }
}
