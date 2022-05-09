//! Airflow simulator for my 3rd-year computer science project
//!
//! The main entrypoint is actually in [`cli::run`] ('src/cli.rs'), which handles the various
//! things that we may want to do -- that in turn calls the `run` method on [`AppSettings`]

use serde::Deserialize;
use std::borrow::Cow;
use std::fmt::{self, Display, Formatter, Write as _};
use std::fs::File;
use std::io::{self, Write as _};
use std::iter;
use std::path::Path;
use std::process::exit;
use std::time::{Duration, Instant};

mod branch_id;
mod cli;
mod float;
mod gen;
mod img;
mod point;
mod sim;
mod tree;
mod vec_utils;

pub use branch_id::BranchId;

use float::Float;
use gen::{EqualChildGenerator, FromJsonGenerator};
use img::{rgb, ImageConfig, PixelCount};
use point::Point;
use sim::SimulationEnvironment;
use tree::BranchTree;

struct AppSettings<'cli> {
    total_time: Float,
    timestep: Float,
    display_method: cli::DisplayMethod<'cli>,
    model: cli::Model,
}

/// The "direction" of a child node -- it's either the left or right child
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum ChildDirection {
    Left,
    Right,
}

/// A path through the tree to a node (subtree)
pub type TreePath = Vec<ChildDirection>;

impl Display for ChildDirection {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            ChildDirection::Left => f.write_str(".left"),
            ChildDirection::Right => f.write_str(".right"),
        }
    }
}

struct DisplayTreePath<'p>(&'p [ChildDirection]);

impl Display for DisplayTreePath<'_> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        for c in self.0 {
            c.fmt(f)?;
        }
        Ok(())
    }
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

    /// The radius of the tube, in m
    radius: Float,
    /// The length of the tube, in m
    length: Float,
    /// The velocity at which air is flowing into the tube, away from the trachea (positive if
    /// breathing in)
    ///
    /// Units of m³/s
    flow_in: Float,
    /// The pressure at the distal end of the tube (i.e. away from the trachea), in Pascals
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

/// We'll say that the default functional residual capacity is about 2.5 liters. It's not
/// *necessarily* accurate but close enough so that we're talking about measurements that probably
/// correspond to a real person.
const DEFAULT_LUNG_FRC: Float = 0.0025;

/// Length of the trachea, in meters. Average length is around 12cm, which is 0.12m
const TRACHEA_LENGTH: Float = 0.12;

/// Radius of the trachea, in meters. Seems like average diameter is around 1.5-2cm. Picked 2cm
/// because it's a simpler number & this doesn't need to be perfectly accurate yet.
const TRACHEA_RADIUS: Float = 0.01;

#[derive(Debug, Copy, Clone, Default, Deserialize)]
#[serde(default)]
pub struct EnvConfig {
    pleural_pressure: PleuralPressureConfig,
    /// The initial total volume of the lungs, in m³
    initial_volume: Option<Float>,
}

/// Values used for the pleural pressure, relative to atmospheric pressure
///
/// These values are later converted when passed to the simulator itself.
#[derive(Debug, Copy, Clone, Deserialize)]
#[serde(default)]
struct PleuralPressureConfig {
    init: Float,
    lo: Float,
    hi: Float,
    period: Float,
    /// Maximum pleural pressure during normal breathing; used for calculating compliance for
    /// desired FRC
    normal_hi: Float,
}

impl Default for PleuralPressureConfig {
    fn default() -> Self {
        PleuralPressureConfig {
            init: -1000.0,
            lo: -2000.0,
            hi: -1000.0,
            period: 4.0,
            normal_hi: -1000.0,
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct Schedule {
    interpolate: Interpolate,
    keyframes: Vec<KeyFrame>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Deserialize)]
enum Interpolate {
    /// Linear interpolation, i.e. `f(x) = x`
    #[serde(alias = "linear")]
    Linear,
    /// Smoothed interpolation via tanh: `f(x) = .5 * tanh(6x - 3)/tanh(3) + .5
    ///
    /// The correction factor `1/tanh(3)` is there to avoid a sudden jump from 0.997 to 1, which
    /// occurs without it
    ///
    /// This is the most curved of the available interpolation methods, about twice as much as
    /// `Trig` and `Cubic`.
    #[serde(alias = "tanh")]
    Tanh,
    /// Smoothed interpolation via cosine: `f(x) = (1 - cos(pi*x)) / 2`
    ///
    /// This is very similar in shape to `Cubic`, and only ever so slightly more curved. It is
    /// about half as curved as `Tanh`.
    #[serde(alias = "trig", alias = "sin", alias = "cos")]
    Trig,
    /// Smoothed interpolation with a cubic function: `f(x) = 3 x^2 - 2 x^3`
    ///
    /// This is very similar in shape to `Trig`, and only ever so slightly less curved. It is about
    /// half as curved as `Tanh`.
    #[serde(alias = "cubic")]
    Cubic,
}

#[derive(Debug, Copy, Clone, Deserialize)]
struct KeyFrame {
    time: Float,
    is_degraded: bool,
}

// Default schedule is just degraded
impl Default for Schedule {
    fn default() -> Self {
        Schedule {
            interpolate: Interpolate::Linear,
            keyframes: vec![KeyFrame {
                time: 0.0,
                is_degraded: true,
            }],
        }
    }
}

impl Schedule {
    fn degradation_at_time(&self, time: Float) -> Float {
        let frame_to_factor = |f: &KeyFrame| f.is_degraded as u8 as Float;

        match self
            .keyframes
            .binary_search_by(|f| f.time.partial_cmp(&time).unwrap())
        {
            // Simple case that occurs very rarely.
            Ok(i) => frame_to_factor(&self.keyframes[i]),

            // Before the first keyframe; return it:
            Err(0) => frame_to_factor(&self.keyframes[0]),
            // After last keyframe; continue with it:
            Err(i) if i == self.keyframes.len() => frame_to_factor(self.keyframes.last().unwrap()),
            // Between two keyframes, mix the factors:
            Err(i) => {
                let before = &self.keyframes[i - 1];
                let after = &self.keyframes[i];

                let f = (time - before.time) / (after.time - before.time);
                Self::interpolate(
                    self.interpolate,
                    f,
                    frame_to_factor(before),
                    frame_to_factor(after),
                )
            }
        }
    }

    fn interpolate(int: Interpolate, f: Float, before: Float, after: Float) -> Float {
        if before == after {
            return before;
        }

        match int {
            Interpolate::Linear => f * after + (1.0 - f) * before,
            Interpolate::Tanh => {
                let factor = 1.0 / Float::tanh(3.0);
                let lerp = 0.5 * factor * Float::tanh(6.0 * f - 3.0) + 0.5;
                Self::interpolate(Interpolate::Linear, lerp, before, after)
            }
            Interpolate::Trig => {
                let lerp = (1.0 - Float::cos(float::PI * f)) / 2.0;
                Self::interpolate(Interpolate::Linear, lerp, before, after)
            }
            Interpolate::Cubic => {
                let lerp = 3.0 * f.powi(2) - 2.0 * f.powi(3);
                Self::interpolate(Interpolate::Linear, lerp, before, after)
            }
        }
    }
}

struct TreePair {
    nominal: BranchTree,
    degraded: BranchTree,
}

fn main() {
    // Internally calls `AppSettings::run`
    cli::run()
}

// Callback to output the state of the tree, given the current time, state of the tree, degradation
// factor, and time spent on the simulation tick
type DisplayCallback<'a> =
    Box<dyn 'a + FnMut(Float, &BranchTree, &SimulationEnvironment, Float, Duration)>;

impl AppSettings<'_> {
    /// Runs the app until completion, using the settings filled by the `cli` module
    fn run(&self) {
        let (tree_pair, bounds, lung_frc, env_cfg, sched) =
            self.make_tree_and_bounds().unwrap_or_else(|e| {
                eprintln!("{:?}", e.wrap_err("failed to construct model"));
                exit(1)
            });

        let mut sim_env = SimulationEnvironment {
            // Shift the pleural pressure here so that it's no lnger relative to atmospheric
            // pressure
            pleural_pressure: ATMOSPHERIC_PRESSURE + env_cfg.pleural_pressure.at_time(0.0),
            tracheal_pressure: ATMOSPHERIC_PRESSURE,
            air_viscosity: 1.8e-5,
            air_density: 1.225,
        };

        let nominal_state = tree_pair.nominal.structure_state();
        let degraded_state = tree_pair.degraded.structure_state();

        let mut tree = tree_pair.degraded;
        let mut deg_factor = sched.degradation_at_time(0.0);
        tree.set_structure(&lerp(&nominal_state, &degraded_state, deg_factor));

        let initial_volume = env_cfg.initial_volume.or(lung_frc);

        sim_env.reset(&mut tree, initial_volume);

        let mut callback = match &self.display_method {
            cli::DisplayMethod::Csv { file, add_paths } => self.csv_callback(*file, &add_paths),
            cli::DisplayMethod::Png { file_pattern } => self.png_callback(bounds, file_pattern),
        };

        callback(0.0, &tree, &sim_env, deg_factor, Duration::ZERO);

        for i in 1.. {
            // Better to recompute than += timestep because the latter is less precise.
            let current_time = i as Float * self.timestep;

            if self.total_time != 0.0 && i as Float * self.timestep > self.total_time {
                break;
            }

            // The pleural pressure needs to update at each timestep. Same shift by atmospheric
            // pressure as above so that it's absolute pressure not relative
            sim_env.pleural_pressure = ATMOSPHERIC_PRESSURE + env_cfg.pleural_pressure.at_time(current_time);

            // And possibly we need to degrade (or un-degrade) the tree
            let old_deg_factor = deg_factor;
            deg_factor = sched.degradation_at_time(current_time);
            if deg_factor != old_deg_factor {
                tree.set_structure(&lerp(&nominal_state, &degraded_state, deg_factor));
            }

            let start = Instant::now();

            if let Err(msg) = sim_env.do_tick(&mut tree, self.timestep) {
                println!("\nfailed to do simulation tick: {}", msg);
                return;
            }

            let elapsed = start.elapsed();

            callback(current_time, &tree, &sim_env, deg_factor, elapsed);
        }
    }

    /// Creates the `BranchTree` and returns a few relevant items:
    ///
    /// * The nominal and degraded `BranchTree`s;
    /// * The upper-right corner of the rectangle of the minimum view into the tree(s);
    /// * Configuration of the simulation environment;
    /// * The expected functional residual capacity of the lungs, if provided; and
    /// * A schedule for determining the degradation amount at each point in time
    ///
    /// It is assumed that the tree shouldn't display further down or left than (0,0).
    fn make_tree_and_bounds(
        &self,
    ) -> eyre::Result<(TreePair, Point, Option<Float>, EnvConfig, Schedule)> {
        match &self.model {
            cli::Model::FromJson { file } => Self::make_json_tree(&file),
            cli::Model::Symmetric { depth } => Ok(Self::make_symmetric_tree(*depth)),
        }
    }

    fn make_json_tree(
        file: &Path,
    ) -> eyre::Result<(TreePair, Point, Option<Float>, EnvConfig, Schedule)> {
        let gen = FromJsonGenerator::from_file(file)?;
        let start = gen.start_parent_info();
        let pair = TreePair {
            nominal: BranchTree::from_generator(start, &gen, false),
            degraded: BranchTree::from_generator(start, &gen, true),
        };
        let lung_frc = Some(gen.lung_frc());
        Ok((
            pair,
            gen.upper_right(),
            lung_frc,
            gen.env_config(),
            gen.schedule(),
        ))
    }

    fn make_symmetric_tree(depth: usize) -> (TreePair, Point, Option<Float>, EnvConfig, Schedule) {
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

        let env_config = EnvConfig::default();

        // Calculate the compliance so that it produces the "correct" total lung volume at
        // atmospheric pressure
        let compliance = {
            let n_acinar = 1 << (depth + 1);
            let volume_per = DEFAULT_LUNG_FRC / n_acinar as Float;

            // volume = compliance * (pressure - pleural pressure), therefore:
            //
            //   compliance = volume / (pressure - pleural pressure)
            //
            // FRC occurs when pleural pressure is at maximum (for normal breaths), so we use that
            // here.
            volume_per / -env_config.pleural_pressure.normal_hi
        };

        println!("info: calculated compliance as {}", compliance);

        let generator = EqualChildGenerator {
            compliance,
            // +1 here because EqualChildGenerator counts depth from 1, and we want to count from
            // zero.
            max_depth: depth + 1,
        };
        let tree = BranchTree::from_generator(start, &generator, true);

        // Realistically, this doesn't *quite* need to be a square, but the actual bounding width
        // is a bit more complicated and I don't really want to deal with that right now.
        let upper_right = Point {
            x: TOTAL_WIDTH,
            y: TOTAL_HEIGHT,
        };

        // No difference.
        let pair = TreePair {
            nominal: tree.clone(),
            degraded: tree,
        };

        (pair, upper_right, None, env_config, Schedule::default())
    }

    fn csv_callback<'p>(&self, file: Option<&str>, paths: &'p [TreePath]) -> DisplayCallback<'p> {
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
        let mut header =
            "time,tick duration,pleural pressure,degradation ratio,flow out,total volume"
                .to_owned();
        for p in paths {
            let s = p
                .iter()
                .map(|d| match d {
                    ChildDirection::Left => ".left",
                    ChildDirection::Right => ".right",
                })
                .collect::<String>();

            header.push(',');
            header.push_str(&s);
            header.push_str(" flow out");
            header.push(',');
            header.push_str(&s);
            header.push_str(" total volume");
        }

        writeln!(writer, "{}", header).expect("failed to write CSV header");

        Box::new(
            move |time: Float,
                  tree: &BranchTree,
                  sim_env: &SimulationEnvironment,
                  deg: Float,
                  dur: Duration| {
                let root_path = Vec::new();

                let mut flow_and_volume_values = Vec::new();

                for p in iter::once(&root_path).chain(paths) {
                    // Get the root node for this path
                    let mut path_root = &tree[tree.root_id()];
                    for (i, c) in p.iter().enumerate() {
                        match path_root {
                            Branch::Bifurcation(b) => match c {
                                ChildDirection::Left => path_root = &tree[b.left_child],
                                ChildDirection::Right => path_root = &tree[b.right_child],
                            },
                            Branch::Acinar(_) => {
                                eprintln!(
                                    "bad path '{}', ends early at '{}'",
                                    DisplayTreePath(&p),
                                    DisplayTreePath(&p[..i + 1])
                                );
                                exit(1);
                            }
                        }
                    }

                    let mut flow_out = -path_root.tube().flow_in;
                    // Convert flow m³/s -> L/s
                    flow_out *= 1e3;
                    flow_and_volume_values.push(flow_out);

                    // Find the total volume:
                    let mut stack = vec![path_root];
                    let mut volume = 0.0; // Units of m³/s
                    while let Some(n) = stack.pop() {
                        match n {
                            Branch::Acinar(a) => volume += a.volume,
                            Branch::Bifurcation(b) => {
                                stack.push(&tree[b.right_child]);
                                stack.push(&tree[b.left_child]);
                            }
                        }
                    }

                    // Convert volume m³ -> L
                    volume *= 1e3;
                    flow_and_volume_values.push(volume);
                }

                // We ensure there's at least two decimal places because pasting into Google Sheets
                // will sometimes get messed up if there's no decimal place in the number.
                let mut line = format!(
                    "{:.2},{:.3},{:.2},{:.2}",
                    time,
                    dur.as_secs_f64() * 1000.0,
                    sim_env.pleural_pressure,
                    deg
                );
                for v in flow_and_volume_values {
                    line.push(',');
                    write!(line, "{:.4}", v).unwrap();
                }

                writeln!(writer, "{}", line).expect("failed to write CSV line");
                writer.flush().expect("Failed to flush CSV line");
            },
        )
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
            background: rgb(0x000000),
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

        Box::new(
            move |_time: Float,
                  tree: &BranchTree,
                  _: &SimulationEnvironment,
                  _: Float,
                  _: Duration| {
                img_config
                    .make_image(tree)
                    .save(cli::substitute_png_file_pattern(&pat, n, max_imgs))
                    .expect("failed to write an image");

                n += 1;
            },
        )
    }
}

fn lerp<'a>(xs: &'a [Float], ys: &'a [Float], factor: Float) -> Cow<'a, [Float]> {
    if factor == 0.0 {
        return Cow::Borrowed(xs);
    } else if factor == 1.0 {
        return Cow::Borrowed(ys);
    }

    let x_f = 1.0 - factor;
    let y_f = factor;
    Cow::Owned(xs.iter().zip(ys).map(|(x, y)| x * x_f + y * y_f).collect())
}

impl PleuralPressureConfig {
    /// Returns the pleural pressure (i.e. the "pressure" exerted on the lungs by the diaphragm) at
    /// the given point in time
    ///
    /// The pressure is determined by a sinusoidal function from `self.lo` to `self.`hi` with the
    /// period given by `self.period`. It's not the *most* accurate representation in the world,
    /// but it works to demonstrate that the model is performing sensibly.
    ///
    /// We start with a pressure of `self.init` at time 0s, and then increase from there (i.e.
    /// breathe out).
    fn at_time(&self, time: Float) -> Float {
        assert!((self.lo..=self.hi).contains(&self.init));

        let amplitude = (self.hi - self.lo) / 2.0;
        let center = (self.hi + self.lo) / 2.0;

        let start_phase_height = (center - self.init) / amplitude;
        let shift = start_phase_height.acos();

        let stretch = float::PI * 2.0 / self.period;

        amplitude + self.lo - amplitude * (time * stretch + shift).cos()
    }
}
