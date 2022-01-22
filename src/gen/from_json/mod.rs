//! The "JSON" model generator, instantiating the tree based on a variety of structured generation
//! options

use super::*;
use crate::{float, ATMOSPHERIC_PRESSURE, TOTAL_LUNG_VOLUME, TRACHEA_LENGTH, TRACHEA_RADIUS};
use eyre::{eyre, Context};
use std::fs;
use std::path::Path;

mod config;
mod ctx;

pub use config::ParsedConfig;
use ctx::{BranchSize, GenCtx};

/// The maximum allowed depth. This is a practical limit so that certain calculations can be made
/// easy -- it shouldn't cause issues in typical use.
const MAX_DEPTH: usize = 64;

/// A minimum depth of 2 is required so that we can have at least one branch
///
/// This requirement may be lifted at some point, if the interface for `BranchGenerator` changes.
const REQUIRED_MIN_DEPTH: usize = 2;

/// The fraction to pad by
const SIZE_PADDING: Float = 0.3;

/// A tree generator that uses a configuration parsed from a chosen JSON file
#[derive(Debug)]
pub struct FromJsonGenerator {
    by_parent_id: Vec<Branch>,
    upper_right: Point,
    start_parent_info: ParentInfo,
}

#[derive(Debug)]
struct Branch {
    left: ChildInfo,
    right: ChildInfo,
}

// The implementation of `BranchGenerator` just uses the values that we've previously constructed
// from a dry-run.
impl BranchGenerator for FromJsonGenerator {
    fn make_children(&self, parent: ParentInfo, _depth: usize) -> (ChildInfo, ChildInfo) {
        let branch = &self.by_parent_id[parent.id];
        (branch.left, branch.right)
    }
}

impl FromJsonGenerator {
    /// Produces a new `FromJsonGenerator` from the JSON at the given path
    pub fn from_file(file: &Path) -> eyre::Result<Self> {
        let file_content = fs::read_to_string(file)
            .wrap_err_with(|| format!("failed to read file at {:?}", file.to_string_lossy()))?;

        let parsed =
            serde_json::from_str(&file_content).wrap_err("could not deserialize JSON structure")?;

        Self::from_parsed(parsed)
    }

    /// Returns the starting `ParentInfo` that should be used with this generator
    pub fn start_parent_info(&self) -> ParentInfo {
        self.start_parent_info
    }

    /// Returns a sensible bounding upper-right corner for displaying the result of generating this
    /// tree
    pub fn upper_right(&self) -> Point {
        self.upper_right
    }

    // Internally:
    //
    // This function needs to generate the entire tree as a first-pass in order to determine some
    // of the scaling information (e.g. distribution of compliance, center point, etc.)
    fn from_parsed(parsed: ParsedConfig) -> eyre::Result<Self> {
        let total_volume = parsed.config.total_volume.unwrap_or(TOTAL_LUNG_VOLUME);
        let trachea_length = parsed.config.trachea_length.unwrap_or(TRACHEA_LENGTH);
        let trachea_radius = parsed.config.trachea_radius.unwrap_or(TRACHEA_RADIUS);

        // Check the maximum depth:
        if parsed.config.common.max_depth.get() < REQUIRED_MIN_DEPTH {
            return Err(eyre!("maximum depth must be at least 2"))
                .context("invalid value at .config.max_depth in JSON model spec");
        }

        // We'll now construct the full tree.
        let mut by_parent_id = Vec::new();

        let initial_branch = BranchSize {
            length: trachea_length,
            nominal_radius: trachea_radius,
            abnormal_radius: trachea_radius,
        };

        let root_start_pos = Point { x: 0.0, y: 0.0 };
        let initial_ctx = GenCtx::root(
            &parsed.root,
            &parsed.config.common,
            initial_branch,
            root_start_pos,
        );

        let mut stack = vec![(initial_ctx, true)];
        let mut child_stack = Vec::new();

        let mut most_negative_point = root_start_pos;
        let mut most_positive_point = root_start_pos;

        while let Some((ctx, first_pass)) = stack.pop() {
            let this_branch = ctx.branch_size();
            let angle_from_parent = ctx.angle_from_parent();

            // If this branch is on its first pass (i.e. first time popped off the stack), then we
            // should process it by (eventually) adding to the child stack
            if first_pass {
                // If this is an acinar region, then there's not much more we can do:
                if ctx.is_acinar() {
                    let (nominal_c, abnormal_c) = ctx.compliance();
                    let base_compliance =
                        Self::guess_compliance_for_depth(total_volume, ctx.depth());

                    // Update extreme points:
                    let end_pos = ctx.end_pos();
                    let radius_guess = Self::guess_radius_for_compliance(base_compliance);

                    most_positive_point.x = most_positive_point.x.max(end_pos.x + radius_guess);
                    most_positive_point.y = most_positive_point.y.max(end_pos.y + radius_guess);
                    most_negative_point.x = most_negative_point.x.min(end_pos.x - radius_guess);
                    most_negative_point.y = most_negative_point.y.min(end_pos.y - radius_guess);

                    let nominal_compliance = nominal_c * base_compliance;
                    let info = ChildInfo {
                        angle_from_parent,
                        length: this_branch.length,
                        tube_radius: this_branch.abnormal_radius,
                        compliance: Some(abnormal_c * base_compliance),
                    };

                    child_stack.push((info, nominal_compliance));
                    continue;
                }

                // Otherwise, this must be a non-terminal branch (i.e. with children)
                let (left_ctx, right_ctx) = ctx.children();

                stack.push((ctx, false));
                // We push left before right because we want to generate the right-hand side first.
                // The reason for that is that our `by_parent_id` array is actually constructed in
                // reverse order -- the last branch added will be the root. So we have to
                // additionally reverse the left-then-right ordering that's normally present.
                stack.push((left_ctx, true));
                stack.push((right_ctx, true));
            } else {
                // Non-terminal branch that's already made its children; we can grab them from the
                // child stack.
                //
                // Note that the the order is the *same* as when we pushed to 'stack'; this is
                // because "pushed to stack first" -> "handled first" -> "lower in 'child_stack'";
                // contrary to what you might otherwise assume.
                let (left, l_compliance) = child_stack.pop().unwrap();
                let (right, r_compliance) = child_stack.pop().unwrap();

                by_parent_id.push(Branch { left, right });

                let this_info = ChildInfo {
                    angle_from_parent,
                    length: this_branch.length,
                    tube_radius: this_branch.abnormal_radius,
                    // non-terminal branches are marked as such by not having a compliance.
                    compliance: None,
                };

                child_stack.push((this_info, l_compliance + r_compliance));
            }
        }

        // There's a note inside the loop above about constructing this in reverse (by the calls to
        // `stack.push(...)`. Refer there for insight into why we're now reversing this.
        by_parent_id.reverse();

        // At the end, there should be exactly the root node left:
        assert!(child_stack.len() == 1);

        // After generating the tree, we now need to do a full pass to rescale compliances to match
        // the expected total volume.
        let generated_total_compliance = child_stack[0].1;
        let expected_total_compliance = Self::guess_compliance_for_depth(total_volume, 1);

        // If the generated compliance is not within a small amount of the actual compliance, then
        // we should rescale:
        let compliance_rescale = expected_total_compliance / generated_total_compliance;
        let epsilon = 1e-6;
        if (1.0 - compliance_rescale).abs() > epsilon {
            by_parent_id
                .iter_mut()
                .map(|b| [&mut b.left, &mut b.right])
                .flatten()
                .filter_map(|info| info.compliance.as_mut())
                .for_each(|compliance| *compliance *= compliance_rescale);
        }

        // And now, we need to adjust the position of the root so that the bottom-most and
        // left-most points in the tree align with y = 0 and x = 0, respectively.
        //
        // This *could* be designed around elsewhere, but it's just easier to do it here.

        // We want a little padding to account for acinar regions (because those aren't counted in
        // the "most positive/negative" points:
        most_negative_point *= 1.0 + SIZE_PADDING;
        most_positive_point *= 1.0 + SIZE_PADDING;

        let start_parent_info = ParentInfo {
            id: 0,
            total_angle: -float::FRAC_PI_2, // point downwards
            pos: -most_negative_point,
            length: initial_branch.length,
            tube_radius: initial_branch.abnormal_radius,
        };

        let upper_right = most_positive_point - most_negative_point;

        Ok(FromJsonGenerator {
            by_parent_id,
            upper_right,
            start_parent_info,
        })
    }

    /// Returns a guess for the compliance of an acinar region at the given depth
    ///
    /// If there are no modified compliances, this will be correct. Otherwise, it may require some
    /// scaling after the initial generation.
    fn guess_compliance_for_depth(total_volume: Float, depth: usize) -> Float {
        // NOTE: This requires depth <= 64. This is accounted for by the `MAX_DEPTH` constant at the
        // top of this module.
        let n_acinar = 1 << (depth as u64 - 1);

        let volume_per = total_volume / n_acinar as Float;

        volume_per / ATMOSPHERIC_PRESSURE
    }

    fn guess_radius_for_compliance(compliance: Float) -> Float {
        // Just use atmospheric pressure
        let volume = compliance * ATMOSPHERIC_PRESSURE;
        (volume * 3.0 / (4.0 * float::PI)).powf(1.0 / 3.0)
    }
}
