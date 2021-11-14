//! "Pretty" image generation to display branch trees

use crate::{float, Branch, BranchId, BranchTree, Float, Point};

use image::{ImageBuffer, Rgba};
use imageproc::drawing::{self, Blend};

/// Creates a color from the provided RGBA integer
///
/// See also: [`rgb`].
///
/// ## Examples
///
/// ```
/// let transparent = rgba(0x00000000);
/// let red_tint = rgba(0xff000077);
/// let teal = rgba(0x00ffffff);
/// ```
///
/// ## Panics
///
/// This function does not panic; every `u32` is a valid RGBA color.
pub fn rgba(int: u32) -> Color {
    let r = ((int >> 24) & 0xff) as u8;
    let g = ((int >> 16) & 0xff) as u8;
    let b = ((int >> 8) & 0xff) as u8;
    let a = ((int >> 0) & 0xff) as u8;

    Rgba([r, g, b, a])
}

/// Creates a color from the provided RGB integer
///
/// See also: [`rgba`].
///
/// ## Examples
///
/// ```
/// let red = rgb(0xff0000);
/// let teal = rgb(0x00ffff);
/// let orange = rgb(0xffaa00);
/// ```
///
/// ## Panics
///
/// This function panics if the provided integer is greater than `0xffffff`.
pub fn rgb(int: u32) -> Color {
    assert!(int <= 0xffffff);

    rgba((int << 8) + 0xff)
}

/// The color type we're using
pub type Color = Rgba<u8>;

/// Type alias to represent a number of pixels. This is only provided so that the meaning behind
/// types can be more clear.
pub type PixelCount = u32;

/// Configuration items for making an image
#[derive(Debug)]
pub struct ImageConfig {
    /// The center point of the region to use for the image
    pub centered_at: Point,

    /// The width of the produced image
    pub width: PixelCount,
    /// The height of the produced image
    pub height: PixelCount,

    /// The scale at which to produce the image
    ///
    /// This value is equivalent to the pixel size that a region size of 1 unit is converted to.
    /// Essentially, if the image has dimensions 1000x500, setting `scale = 500` will make the
    /// image focus on a region that's 2.0x1.0
    pub scale: Float,

    /// Grid lines, given by their color, value, and pixel width on a particular axis
    ///
    /// Lines that are provided earlier will be given precedence (i.e. drawn on top of other lines)
    pub grid_lines: Vec<(Color, Axis, Float, PixelCount)>,

    /// Background color of the image
    pub background: Color,
    /// Color of stems
    pub stem_color: Color,
    /// Color of individual sacks
    pub sack_color: Color,
}

/// The axis on which to draw grid lines -- either X or Y
#[derive(Copy, Clone, Debug)]
pub enum Axis {
    X,
    Y,
}

/// Helper struct to record position of the endpoints of a branch's stem
#[derive(Copy, Clone, Debug)]
struct StemPos {
    /// The position of the end of the stem that connects to this branch's parent
    start: Point,
    /// The position other end
    end: Point,
    /// The angle of the stem - i.e. the rotation upwards from the X axis that would be required to
    /// move `start` and `end` from the X axis into their current position
    angle: Float,
}

#[rustfmt::skip]
impl StemPos {
    /// Gives the top-left corner of the rotated rectangle representing a stem at this position and
    /// angle with a radius of `r`
    fn top_left(&self, r: Float) -> Point {
        Point { x: self.start.x + r * self.angle.sin(), y: self.start.y - r * self.angle.cos() }
    }

    /// Gives the top-right corner of the rotated rectangle representing a stem at this position
    /// and angle with a radius of `r`
    fn top_right(&self, r: Float) -> Point {
        Point { x: self.start.x - r * self.angle.sin(), y: self.start.y + r * self.angle.cos() }
    }

    /// Gives the bottom-left corner of the rotated rectangle representing a stem at this position
    /// and angle with a radius of `r`
    fn bot_left(&self, r: Float) -> Point {
        Point { x: self.end.x + r * self.angle.sin(), y: self.end.y - r * self.angle.cos() }
    }

    /// Gives the bottom-right corner of the rotated rectangle representing a stem at this position
    /// and angle with a radius of `r`
    fn bot_right(&self, r: Float) -> Point {
        Point { x: self.end.x - r * self.angle.sin(), y: self.end.y + r * self.angle.cos() }
    }
}

/// Helper type alias
type ImageCanvas = Blend<ImageBuffer<Color, Vec<u8>>>;

impl ImageConfig {
    /// Creates an `ImageBuffer` representing the provided tree, using the available configuration
    /// options
    pub fn make_image(&self, tree: &BranchTree) -> ImageBuffer<Color, Vec<u8>> {
        // Blending from the `image` crate generally assumes that `self` is the background and
        // `other` is the foreground. This basicaly means that we can get proper relative depths
        // (mostly) by just drawing things in order.
        //
        // The implementation of `Canvas` for `Blend` will appropriately use `image`'s blending.
        let mut buf = Blend(ImageBuffer::from_pixel(
            self.width,
            self.height,
            self.background,
        ));

        let ctx = DrawContext {
            bot_left: Point {
                x: self.centered_at.x - (self.width as Float / self.scale / 2.0),
                y: self.centered_at.y - (self.height as Float / self.scale / 2.0),
            },
            scale: self.scale,
            height: self.height,
        };

        self.draw_axes(&mut buf, ctx);
        self.draw_full_tree(&mut buf, tree, ctx);

        buf.0
    }

    /// Draws the required axes onto the image
    fn draw_axes(&self, canvas: &mut ImageCanvas, ctx: DrawContext) {
        // Convenience function to process a pair into a Point
        let map = |x, y| imageproc::point::Point { x, y };

        // The iterator here is reversed so that earlier grid lines get drawn last (i.e. on top)
        for &(color, axis, value, px_width) in self.grid_lines.iter().rev() {
            let px_radius = (px_width / 2 + (if px_width % 2 == 1 { 1 } else { 0 })) as i32;

            let points = match axis {
                Axis::X => {
                    let (img_coord, _) = ctx.point_to_coords(Point { x: value, y: 0.0 });
                    if !(0..self.width as i32).contains(&img_coord) {
                        continue;
                    }

                    [
                        map(img_coord - px_radius, 0),
                        map(img_coord + px_radius, 0),
                        map(img_coord + px_radius, self.height as i32),
                        map(img_coord - px_radius, self.height as i32),
                    ]
                }
                Axis::Y => {
                    let (_, img_coord) = ctx.point_to_coords(Point { x: 0.0, y: value });
                    if !(0..self.height as i32).contains(&img_coord) {
                        continue;
                    }

                    [
                        map(0, img_coord - px_radius),
                        map(0, img_coord + px_radius),
                        map(self.width as i32, img_coord + px_radius),
                        map(self.width as i32, img_coord - px_radius),
                    ]
                }
            };

            drawing::draw_polygon_mut(canvas, &points, color);
        }
    }

    /// Draws the full tree to the canvas, using the existing configuration
    ///
    /// In order to get proper relative depths, we have to draw from the bottom of the tree
    /// upwards, so we have a pre-processing step where we sort every branch (including terminals)
    /// into buckets based on their depth.
    fn draw_full_tree(&self, canvas: &mut ImageCanvas, tree: &BranchTree, ctx: DrawContext) {
        type DepthBuckets<'tree> = Vec<Vec<(StemPos, &'tree Branch)>>;

        let mut depth_buckets: DepthBuckets = Vec::new();

        // Recursive helper function to add a branch to its bucket
        fn add_to_bucket<'tree>(
            buckets: &mut DepthBuckets<'tree>,
            id: BranchId,
            stem_start: Point,
            parent_angle: Float,
            tree: &'tree BranchTree,
            depth: usize,
        ) {
            let branch = &tree[id];

            // Calculate this entry:
            let angle = parent_angle + branch.tube().angle_from_parent;
            let length = branch.tube().length;

            let stem_end = Point {
                x: stem_start.x + length * angle.cos(),
                y: stem_start.y + length * angle.sin(),
            };

            // And then add it to its bucket:
            let stem = StemPos {
                start: stem_start,
                end: stem_end,
                angle,
            };

            // Becuase this was called by its parent (or the initial call), buckets.len() will
            // always be <= depth
            if buckets.len() == depth {
                buckets.push(Vec::new());
            }
            buckets[depth].push((stem, branch));

            // Add the children if there are any
            if let Branch::Bifurcation(b) = branch {
                // We use `stem_end` as the child's `stem_start` because that's how they're linked
                // together.
                add_to_bucket(buckets, b.left_child, stem_end, angle, tree, depth + 1);
                add_to_bucket(buckets, b.right_child, stem_end, angle, tree, depth + 1);
            }
        }

        add_to_bucket(
            &mut depth_buckets,
            tree.root_id(),
            tree.root_start_pos(),
            0.0,
            tree,
            0,
        );

        // After we've added everything to the buckets, we draw in reverse-depth order
        for bucket in depth_buckets.into_iter().rev() {
            for (stem, branch) in bucket {
                self.draw_branch(canvas, &ctx, stem, branch);
            }
        }
    }

    /// Draws the branch itself, without referring to sub-branches
    fn draw_branch(
        &self,
        canvas: &mut ImageCanvas,
        ctx: &DrawContext,
        stem: StemPos,
        branch: &Branch,
    ) {
        // Draw the stem.
        let radius = branch.tube().radius;

        // We can't actually draw the stem as a wide line - that isn't provided so instead we'll do
        // it with a rotated rectangle. This gets drawn as a polygon, so there unfortunately isn't
        // any antialiasing that we can do here.

        let map = |p| {
            let (x, y) = ctx.point_to_coords(p);
            imageproc::point::Point { x, y }
        };

        // Keep the points as an array, so that we don't have the overhead from vector allocations.
        //
        // Realistically, it doesn't slow it down *too* much, but it's enough of a difference that
        // it's worth doing.
        let poly_points = [
            map(stem.top_left(radius)),
            map(stem.top_right(radius)),
            map(stem.bot_right(radius)),
            map(stem.bot_left(radius)),
        ];

        drawing::draw_polygon_mut(canvas, &poly_points, self.stem_color);

        // If this is a terminal branch, we also need to draw the sack. That one's pretty easy -
        // it's just a circle.
        if let Branch::Acinar(b) = branch {
            // Volume of a sphere: V = 4/3 π r³, so:
            // r = (V 3/4π)^(1/3)
            let radius = (b.volume * 3.0 / (4.0 * float::PI)).powf(1.0 / 3.0);

            // Coordinates of the center of the sack
            let center_coords = ctx.point_to_coords(Point {
                x: stem.end.x + radius * stem.angle.cos(),
                y: stem.end.y + radius * stem.angle.sin(),
            });

            let radius_px = (radius * self.scale).round() as i32;

            drawing::draw_filled_circle_mut(canvas, center_coords, radius_px, self.sack_color);
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct DrawContext {
    /// The point corresponding to the bottom-left corner of the image
    bot_left: Point,

    /// The amount scale is multiplied, from the "true" coordinates to the image itself - as a
    /// float
    scale: Float,

    /// The height of the image. We need this because drawing has the origin at the top-left
    /// corner, so we need to flip the image to get it at the bottom-left.
    height: PixelCount,
}

impl DrawContext {
    /// Converts a point to its corresponding location in the image
    ///
    /// The values returned are signed because it's possible for points outside the image to still
    /// provide value (e.g. as the vertices in a polygon).
    ///
    /// Another piece is of note: when drawing, the origin is at the top-left corner. We want it at
    /// the bottom-left, so we have to flip the entire image.
    fn point_to_coords(&self, p: Point) -> (i32, i32) {
        let x = ((p.x - self.bot_left.x) * self.scale).round() as i32;
        let y = ((p.y - self.bot_left.y) * self.scale).round() as i32;
        (x, self.height as i32 - y)
    }
}
