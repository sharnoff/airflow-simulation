//! Wrapper module around the `Point` type

use crate::Float;
use std::ops::*;

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
    pub x: Float,
    pub y: Float,
}

impl Add<Point> for Point {
    type Output = Self;

    fn add(self, other: Point) -> Self {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub<Point> for Point {
    type Output = Self;

    fn sub(self, other: Point) -> Self {
        self + -1.0 * other
    }
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self {
        self * -1.0
    }
}

impl Mul<Point> for Float {
    type Output = Point;

    fn mul(self, point: Point) -> Point {
        point * self
    }
}

impl Mul<Float> for Point {
    type Output = Self;

    fn mul(self, scale: Float) -> Self {
        Point {
            x: scale * self.x,
            y: scale * self.y,
        }
    }
}

impl MulAssign<Float> for Point {
    fn mul_assign(&mut self, scale: Float) {
        *self = *self * scale;
    }
}
