//! Wrapper module for [`BranchId`]

use crate::BranchKind;
use std::fmt::{self, Debug, Formatter};

/// Unique identifier for a [`Branch`]
///
/// The type of the branch (either [`Bifurcation`] or [`AcinarRegion`]) is encoded in the least
/// significant bit of the value, which means that the `i`th `Bifurcation` is actually given the id
/// `2*i` and the `i`th `AcinarRegion` is `2*i + 1`.
///
/// For this reason, the API of `BranchId`s is restricted to keep this invariant held correctly.
///
/// One side effect is that we only allow up to `isize::MAX` of either branch type, instead of
/// `usize::MAX` of the two combined. This really shouldn't matter in practice.
///
/// [`Branch`]: super::Branch
/// [`Bifurcation`]: super::Bifurcation
/// [`AcinarRegion`]: super::AcinarRegion
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BranchId(usize);

impl BranchId {
    /// Creates a new `BranchId` for the type of branch with given index
    ///
    /// These values can be recovered via calls to the [`deconstruct`] method.
    ///
    /// [`deconstruct`]: Self::deconstruct
    pub(crate) fn new(kind: BranchKind, idx: usize) -> Self {
        if idx > isize::MAX as usize {
            panic!("BranchId cannot be constructed with index greater than isize::MAX")
        }

        BranchId(idx << 1 | kind as usize)
    }

    /// Returns the `BranchId` as the pair of "branch kind, index" that it represents
    pub(crate) fn deconstruct(&self) -> (BranchKind, usize) {
        let kind = match self.0 % 2 == 0 {
            true => BranchKind::Bifurcation,
            false => BranchKind::Acinar,
        };

        (kind, self.0 >> 1)
    }
}

impl Debug for BranchId {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self.deconstruct() {
            (BranchKind::Bifurcation, idx) => write!(f, "Bifurcation#{}", idx),
            (BranchKind::Acinar, idx) => write!(f, "Acinar#{}", idx),
        }
    }
}
