//! Wrapper module for the config structure, rooted at [`ParsedConfig`]

use crate::float::Float;
use serde::Deserialize;
use std::num::NonZeroUsize;

#[derive(Debug, Deserialize)]
pub struct ParsedConfig {
    pub root: GeneratorNode,
    pub config: RootGenerationConfig,
}

#[derive(Debug, Deserialize)]
pub struct RootGenerationConfig {
    #[serde(flatten)]
    pub common: BaseGenerationConfig,

    /// The total volume of all acinar regions in the lungs, in m³
    ///
    /// If present, overrides the default given by [`TOTAL_LUNG_VOLUME`]
    pub total_volume: Option<Float>,

    /// The length of the trachea, in meters
    ///
    /// If present, overrides the default given by [`TRACHEA_LENGTH`]
    pub trachea_length: Option<Float>,

    /// The radius of the trachea, in meters
    ///
    /// If present, overrides the default given by [`TRACHEA_RADIUS`]
    pub trachea_radius: Option<Float>,
}

/// The fallback configuration for the parameters in each [`GenerationConfig`]
#[derive(Debug, Deserialize)]
pub struct BaseGenerationConfig {
    /// The ratio of decrease in length from a branch to its children
    pub branch_length_decrease: Float,
    /// The ratio of decrease in radius from a branch to its children
    pub branch_radius_decrease: Float,
    /// The angle (in 0..π) that child branches are rotated away from the parent by default
    pub split_angle: Float,
    /// The maximum depth of child branches to generate, with a minimum of 1
    pub max_depth: NonZeroUsize,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
pub enum GeneratorNode {
    #[serde(rename = "manual")]
    Manual {
        left: Option<Box<GeneratorNode>>,
        right: Option<Box<GeneratorNode>>,
        left_override: Option<ManualSettings>,
        right_override: Option<ManualSettings>,
    },
    #[serde(rename = "auto")]
    Auto(Option<GenerationConfig>),
}

/// Overriden settings for *just one* node - i.e. not any of its children/descendants
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ManualSettings {
    /// If present, overrides the angle of the child, relative to its parent
    ///
    /// The angle is automatically sign-adjusted for left vs. right children, so only positive
    /// angles less than π should be supplied.
    pub relative_angle: Option<Float>,

    /// If present, overrides the ratio from the parent's branch length to this one
    pub relative_length: Option<Float>,

    /// If present, overrides the "nominal" ratio from the parent's branch radius to this one
    #[serde(alias = "relative_radius")]
    pub relative_radius_nominal: Option<Float>,

    /// If present, overrides the "abnormal" ratio from the parent's branch radius to this one
    pub relative_radius_abnormal: Option<Float>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenerationConfig {
    /// Sets or overrides the decrease factor for child branches
    pub branch_length_decrease: Option<(Option<Float>, Option<Float>, ChildInheritConfig)>,

    /// Sets or overrides the distribution of radius sizes for child branches
    pub branch_radius_decrease: Option<(Option<Float>, Option<Float>, ChildInheritConfig)>,

    /// Sets or overrides the pair of angles to offset child branches by
    ///
    /// Like in `ManualSettings`, both of these values should be positive.
    pub child_angles: Option<(Float, Float, ChildInheritConfig)>,

    /// Sets a maximum *additional* depth to generate all descendants to
    ///
    /// Note: The depth is relative to *this* node, i.e. it is not aware of its own position of the
    /// tree.
    pub max_depth: Option<usize>,

    /// Increases or decreases the compliance of descendants, relative to the parent configuration
    pub relative_compliance: Option<ComplianceScaling>,
}

#[derive(Debug, Deserialize)]
pub struct ComplianceScaling {
    /// "abnormal" compliance meaning that compliance is distributed according to the volume of the
    /// lungs *as if this were not present*.
    ///
    /// Generally this is used to distinguish between lung structure vs. temporary changes in lung
    /// morphology.
    pub abnormal: Option<Float>,

    /// Opposite of "abnormal" -- scales the compliance that would be expcted under typical lung
    /// function
    pub nominal: Option<Float>,
}

#[derive(Debug, Default, Deserialize)]
pub struct ChildInheritConfig {
    /// If present, marks the setting as temporary -- i.e. it should reset to the next applicable
    /// value once the count is provided
    ///
    /// A value of zero would mean that the setting is never applied; as such, it is disallowed.
    pub reset_after: Option<NonZeroUsize>,
}
