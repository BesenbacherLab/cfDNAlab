#[cfg(feature = "cli")]
use clap::ValueEnum;
#[cfg(feature = "cmd_fragment_kmers")]
use serde::{Deserialize, Serialize};
use strum_macros::{EnumCount as EnumCountMacro, EnumIter};

/// Group positional selections by which fragment-side coordinate system they use.
///
/// Left and right selections are anchored to the fragment ends, while `Mid`
/// selections are centered around the midpoint frame. This enum is shared by
/// the positional selection code, visualization helpers, and the k-mer codec
/// orientation logic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[cfg(feature = "cmd_fragment_kmers")]
pub enum PositionGroup {
    Left,
    Right,
    Mid,
}

/// Whether coordinates should be interpreted as forward or reverse-oriented.
///
/// Left and midpoint selections are naturally forward-oriented, while right-end
/// selections are reverse-oriented because they are interpreted from the right
/// fragment end inward.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg(feature = "cmd_fragment_kmers")]
pub enum PositionOrientation {
    Forward,
    Reverse,
}

#[cfg(feature = "cmd_fragment_kmers")]
impl PositionOrientation {
    /// Convert a position group into its natural orientation.
    ///
    /// Parameters
    /// ----------
    /// - `group`:
    ///   Positional group being interpreted
    ///
    /// Returns
    /// -------
    /// - `PositionOrientation`:
    ///   Forward for left and midpoint groups, reverse for right-end groups
    pub fn from_position_group(group: PositionGroup) -> PositionOrientation {
        match group {
            PositionGroup::Left | PositionGroup::Mid => PositionOrientation::Forward,
            PositionGroup::Right => PositionOrientation::Reverse,
        }
    }
}

/// Describe which coordinate frame a positional selection should use.
///
/// These frames are shared by `fragment-kmers` and the visualization helper so
/// the user can prototype the exact same positional logic that the counting
/// command later uses on real fragments.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash, EnumCountMacro, EnumIter)]
pub enum ReferenceFrame {
    #[default]
    Left,
    Right,
    PerEnd,
    Nearest,
    Mid,
}

impl ReferenceFrame {
    /// Convert the frame into its stable CLI/config string.
    ///
    /// Parameters
    /// ----------
    /// - `self`:
    ///   Frame value to stringify
    ///
    /// Returns
    /// -------
    /// - `&'static str`:
    ///   Stable lower-case frame name
    pub fn as_str(self) -> &'static str {
        match self {
            ReferenceFrame::Left => "left",
            ReferenceFrame::Right => "right",
            ReferenceFrame::PerEnd => "per-end",
            ReferenceFrame::Nearest => "nearest",
            ReferenceFrame::Mid => "mid",
        }
    }
}

/// Decide whether positions should come from reads or the reference span.
///
/// This is shared between `fragment-kmers` and the visualization helper so the
/// same coordinate choice can be described and previewed with a single enum.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BasesFrom {
    /// Always use reference positions regardless of read coverage.
    #[default]
    Reference,
    /// Prefer observed read coordinates, but fall back to the reference span.
    PreferReads,
    /// Only include positions covered by either read.
    Reads,
    /// Clamp to the read nearest to the frame origin.
    NearestRead,
}

impl BasesFrom {
    /// Convert the source mode into its stable CLI/config string.
    ///
    /// Parameters
    /// ----------
    /// - `self`:
    ///   Source mode to stringify
    ///
    /// Returns
    /// -------
    /// - `&'static str`:
    ///   Stable lower-case mode name
    pub fn as_str(self) -> &'static str {
        match self {
            BasesFrom::Reference => "reference",
            BasesFrom::PreferReads => "prefer-reads",
            BasesFrom::Reads => "reads",
            BasesFrom::NearestRead => "nearest-read",
        }
    }
}

/// Choose how overlapping read mismatches should be resolved.
///
/// When paired reads disagree about a base in an overlapping region, this enum
/// records which source should win so positional extraction stays deterministic.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MismatchBasesFrom {
    #[default]
    NearestRead,
    BaseQuality,
    Reference,
}

impl MismatchBasesFrom {
    /// Convert the mismatch mode into its stable CLI/config string.
    ///
    /// Parameters
    /// ----------
    /// - `self`:
    ///   Mismatch-resolution mode to stringify
    ///
    /// Returns
    /// -------
    /// - `&'static str`:
    ///   Stable lower-case mode name
    pub fn as_str(self) -> &'static str {
        match self {
            MismatchBasesFrom::NearestRead => "nearest-read",
            MismatchBasesFrom::BaseQuality => "base-quality",
            MismatchBasesFrom::Reference => "reference",
        }
    }
}
