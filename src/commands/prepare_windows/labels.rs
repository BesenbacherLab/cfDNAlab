/// Label validation helpers and reserved words.

pub const RESERVED_LABEL_NONE: &str = "none";
pub const MISSING_GROUP_LABEL: &str = "[NA]";
pub const NO_NEAR_LABEL: &str = "[NONE]";
pub const NO_NEAR_BIN_LABEL: &str = "[NO-NEAR]";

const RESERVED_LABELS: [&str; 6] = [
    "input",
    "near-side",
    "near-name",
    "bin",
    "cluster",
    RESERVED_LABEL_NONE,
];

/// Return true if the label token is reserved.
pub fn is_reserved_label(token: &str) -> bool {
    RESERVED_LABELS
        .iter()
        .any(|reserved| token.eq_ignore_ascii_case(reserved))
}

/// Validate a label token for user-defined labels.
///
/// Tokens must be ASCII alphanumerics and cannot be reserved.
pub fn validate_label_token(token: &str, context: &str) -> Result<(), String> {
    if token.is_empty() {
        return Err(format!("{context} cannot be empty"));
    }
    if is_reserved_label(token) {
        return Err(format!("{context} cannot be '{token}'"));
    }
    if !token.chars().all(|ch| ch.is_ascii_alphanumeric()) {
        return Err(format!("{context} must be ASCII alphanumeric: '{token}'"));
    }
    Ok(())
}

use crate::commands::prepare_windows::config::ComposeSpec;
use anyhow::{Result as AnyResult, bail};
use fxhash::FxHashMap;
use std::collections::BTreeSet;

/// Atomic label parts carried by each window.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum AtomicLabelPart {
    Input,
    NearSide,
    NearName,
    Bin,
    Cluster,
}

impl AtomicLabelPart {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Input => "input",
            Self::NearSide => "near-side",
            Self::NearName => "near-name",
            Self::Bin => "bin",
            Self::Cluster => "cluster",
        }
    }

    pub fn from_token(token: &str) -> Option<Self> {
        match token {
            "input" => Some(Self::Input),
            "near-side" => Some(Self::NearSide),
            "near-name" => Some(Self::NearName),
            "bin" => Some(Self::Bin),
            "cluster" => Some(Self::Cluster),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum LabelKey {
    /// Reference an atomic label part like `input` or `near-side`.
    Atomic(AtomicLabelPart),
    /// Reference a named composition by index in the schema composition list.
    ///
    /// The index is the position of the composition in the `--compose` list
    /// after validation. The schema stores compositions in declaration order,
    /// so earlier `--compose` entries have lower indices.
    Composition(usize),
}

/// A composition part used when defining a named label.
///
/// Each part is either an atomic label or a reference to an earlier
/// named composition in the schema.
#[derive(Clone, Debug)]
pub enum LabelPartRef {
    Atomic(AtomicLabelPart),
    Composition(usize),
}

/// A named label composition built from ordered parts.
///
/// The `parts` are stored in the user-provided order so composition output
/// stays stable. The `depends_on` set tracks which atomic parts are required
/// for optimization decisions during rendering.
#[derive(Clone, Debug)]
pub struct LabelComposition {
    pub name: String,
    pub parts: Vec<LabelPartRef>,
    pub depends_on: BTreeSet<AtomicLabelPart>,
}

/// Resolved schema for named compositions and their indices.
///
/// This stores compositions in declaration order and provides a
/// name-to-index lookup for fast resolution.
#[derive(Clone, Debug, Default)]
pub struct LabelSchema {
    compositions: Vec<LabelComposition>,
    composition_by_name: FxHashMap<String, usize>,
}

impl LabelSchema {
    /// Build a validated label schema from `--compose` specifications.
    ///
    /// Compositions can reference atomic parts or earlier compositions.
    /// Unknown or cyclic references are rejected.
    ///
    /// Parameters
    /// ----------
    /// - `specs`:
    ///     Composition specs in the order provided by the user.
    ///
    /// Returns
    /// -------
    /// - `schema`:
    ///     Schema with resolved part references and dependency sets.
    pub fn new(specs: &[ComposeSpec]) -> AnyResult<Self> {
        let mut schema = Self::default();
        for spec in specs {
            if AtomicLabelPart::from_token(&spec.name).is_some() {
                bail!("compose name '{}' conflicts with atomic part", spec.name);
            }
            if schema.composition_by_name.contains_key(&spec.name) {
                bail!("compose name '{}' is defined more than once", spec.name);
            }

            let mut parts: Vec<LabelPartRef> = Vec::with_capacity(spec.parts.len());
            let mut depends_on: BTreeSet<AtomicLabelPart> = BTreeSet::new();

            for part in &spec.parts {
                if let Some(atomic) = AtomicLabelPart::from_token(part.as_str()) {
                    parts.push(LabelPartRef::Atomic(atomic));
                    depends_on.insert(atomic);
                    continue;
                }
                if let Some(idx) = schema.composition_by_name.get(part).copied() {
                    parts.push(LabelPartRef::Composition(idx));
                    depends_on.extend(schema.compositions[idx].depends_on.iter().copied());
                    continue;
                }
                bail!("compose '{}' references unknown part '{}'", spec.name, part);
            }

            let idx = schema.compositions.len();
            schema.compositions.push(LabelComposition {
                name: spec.name.clone(),
                parts,
                depends_on,
            });
            schema.composition_by_name.insert(spec.name.clone(), idx);
        }
        Ok(schema)
    }

    /// Resolve a label key token into an atomic or composition reference.
    ///
    /// Parameters
    /// ----------
    /// - `token`:
    ///     Label key token from CLI arguments.
    ///
    /// Returns
    /// -------
    /// - `key`:
    ///     Resolved label key.
    pub fn resolve_key(&self, token: &str) -> AnyResult<LabelKey> {
        if let Some(atomic) = AtomicLabelPart::from_token(token) {
            return Ok(LabelKey::Atomic(atomic));
        }
        if let Some(idx) = self.composition_by_name.get(token).copied() {
            return Ok(LabelKey::Composition(idx));
        }
        let mut known: Vec<String> = vec![
            AtomicLabelPart::Input.as_str().to_string(),
            AtomicLabelPart::NearSide.as_str().to_string(),
            AtomicLabelPart::NearName.as_str().to_string(),
            AtomicLabelPart::Bin.as_str().to_string(),
            AtomicLabelPart::Cluster.as_str().to_string(),
        ];
        let mut comp_names: Vec<String> =
            self.compositions.iter().map(|c| c.name.clone()).collect();
        comp_names.sort();
        known.extend(comp_names);
        bail!(
            "unknown label key '{}'; expected one of {}",
            token,
            known.join(", ")
        );
    }

    /// Resolve an ordered list of label keys from CLI tokens.
    ///
    /// Parameters
    /// ----------
    /// - `tokens`:
    ///     List of label key strings.
    ///
    /// Returns
    /// -------
    /// - `keys`:
    ///     Resolved label keys in the same order.
    pub fn resolve_keys(&self, tokens: &[String]) -> AnyResult<Vec<LabelKey>> {
        let mut keys = Vec::with_capacity(tokens.len());
        for token in tokens {
            keys.push(self.resolve_key(token)?);
        }
        Ok(keys)
    }

    /// Return composition definitions in declaration order.
    ///
    /// The order matches the user input and is stable for rendering.
    ///
    /// Returns
    /// -------
    /// - `compositions`:
    ///     Slice of composition definitions.
    pub fn compositions(&self) -> &[LabelComposition] {
        &self.compositions
    }
}

/// Atomic label tuple attached to a window.
#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct LabelTuple {
    pub input: String,
    pub near_side: Option<String>,
    pub near_name: Option<String>,
    pub bin: Option<String>,
    pub cluster: Option<String>,
}

impl LabelTuple {
    pub fn new(input: String) -> Self {
        Self {
            input,
            near_side: None,
            near_name: None,
            bin: None,
            cluster: None,
        }
    }
}

#[inline]
fn atomic_value<'a>(tuple: &'a LabelTuple, part: AtomicLabelPart) -> &'a str {
    match part {
        AtomicLabelPart::Input => tuple.input.as_str(),
        AtomicLabelPart::NearSide => tuple.near_side.as_deref().unwrap_or(""),
        AtomicLabelPart::NearName => tuple.near_name.as_deref().unwrap_or(""),
        AtomicLabelPart::Bin => tuple.bin.as_deref().unwrap_or(""),
        AtomicLabelPart::Cluster => tuple.cluster.as_deref().unwrap_or(""),
    }
}

fn build_composition_values(tuple: &LabelTuple, compositions: &[LabelComposition]) -> Vec<String> {
    let mut values: Vec<String> = Vec::with_capacity(compositions.len());
    for composition in compositions {
        let mut parts: Vec<String> = Vec::with_capacity(composition.parts.len());
        for part in &composition.parts {
            match part {
                LabelPartRef::Atomic(atomic) => {
                    parts.push(atomic_value(tuple, *atomic).to_string());
                }
                LabelPartRef::Composition(idx) => {
                    parts.push(values[*idx].clone());
                }
            }
        }
        values.push(parts.join("."));
    }
    values
}

/// Build composition values for each tuple in order.
///
/// Produces one vector per tuple, where each entry matches the composition
/// index in the schema.
///
/// Parameters
/// ----------
/// - `tuples`:
///     Label tuples to expand.
/// - `schema`:
///     Resolved composition schema.
///
/// Returns
/// -------
/// - `values`:
///     Per-tuple composition values in schema order.
#[inline]
pub fn build_tuple_compositions(tuples: &[LabelTuple], schema: &LabelSchema) -> Vec<Vec<String>> {
    tuples
        .iter()
        .map(|tuple| build_composition_values(tuple, schema.compositions()))
        .collect()
}

/// Sort and deduplicate label tuples for stable output.
///
/// This keeps tuple order deterministic so comma-separated lists preserve
/// pairings across label columns.
///
/// Parameters
/// ----------
/// - `tuples`:
///     Label tuples to normalize.
pub fn normalize_label_tuples(tuples: &mut Vec<LabelTuple>) {
    tuples.sort();
    tuples.dedup();
}

#[inline]
fn all_parts_match(tuples: &[LabelTuple], parts: &BTreeSet<AtomicLabelPart>) -> bool {
    for part in parts {
        if matches!(part, AtomicLabelPart::Input) {
            continue;
        }
        let mut it = tuples.iter();
        let Some(first) = it.next() else {
            return true;
        };
        let first_value = atomic_value(first, *part);
        if it.any(|tuple| atomic_value(tuple, *part) != first_value) {
            return false;
        }
    }
    true
}

#[inline]
fn combined_input_value(tuples: &[LabelTuple]) -> String {
    let mut values: Vec<&str> = tuples.iter().map(|tuple| tuple.input.as_str()).collect();
    values.sort();
    values.dedup();
    values.join("__")
}

fn render_atomic_label(tuples: &[LabelTuple], part: AtomicLabelPart) -> String {
    if tuples.is_empty() {
        return String::new();
    }

    if matches!(part, AtomicLabelPart::Input) {
        let mut it = tuples.iter();
        let first = it.next().unwrap();
        let first_value = first.input.as_str();
        if it.all(|tuple| tuple.input.as_str() == first_value) {
            return first_value.to_string();
        }
        let non_input_parts: BTreeSet<AtomicLabelPart> = [
            AtomicLabelPart::NearSide,
            AtomicLabelPart::NearName,
            AtomicLabelPart::Bin,
            AtomicLabelPart::Cluster,
        ]
        .iter()
        .copied()
        .collect();
        let only_input_differs = all_parts_match(tuples, &non_input_parts);
        if only_input_differs {
            return combined_input_value(tuples);
        }
        return tuples
            .iter()
            .map(|tuple| tuple.input.as_str())
            .collect::<Vec<&str>>()
            .join(",");
    }

    let mut it = tuples.iter();
    let first = it.next().unwrap();
    let first_value = atomic_value(first, part);
    if it.all(|tuple| atomic_value(tuple, part) == first_value) {
        return first_value.to_string();
    }
    tuples
        .iter()
        .map(|tuple| atomic_value(tuple, part))
        .collect::<Vec<&str>>()
        .join(",")
}

fn render_composition_label(
    tuples: &[LabelTuple],
    tuple_compositions: &[Vec<String>],
    composition_idx: usize,
    schema: &LabelSchema,
) -> String {
    if tuples.is_empty() {
        return String::new();
    }

    let mut values: Vec<&str> = Vec::with_capacity(tuples.len());
    for (tuple_idx, _tuple) in tuples.iter().enumerate() {
        values.push(tuple_compositions[tuple_idx][composition_idx].as_str());
    }

    if values.iter().all(|value| *value == values[0]) {
        return values[0].to_string();
    }

    let composition = &schema.compositions()[composition_idx];
    let depends_on_input = composition.depends_on.contains(&AtomicLabelPart::Input);
    if depends_on_input && all_parts_match(tuples, &composition.depends_on) {
        let combined_input = combined_input_value(tuples);
        let mut synthetic = tuples[0].clone();
        synthetic.input = combined_input;
        let synthetic_values = build_composition_values(&synthetic, schema.compositions());
        return synthetic_values[composition_idx].clone();
    }

    values.join(",")
}

/// Render a label value for a window based on its tuples.
///
/// Uses compact forms when only the input value varies, otherwise it emits
/// comma-separated lists that preserve tuple order.
///
/// Parameters
/// ----------
/// - `tuples`:
///     Label tuples for the window.
/// - `tuple_compositions`:
///     Precomputed composition values for those tuples.
/// - `key`:
///     Atomic or composition key to render.
/// - `schema`:
///     Resolved composition schema.
///
/// Returns
/// -------
/// - `label`:
///     Rendered label string for output or grouping.
pub fn render_label_for_key(
    tuples: &[LabelTuple],
    tuple_compositions: &[Vec<String>],
    key: &LabelKey,
    schema: &LabelSchema,
) -> String {
    match key {
        LabelKey::Atomic(part) => render_atomic_label(tuples, *part),
        LabelKey::Composition(idx) => {
            render_composition_label(tuples, tuple_compositions, *idx, schema)
        }
    }
}
