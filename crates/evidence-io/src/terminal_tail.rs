//! Sparse, molecule-attached terminal poly(A)-tail evidence.
//!
//! This optional section losslessly retains every explicitly uniquely mapped (`NH=1`),
//! deduplicated cleavage-anchor event selected by the versioned extraction rule below. It is
//! deliberately not a sequence archive: clipped bases, qualities, and exact signal counts are
//! discarded after a bounded signal summary is retained.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};

use crate::archive::{put_svarint, put_varint};
use crate::format::Cursor;

pub const TERMINAL_TAIL_LAYOUT: &str = "chunk-event-list-v1";
pub const TERMINAL_TAIL_RULE: &str = "forward-cdna-terminal-softclip-v1";
pub const TERMINAL_TAIL_INDEX_SECTION: &str = "index.tail";
pub const TERMINAL_TAIL_INDEX_MAGIC: &[u8; 8] = b"TAILIDX1";
pub const TERMINAL_TAIL_CHUNK_MAGIC: &[u8; 8] = b"TAILCHN1";

/// Frozen selection thresholds for [`TERMINAL_TAIL_RULE`].
pub const MIN_CLIP: u8 = 6;
pub const MIN_TAIL_NUMERATOR: u8 = 4;
pub const MIN_TAIL_DENOMINATOR: u8 = 5;
pub const MIN_TERMINAL_RUN: u8 = 4;

/// Capability declaration stored in `meta.terminal_tail`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TerminalTailMetadata {
    pub layout: String,
    pub extraction_rule: String,
    pub alignment_scope: String,
    pub solo_strand: String,
    pub forward_edge: String,
    pub reverse_edge: String,
    pub coordinate: String,
    pub min_clip: u8,
    pub min_tail_numerator: u8,
    pub min_tail_denominator: u8,
    pub min_terminal_run: u8,
    pub selected_molecules: u64,
    pub events: u64,
    pub chunks: u32,
    pub sequence_retained: bool,
    pub quality_retained: bool,
}

impl TerminalTailMetadata {
    pub fn new(selected_molecules: u64, events: u64, chunks: u32) -> Self {
        Self {
            layout: TERMINAL_TAIL_LAYOUT.into(),
            extraction_rule: TERMINAL_TAIL_RULE.into(),
            alignment_scope: "mapped-primary-nonsupplementary-explicit-nh1".into(),
            solo_strand: "forward".into(),
            forward_edge: "trailing-soft-clip-A".into(),
            reverse_edge: "leading-soft-clip-T".into(),
            coordinate: "0-based-exclusive-end-plus;0-based-inclusive-start-minus".into(),
            min_clip: MIN_CLIP,
            min_tail_numerator: MIN_TAIL_NUMERATOR,
            min_tail_denominator: MIN_TAIL_DENOMINATOR,
            min_terminal_run: MIN_TERMINAL_RUN,
            selected_molecules,
            events,
            chunks,
            sequence_retained: false,
            quality_retained: false,
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.layout != TERMINAL_TAIL_LAYOUT {
            bail!("unsupported terminal-tail layout {}", self.layout);
        }
        if self.extraction_rule != TERMINAL_TAIL_RULE
            || self.alignment_scope != "mapped-primary-nonsupplementary-explicit-nh1"
            || self.solo_strand != "forward"
            || self.forward_edge != "trailing-soft-clip-A"
            || self.reverse_edge != "leading-soft-clip-T"
            || self.coordinate != "0-based-exclusive-end-plus;0-based-inclusive-start-minus"
            || self.min_clip != MIN_CLIP
            || self.min_tail_numerator != MIN_TAIL_NUMERATOR
            || self.min_tail_denominator != MIN_TAIL_DENOMINATOR
            || self.min_terminal_run != MIN_TERMINAL_RUN
        {
            bail!("terminal-tail metadata does not match extraction rule {TERMINAL_TAIL_RULE}");
        }
        if self.sequence_retained || self.quality_retained {
            bail!("terminal-tail v1 must not claim retained sequence or quality");
        }
        let all_zero = self.events == 0 && self.selected_molecules == 0 && self.chunks == 0;
        let all_nonzero = self.events > 0 && self.selected_molecules > 0 && self.chunks > 0;
        if !all_zero && !all_nonzero {
            bail!("terminal-tail metadata has inconsistent zero cardinalities");
        }
        if self.events < self.selected_molecules || u64::from(self.chunks) > self.selected_molecules
        {
            bail!("terminal-tail metadata cardinalities are inconsistent");
        }
        Ok(())
    }
}

/// The bounded observable retained for one selected read-level tail witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TerminalTailSignal {
    pub clip_len: u8,
    pub tail_bases: u8,
    pub terminal_run: u8,
}

impl TerminalTailSignal {
    /// Construct the exact v1 archive observable. Counts above 31 saturate because the public
    /// evidence contract is the packed five-bit triple, not the discarded clip sequence.
    pub fn saturated(clip_len: usize, tail_bases: usize, terminal_run: usize) -> Self {
        Self {
            clip_len: clip_len.min(31) as u8,
            tail_bases: tail_bases.min(31) as u8,
            terminal_run: terminal_run.min(31) as u8,
        }
    }

    pub fn passes_v1(self) -> bool {
        self.clip_len >= MIN_CLIP
            && u16::from(self.tail_bases) * u16::from(MIN_TAIL_DENOMINATOR)
                >= u16::from(self.clip_len) * u16::from(MIN_TAIL_NUMERATOR)
            && self.terminal_run >= MIN_TERMINAL_RUN
    }

    /// Deterministic preference among bounded signals already stored in an archive.
    pub fn stronger_than(self, other: Self) -> bool {
        let left = u16::from(self.tail_bases) * u16::from(other.clip_len.max(1));
        let right = u16::from(other.tail_bases) * u16::from(self.clip_len.max(1));
        (left, self.terminal_run, self.clip_len) > (right, other.terminal_run, other.clip_len)
    }

    pub fn pack(self, reverse: bool) -> Result<u16> {
        if self.clip_len > 31 || self.tail_bases > 31 || self.terminal_run > 31 {
            bail!("terminal-tail signal exceeds its five-bit field");
        }
        if self.tail_bases > self.clip_len
            || self.terminal_run > self.tail_bases
            || self.terminal_run > self.clip_len
        {
            bail!("terminal-tail signal counts are internally inconsistent");
        }
        Ok(((reverse as u16) << 15)
            | (u16::from(self.clip_len) << 10)
            | (u16::from(self.tail_bases) << 5)
            | u16::from(self.terminal_run))
    }

    pub fn unpack(word: u16) -> Result<(bool, Self)> {
        let signal = Self {
            clip_len: ((word >> 10) & 31) as u8,
            tail_bases: ((word >> 5) & 31) as u8,
            terminal_run: (word & 31) as u8,
        };
        signal.pack(word & 0x8000 != 0)?;
        // Saturation can only increase the stored tail fraction of a qualifying raw signal, so
        // every valid packed witness still satisfies the same lower bound. Enforcing it also
        // prevents malformed saturated words from bypassing fail-closed validation.
        if !signal.passes_v1() {
            bail!("terminal-tail signal does not satisfy the declared extraction rule");
        }
        Ok((word & 0x8000 != 0, signal))
    }
}

/// Coordinate envelope and cardinalities for one tail-bearing molecule chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TerminalTailRoute {
    pub chunk: u32,
    pub chrom: u32,
    pub min_anchor: u32,
    pub max_anchor: u32,
    pub selected_molecules: u32,
    pub events: u32,
}

pub fn encode_index(routes: &[TerminalTailRoute]) -> Result<Vec<u8>> {
    validate_routes(routes, None)?;
    let mut raw = TERMINAL_TAIL_INDEX_MAGIC.to_vec();
    put_varint(&mut raw, routes.len() as u64);
    let mut previous = 0u32;
    for (index, route) in routes.iter().enumerate() {
        let delta = if index == 0 {
            route.chunk
        } else {
            route.chunk - previous
        };
        put_varint(&mut raw, u64::from(delta));
        put_varint(&mut raw, u64::from(route.chrom));
        put_varint(&mut raw, u64::from(route.min_anchor));
        put_varint(&mut raw, u64::from(route.max_anchor - route.min_anchor));
        put_varint(&mut raw, u64::from(route.selected_molecules));
        put_varint(&mut raw, u64::from(route.events));
        previous = route.chunk;
    }
    Ok(raw)
}

pub fn decode_index(raw: &[u8], chunk_count: usize) -> Result<Vec<TerminalTailRoute>> {
    if raw.get(..8) != Some(TERMINAL_TAIL_INDEX_MAGIC) {
        bail!("terminal-tail index has bad magic or is truncated");
    }
    let mut cursor = Cursor::new(&raw[8..]);
    let count =
        usize::try_from(cursor.varint()?).context("terminal-tail route count exceeds usize")?;
    if count > chunk_count {
        bail!("terminal-tail index has {count} routes for only {chunk_count} molecule chunks");
    }
    let mut routes = Vec::with_capacity(count);
    let mut previous = 0u32;
    for index in 0..count {
        let delta =
            u32::try_from(cursor.varint()?).context("terminal-tail chunk delta exceeds u32")?;
        if index > 0 && delta == 0 {
            bail!("terminal-tail chunk ids are not strictly increasing");
        }
        let chunk = if index == 0 {
            delta
        } else {
            previous
                .checked_add(delta)
                .context("terminal-tail chunk id overflow")?
        };
        let chrom =
            u32::try_from(cursor.varint()?).context("terminal-tail chromosome exceeds u32")?;
        let min_anchor =
            u32::try_from(cursor.varint()?).context("terminal-tail anchor exceeds u32")?;
        let span =
            u32::try_from(cursor.varint()?).context("terminal-tail anchor span exceeds u32")?;
        let max_anchor = min_anchor
            .checked_add(span)
            .context("terminal-tail anchor extent overflow")?;
        let selected_molecules = u32::try_from(cursor.varint()?)
            .context("terminal-tail selected-molecule count exceeds u32")?;
        let events =
            u32::try_from(cursor.varint()?).context("terminal-tail event count exceeds u32")?;
        routes.push(TerminalTailRoute {
            chunk,
            chrom,
            min_anchor,
            max_anchor,
            selected_molecules,
            events,
        });
        previous = chunk;
    }
    if !cursor.is_empty() {
        bail!("terminal-tail index has trailing bytes");
    }
    validate_routes(&routes, Some(chunk_count))?;
    Ok(routes)
}

fn validate_routes(routes: &[TerminalTailRoute], chunk_count: Option<usize>) -> Result<()> {
    let mut previous = None;
    for route in routes {
        if previous.is_some_and(|value| route.chunk <= value) {
            bail!("terminal-tail chunk ids are not strictly increasing");
        }
        if chunk_count.is_some_and(|count| route.chunk as usize >= count) {
            bail!(
                "terminal-tail route references absent molecule chunk {}",
                route.chunk
            );
        }
        if route.min_anchor > route.max_anchor || route.selected_molecules == 0 || route.events == 0
        {
            bail!(
                "terminal-tail route {} has invalid cardinality or extent",
                route.chunk
            );
        }
        if route.events < route.selected_molecules {
            bail!(
                "terminal-tail route {} has fewer events than selected molecules",
                route.chunk
            );
        }
        previous = Some(route.chunk);
    }
    Ok(())
}

/// One event before its anchor delta is resolved against the selected molecule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EncodedTerminalTailEvent {
    pub anchor_delta: i64,
    pub reverse: bool,
    pub signal: TerminalTailSignal,
}

/// All retained tail events attached to one local molecule ordinal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EncodedTerminalTailMolecule {
    pub local_ordinal: u32,
    pub events: Vec<EncodedTerminalTailEvent>,
}

pub fn encode_chunk(
    molecule_count: u32,
    molecules: &[EncodedTerminalTailMolecule],
) -> Result<Vec<u8>> {
    validate_chunk(molecule_count, molecules)?;
    let event_count: usize = molecules.iter().map(|molecule| molecule.events.len()).sum();
    let mut raw = TERMINAL_TAIL_CHUNK_MAGIC.to_vec();
    put_varint(&mut raw, u64::from(molecule_count));
    put_varint(&mut raw, molecules.len() as u64);
    put_varint(&mut raw, event_count as u64);
    let mut previous = 0u32;
    for (index, molecule) in molecules.iter().enumerate() {
        let delta = if index == 0 {
            molecule.local_ordinal
        } else {
            molecule.local_ordinal - previous
        };
        put_varint(&mut raw, u64::from(delta));
        put_varint(&mut raw, molecule.events.len() as u64);
        for event in &molecule.events {
            put_svarint(&mut raw, event.anchor_delta);
            raw.extend_from_slice(&event.signal.pack(event.reverse)?.to_le_bytes());
        }
        previous = molecule.local_ordinal;
    }
    Ok(raw)
}

pub fn decode_chunk(
    raw: &[u8],
    expected_molecule_count: u32,
) -> Result<Vec<EncodedTerminalTailMolecule>> {
    if raw.get(..8) != Some(TERMINAL_TAIL_CHUNK_MAGIC) {
        bail!("terminal-tail chunk has bad magic or is truncated");
    }
    let mut cursor = Cursor::new(&raw[8..]);
    let molecule_count =
        u32::try_from(cursor.varint()?).context("terminal-tail molecule count exceeds u32")?;
    if molecule_count != expected_molecule_count {
        bail!(
            "terminal-tail chunk declares {molecule_count} molecules; ordinary chunk has {expected_molecule_count}"
        );
    }
    let selected = usize::try_from(cursor.varint()?)
        .context("terminal-tail selected-molecule count exceeds usize")?;
    if selected > molecule_count as usize {
        bail!("terminal-tail selected-molecule count exceeds ordinary chunk size");
    }
    let expected_events =
        usize::try_from(cursor.varint()?).context("terminal-tail event count exceeds usize")?;
    if expected_events < selected {
        bail!("terminal-tail chunk has fewer events than selected molecules");
    }
    // Each selected molecule needs at least an ordinal, event count, signed delta, and u16.
    if selected > cursor_remaining(raw, &cursor) / 5
        || expected_events > cursor_remaining(raw, &cursor) / 3
    {
        bail!("terminal-tail chunk cardinality exceeds its remaining bytes");
    }
    let mut molecules = Vec::with_capacity(selected);
    let mut previous = 0u32;
    let mut decoded_events = 0usize;
    for index in 0..selected {
        let delta =
            u32::try_from(cursor.varint()?).context("terminal-tail ordinal delta exceeds u32")?;
        if index > 0 && delta == 0 {
            bail!("terminal-tail molecule ordinals are not strictly increasing");
        }
        let local_ordinal = if index == 0 {
            delta
        } else {
            previous
                .checked_add(delta)
                .context("terminal-tail molecule ordinal overflow")?
        };
        if local_ordinal >= molecule_count {
            bail!("terminal-tail molecule ordinal {local_ordinal} is out of range");
        }
        let count = usize::try_from(cursor.varint()?)
            .context("terminal-tail per-molecule event count exceeds usize")?;
        if count == 0 || count > cursor_remaining(raw, &cursor) / 3 {
            bail!("terminal-tail molecule {local_ordinal} has invalid event count {count}");
        }
        decoded_events = decoded_events
            .checked_add(count)
            .context("terminal-tail event count overflow")?;
        if decoded_events > expected_events {
            bail!("terminal-tail signal count exceeds declared event count");
        }
        let mut events = Vec::with_capacity(count);
        for _ in 0..count {
            let anchor_delta = cursor.svarint()?;
            let bytes = cursor.take(2)?;
            let (reverse, signal) =
                TerminalTailSignal::unpack(u16::from_le_bytes([bytes[0], bytes[1]]))?;
            events.push(EncodedTerminalTailEvent {
                anchor_delta,
                reverse,
                signal,
            });
        }
        if events
            .windows(2)
            .any(|pair| pair[0].anchor_delta >= pair[1].anchor_delta)
        {
            bail!("terminal-tail molecule {local_ordinal} has noncanonical event order");
        }
        molecules.push(EncodedTerminalTailMolecule {
            local_ordinal,
            events,
        });
        previous = local_ordinal;
    }
    if decoded_events != expected_events {
        bail!("terminal-tail chunk decoded {decoded_events} events; header declares {expected_events}");
    }
    if !cursor.is_empty() {
        bail!("terminal-tail chunk has trailing bytes");
    }
    validate_chunk(molecule_count, &molecules)?;
    Ok(molecules)
}

fn cursor_remaining(raw: &[u8], cursor: &Cursor<'_>) -> usize {
    raw.len().saturating_sub(8 + cursor.position())
}

fn validate_chunk(molecule_count: u32, molecules: &[EncodedTerminalTailMolecule]) -> Result<()> {
    let mut previous = None;
    for molecule in molecules {
        if molecule.local_ordinal >= molecule_count {
            bail!(
                "terminal-tail molecule ordinal {} is out of range",
                molecule.local_ordinal
            );
        }
        if previous.is_some_and(|value| molecule.local_ordinal <= value) {
            bail!("terminal-tail molecule ordinals are not strictly increasing");
        }
        if molecule.events.is_empty() {
            bail!("terminal-tail selected molecule has no events");
        }
        for pair in molecule.events.windows(2) {
            if pair[0].anchor_delta >= pair[1].anchor_delta {
                bail!("terminal-tail events must have strictly increasing anchors per molecule");
            }
        }
        for event in &molecule.events {
            event.signal.pack(event.reverse)?;
            if !event.signal.passes_v1() {
                bail!("terminal-tail signal does not satisfy the declared extraction rule");
            }
        }
        previous = Some(molecule.local_ordinal);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn signal() -> TerminalTailSignal {
        TerminalTailSignal::saturated(10, 9, 6)
    }

    #[test]
    fn signal_roundtrip_and_threshold_are_exact() {
        let signal = signal();
        assert!(signal.passes_v1());
        assert_eq!(
            TerminalTailSignal::unpack(signal.pack(false).unwrap()).unwrap(),
            (false, signal)
        );
        assert_eq!(
            TerminalTailSignal::unpack(signal.pack(true).unwrap()).unwrap(),
            (true, signal)
        );
        assert!(!TerminalTailSignal::saturated(10, 7, 6).passes_v1());
        assert!(TerminalTailSignal::saturated(10, 8, 4).passes_v1());
        let malformed = (31u16 << 10) | 4;
        assert!(TerminalTailSignal::unpack(malformed).is_err());
        let impossible_run = (10u16 << 10) | (8u16 << 5) | 9;
        assert!(TerminalTailSignal::unpack(impossible_run).is_err());
    }

    #[test]
    fn sparse_index_and_one_to_many_chunk_roundtrip() {
        let routes = vec![TerminalTailRoute {
            chunk: 2,
            chrom: 1,
            min_anchor: 100,
            max_anchor: 180,
            selected_molecules: 1,
            events: 2,
        }];
        assert_eq!(
            decode_index(&encode_index(&routes).unwrap(), 3).unwrap(),
            routes
        );

        let molecules = vec![EncodedTerminalTailMolecule {
            local_ordinal: 3,
            events: vec![
                EncodedTerminalTailEvent {
                    anchor_delta: -5,
                    reverse: false,
                    signal: signal(),
                },
                EncodedTerminalTailEvent {
                    anchor_delta: 75,
                    reverse: false,
                    signal: signal(),
                },
            ],
        }];
        assert_eq!(
            decode_chunk(&encode_chunk(5, &molecules).unwrap(), 5).unwrap(),
            molecules
        );
    }

    #[test]
    fn malformed_sparse_layouts_fail_closed() {
        let duplicate = vec![
            EncodedTerminalTailMolecule {
                local_ordinal: 1,
                events: vec![EncodedTerminalTailEvent {
                    anchor_delta: 0,
                    reverse: false,
                    signal: signal(),
                }],
            },
            EncodedTerminalTailMolecule {
                local_ordinal: 1,
                events: vec![EncodedTerminalTailEvent {
                    anchor_delta: 1,
                    reverse: false,
                    signal: signal(),
                }],
            },
        ];
        assert!(encode_chunk(3, &duplicate).is_err());

        let valid = vec![EncodedTerminalTailMolecule {
            local_ordinal: 0,
            events: vec![EncodedTerminalTailEvent {
                anchor_delta: 0,
                reverse: false,
                signal: signal(),
            }],
        }];
        let mut trailing = encode_chunk(1, &valid).unwrap();
        trailing.push(0);
        assert!(decode_chunk(&trailing, 1).is_err());
        assert!(decode_chunk(&encode_chunk(1, &valid).unwrap(), 2).is_err());
    }
}
