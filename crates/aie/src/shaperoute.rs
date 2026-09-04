//! Exact, derived routing from an intron span to local archive shape offsets.
//!
//! A route is only a candidate filter. Molecule chunks remain authoritative, and callers must
//! use [`SpanRoute::matches`] to verify both genomic coordinates for every representative or
//! alternative placement they examine.

use anyhow::{bail, Context, Result};

const ROUTE_MAGIC: &[u8; 4] = b"JSHR";
const ROUTE_HEADER_BYTES: usize = 28;
const MAX_SHAPES: usize = 10_000_000;
const MAX_ROUTE_PAIRS: usize = 10_000_000;

pub(crate) const ROUTE_CODEC_VERSION: u32 = 1;
pub(crate) const MAX_SPANS_PER_BLOCK: usize = 256;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct RoutePair {
    pub(crate) shape_id: u32,
    pub(crate) donor_offset: u32,
}

impl RoutePair {
    pub(crate) fn acceptor_offset(self, exact_span: u32) -> Result<u32> {
        self.donor_offset
            .checked_add(exact_span)
            .context("shape-route acceptor offset overflow")
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SpanRoute {
    pub(crate) span: u32,
    pub(crate) pairs: Vec<RoutePair>,
}

impl SpanRoute {
    pub(crate) fn pairs_for_shape(&self, shape_id: u32) -> &[RoutePair] {
        let begin = self.pairs.partition_point(|pair| pair.shape_id < shape_id);
        let end = self.pairs[begin..].partition_point(|pair| pair.shape_id == shape_id) + begin;
        &self.pairs[begin..end]
    }

    /// Test one decoded placement against this exact-span candidate set. Both additions are
    /// checked independently even though bucket membership makes the second equality redundant
    /// after the first; this preserves explicit two-coordinate verification.
    pub(crate) fn matches(
        &self,
        shape_id: u32,
        position: u32,
        donor: u32,
        acceptor: u32,
    ) -> Result<bool> {
        let requested_span = acceptor
            .checked_sub(donor)
            .context("junction acceptor precedes donor")?;
        if requested_span != self.span {
            return Ok(false);
        }
        for pair in self.pairs_for_shape(shape_id) {
            let acceptor_offset = pair.acceptor_offset(self.span)?;
            let observed_donor = position
                .checked_add(pair.donor_offset)
                .context("shape-route genomic donor overflow")?;
            let observed_acceptor = position
                .checked_add(acceptor_offset)
                .context("shape-route genomic acceptor overflow")?;
            if observed_donor == donor && observed_acceptor == acceptor {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RouteBlock {
    pub(crate) archive_ordinal: u32,
    pub(crate) n_shapes: u32,
    pub(crate) spans: Vec<SpanRoute>,
}

impl RouteBlock {
    pub(crate) fn first_span(&self) -> Option<u32> {
        self.spans.first().map(|row| row.span)
    }

    pub(crate) fn last_span(&self) -> Option<u32> {
        self.spans.last().map(|row| row.span)
    }

    pub(crate) fn span(&self, exact_span: u32) -> Option<&SpanRoute> {
        self.spans
            .binary_search_by_key(&exact_span, |row| row.span)
            .ok()
            .map(|index| &self.spans[index])
    }

    pub(crate) fn pair_count(&self) -> usize {
        self.spans.iter().map(|row| row.pairs.len()).sum()
    }

    pub(crate) fn descriptor(&self) -> Result<RouteBlockDescriptor> {
        validate_block(self)?;
        let first_span = self.first_span().context("shape-route block is empty")?;
        let last_span = self.last_span().context("shape-route block is empty")?;
        Ok(RouteBlockDescriptor {
            first_span,
            last_span,
            section_name: route_section_name(self.archive_ordinal, first_span),
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RouteBlockDescriptor {
    pub(crate) first_span: u32,
    pub(crate) last_span: u32,
    pub(crate) section_name: String,
}

/// Binding stored in a collection manifest. `shapes_digest` is the compressed-payload BLAKE3
/// committed by the rooted source archive's directory entry named `shapes`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ShapeRouteBinding {
    pub(crate) archive_ordinal: u32,
    pub(crate) source_root: [u8; 32],
    pub(crate) shapes_digest: [u8; 32],
    pub(crate) n_shapes: u32,
    pub(crate) blocks: Vec<RouteBlockDescriptor>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DerivedShapeRoutes {
    pub(crate) archive_ordinal: u32,
    pub(crate) n_shapes: u32,
    pub(crate) blocks: Vec<RouteBlock>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct EncodedRouteBlock {
    pub(crate) archive_ordinal: u32,
    pub(crate) n_shapes: u32,
    pub(crate) first_span: u32,
    pub(crate) last_span: u32,
    pub(crate) n_spans: usize,
    pub(crate) n_pairs: usize,
    pub(crate) raw: Vec<u8>,
}

impl EncodedRouteBlock {
    pub(crate) fn descriptor(&self) -> RouteBlockDescriptor {
        RouteBlockDescriptor {
            first_span: self.first_span,
            last_span: self.last_span,
            section_name: route_section_name(self.archive_ordinal, self.first_span),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct EncodedShapeRoutes {
    pub(crate) archive_ordinal: u32,
    pub(crate) n_shapes: u32,
    pub(crate) blocks: Vec<EncodedRouteBlock>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RouteTriple {
    span: u32,
    shape_id: u32,
    donor_offset: u32,
}

#[cfg(test)]
impl DerivedShapeRoutes {
    pub(crate) fn total_spans(&self) -> usize {
        self.blocks.iter().map(|block| block.spans.len()).sum()
    }

    pub(crate) fn total_pairs(&self) -> usize {
        self.blocks.iter().map(RouteBlock::pair_count).sum()
    }

    pub(crate) fn binding(
        &self,
        source_root: [u8; 32],
        shapes_digest: [u8; 32],
    ) -> Result<ShapeRouteBinding> {
        let blocks = self
            .blocks
            .iter()
            .map(RouteBlock::descriptor)
            .collect::<Result<Vec<_>>>()?;
        let binding = ShapeRouteBinding {
            archive_ordinal: self.archive_ordinal,
            source_root,
            shapes_digest,
            n_shapes: self.n_shapes,
            blocks,
        };
        validate_binding(&binding, self.archive_ordinal, source_root, shapes_digest)?;
        Ok(binding)
    }
}

pub(crate) fn route_section_name(archive_ordinal: u32, first_span: u32) -> String {
    format!("s.{archive_ordinal}.{first_span}")
}

pub(crate) fn parse_route_section_name(name: &str) -> Result<(u32, u32)> {
    let mut fields = name.split('.');
    if fields.next() != Some("s") {
        bail!("shape-route section {name} does not begin with s");
    }
    let archive = fields
        .next()
        .context("shape-route section lacks an archive ordinal")?
        .parse::<u32>()
        .with_context(|| format!("shape-route section {name} has an invalid archive ordinal"))?;
    let first_span = fields
        .next()
        .context("shape-route section lacks a first span")?
        .parse::<u32>()
        .with_context(|| format!("shape-route section {name} has an invalid first span"))?;
    if fields.next().is_some() {
        bail!("shape-route section {name} has trailing name fields");
    }
    let canonical = route_section_name(archive, first_span);
    if name != canonical {
        bail!("shape-route section {name} is not canonically named {canonical}");
    }
    Ok((archive, first_span))
}

pub(crate) fn validate_binding(
    binding: &ShapeRouteBinding,
    expected_archive: u32,
    expected_root: [u8; 32],
    expected_shapes_digest: [u8; 32],
) -> Result<()> {
    if binding.archive_ordinal != expected_archive {
        bail!(
            "shape-route binding archive ordinal {} differs from expected {expected_archive}",
            binding.archive_ordinal
        );
    }
    if binding.source_root != expected_root {
        bail!("shape-route source-root binding mismatch");
    }
    if binding.shapes_digest != expected_shapes_digest {
        bail!("shape-route compressed shapes-digest binding mismatch");
    }
    if binding.n_shapes as usize > MAX_SHAPES {
        bail!(
            "shape-route binding declares {} shapes; safety limit is {MAX_SHAPES}",
            binding.n_shapes
        );
    }
    let mut previous_last = None;
    for descriptor in &binding.blocks {
        if descriptor.first_span == 0 || descriptor.last_span < descriptor.first_span {
            bail!(
                "shape-route block {} has invalid span interval {}..{}",
                descriptor.section_name,
                descriptor.first_span,
                descriptor.last_span
            );
        }
        let (archive, named_first) = parse_route_section_name(&descriptor.section_name)?;
        if archive != binding.archive_ordinal || named_first != descriptor.first_span {
            bail!(
                "shape-route block {} disagrees with its archive/span binding",
                descriptor.section_name
            );
        }
        if previous_last.is_some_and(|last| descriptor.first_span <= last) {
            bail!("shape-route block directory overlaps or is not strictly span-ordered");
        }
        previous_last = Some(descriptor.last_span);
    }
    Ok(())
}

fn validate_block(block: &RouteBlock) -> Result<()> {
    if block.n_shapes as usize > MAX_SHAPES {
        bail!(
            "shape-route block declares {} shapes; safety limit is {MAX_SHAPES}",
            block.n_shapes
        );
    }
    if block.spans.is_empty() || block.spans.len() > MAX_SPANS_PER_BLOCK {
        bail!(
            "shape-route block has {} spans; expected 1..={MAX_SPANS_PER_BLOCK}",
            block.spans.len()
        );
    }
    let mut previous_span = None;
    let mut total_pairs = 0usize;
    for row in &block.spans {
        if row.span == 0 || previous_span.is_some_and(|span| row.span <= span) {
            bail!("shape-route spans are zero, duplicate, or not strictly sorted");
        }
        if row.pairs.is_empty() {
            bail!("shape-route span {} has no candidate pairs", row.span);
        }
        total_pairs = total_pairs
            .checked_add(row.pairs.len())
            .context("shape-route pair count overflow")?;
        if total_pairs > MAX_ROUTE_PAIRS {
            bail!("shape-route block exceeds the {MAX_ROUTE_PAIRS}-pair safety limit");
        }
        let mut previous_pair: Option<RoutePair> = None;
        for &pair in &row.pairs {
            if pair.shape_id >= block.n_shapes {
                bail!(
                    "shape-route shape id {} is outside the {}-shape source dictionary",
                    pair.shape_id,
                    block.n_shapes
                );
            }
            if previous_pair.is_some_and(|previous| pair <= previous) {
                bail!("shape-route pairs are duplicate or not strictly sorted");
            }
            pair.acceptor_offset(row.span)?;
            previous_pair = Some(pair);
        }
        previous_span = Some(row.span);
    }
    Ok(())
}

fn visit_shape_introns(
    shapes_raw: &[u8],
    mut visit: impl FnMut(RouteTriple) -> Result<()>,
) -> Result<(u32, usize)> {
    let mut cursor = StrictCursor::new(shapes_raw);
    let mut n_shapes = 0usize;
    let mut total_pairs = 0usize;
    while !cursor.is_empty() {
        if n_shapes >= MAX_SHAPES {
            bail!("source shapes exceed the {MAX_SHAPES}-shape safety limit");
        }
        let n_blocks = usize::try_from(cursor.varint("shape block count")?)
            .context("shape block count exceeds usize")?;
        if n_blocks > cursor.remaining() / 2 {
            bail!("shape block count exceeds the remaining shapes payload");
        }
        let shape_id = u32::try_from(n_shapes).context("shape id exceeds u32")?;
        let mut previous_end = None;
        for block_index in 0..n_blocks {
            let offset = u32::try_from(cursor.varint("shape block offset")?)
                .context("shape block offset exceeds u32")?;
            let length = u32::try_from(cursor.varint("shape block length")?)
                .context("shape block length exceeds u32")?;
            if length == 0 {
                bail!("shape {shape_id} block {block_index} has zero length");
            }
            let end = offset
                .checked_add(length)
                .with_context(|| format!("shape {shape_id} block end overflows u32"))?;
            if block_index == 0 && offset != 0 {
                bail!("shape {shape_id} first block does not begin at offset zero");
            }
            if let Some(donor_offset) = previous_end {
                let span = offset.checked_sub(donor_offset).with_context(|| {
                    format!("shape {shape_id} has reversed or overlapping blocks")
                })?;
                if span == 0 {
                    bail!("shape {shape_id} contains a zero-span intron");
                }
                RoutePair {
                    shape_id,
                    donor_offset,
                }
                .acceptor_offset(span)?;
                total_pairs = total_pairs
                    .checked_add(1)
                    .context("derived shape-route pair count overflow")?;
                if total_pairs > MAX_ROUTE_PAIRS {
                    bail!("derived shape routes exceed the {MAX_ROUTE_PAIRS}-pair safety limit");
                }
                visit(RouteTriple {
                    span,
                    shape_id,
                    donor_offset,
                })?;
            }
            previous_end = Some(end);
        }
        n_shapes += 1;
    }
    Ok((
        u32::try_from(n_shapes).context("source shape count exceeds u32")?,
        total_pairs,
    ))
}

/// Derive the exact canonical route sections while retaining one contiguous sortable table and
/// the encoded output. This is the collection-build representation: unlike a map of per-span
/// vectors, it has no allocation per distinct span and it can be written without reconstructing
/// the decoded route hierarchy.
pub(crate) fn derive_encoded_from_shapes(
    shapes_raw: &[u8],
    archive_ordinal: u32,
) -> Result<EncodedShapeRoutes> {
    let (expected_shapes, expected_pairs) = visit_shape_introns(shapes_raw, |_| Ok(()))?;
    let mut routes = Vec::new();
    routes
        .try_reserve_exact(expected_pairs)
        .context("allocating contiguous shape-route table")?;
    let (n_shapes, observed_pairs) = visit_shape_introns(shapes_raw, |route| {
        routes.push(route);
        Ok(())
    })?;
    if n_shapes != expected_shapes
        || observed_pairs != expected_pairs
        || routes.len() != expected_pairs
    {
        bail!("shape-route source changed between deterministic derivation passes");
    }
    routes.sort_unstable();
    if let Some(pair) = routes.windows(2).find(|pair| pair[0] == pair[1]) {
        bail!(
            "source shapes derive a duplicate route pair for span {}",
            pair[0].span
        );
    }

    let mut blocks = Vec::new();
    let mut begin = 0usize;
    while begin < routes.len() {
        let first_span = routes[begin].span;
        let mut end = begin;
        let mut n_spans = 0usize;
        let mut last_span = first_span;
        while end < routes.len() && n_spans < MAX_SPANS_PER_BLOCK {
            last_span = routes[end].span;
            n_spans += 1;
            end += 1;
            while end < routes.len() && routes[end].span == last_span {
                end += 1;
            }
        }
        let n_spans_u32 = u32::try_from(n_spans).context("shape-route span count exceeds u32")?;
        let mut raw = Vec::new();
        raw.extend_from_slice(ROUTE_MAGIC);
        for value in [
            ROUTE_CODEC_VERSION,
            archive_ordinal,
            n_shapes,
            first_span,
            last_span,
            n_spans_u32,
        ] {
            raw.extend_from_slice(&value.to_le_bytes());
        }
        let mut row_begin = begin;
        let mut previous_span = 0u32;
        while row_begin < end {
            let span = routes[row_begin].span;
            let mut row_end = row_begin + 1;
            while row_end < end && routes[row_end].span == span {
                row_end += 1;
            }
            let span_code = if row_begin == begin {
                span
            } else {
                span.checked_sub(previous_span)
                    .context("shape-route spans are not strictly sorted")?
            };
            put_varint(&mut raw, span_code as u64);
            put_varint(&mut raw, (row_end - row_begin) as u64);
            let mut previous_pair: Option<RoutePair> = None;
            for route in &routes[row_begin..row_end] {
                let pair = RoutePair {
                    shape_id: route.shape_id,
                    donor_offset: route.donor_offset,
                };
                pair.acceptor_offset(span)?;
                let (shape_code, donor_code) = match previous_pair {
                    None => (pair.shape_id, pair.donor_offset),
                    Some(previous) if pair.shape_id == previous.shape_id => (
                        0,
                        pair.donor_offset
                            .checked_sub(previous.donor_offset)
                            .context("shape-route donor offsets are not sorted")?,
                    ),
                    Some(previous) => (
                        pair.shape_id
                            .checked_sub(previous.shape_id)
                            .context("shape-route shape ids are not sorted")?,
                        pair.donor_offset,
                    ),
                };
                put_varint(&mut raw, shape_code as u64);
                put_varint(&mut raw, donor_code as u64);
                previous_pair = Some(pair);
            }
            previous_span = span;
            row_begin = row_end;
        }
        blocks.push(EncodedRouteBlock {
            archive_ordinal,
            n_shapes,
            first_span,
            last_span,
            n_spans,
            n_pairs: end - begin,
            raw,
        });
        begin = end;
    }
    Ok(EncodedShapeRoutes {
        archive_ordinal,
        n_shapes,
        blocks,
    })
}

/// Parse the source `shapes` payload without materializing the shape dictionary. This decoded
/// representation is retained for query verification and focused codec tests; collection builds
/// use [`derive_encoded_from_shapes`] directly.
pub(crate) fn derive_from_shapes(
    shapes_raw: &[u8],
    archive_ordinal: u32,
) -> Result<DerivedShapeRoutes> {
    let encoded = derive_encoded_from_shapes(shapes_raw, archive_ordinal)?;
    let blocks = encoded
        .blocks
        .iter()
        .map(|block| {
            decode_block(
                &block.raw,
                archive_ordinal,
                encoded.n_shapes,
                block.first_span,
                block.last_span,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(DerivedShapeRoutes {
        archive_ordinal,
        n_shapes: encoded.n_shapes,
        blocks,
    })
}

pub(crate) fn encode_block(block: &RouteBlock) -> Result<Vec<u8>> {
    validate_block(block)?;
    let first_span = block.first_span().context("shape-route block is empty")?;
    let last_span = block.last_span().context("shape-route block is empty")?;
    let n_spans = u32::try_from(block.spans.len()).context("shape-route span count exceeds u32")?;
    let mut out = Vec::new();
    out.extend_from_slice(ROUTE_MAGIC);
    for value in [
        ROUTE_CODEC_VERSION,
        block.archive_ordinal,
        block.n_shapes,
        first_span,
        last_span,
        n_spans,
    ] {
        out.extend_from_slice(&value.to_le_bytes());
    }
    let mut previous_span = 0u32;
    for (span_index, row) in block.spans.iter().enumerate() {
        let span_code = if span_index == 0 {
            row.span
        } else {
            row.span - previous_span
        };
        put_varint(&mut out, span_code as u64);
        put_varint(&mut out, row.pairs.len() as u64);
        let mut previous_pair: Option<RoutePair> = None;
        for &pair in &row.pairs {
            let (shape_code, donor_code) = match previous_pair {
                None => (pair.shape_id, pair.donor_offset),
                Some(previous) if pair.shape_id == previous.shape_id => (
                    0,
                    pair.donor_offset
                        .checked_sub(previous.donor_offset)
                        .context("shape-route donor offsets are not sorted")?,
                ),
                Some(previous) => (
                    pair.shape_id
                        .checked_sub(previous.shape_id)
                        .context("shape-route shape ids are not sorted")?,
                    pair.donor_offset,
                ),
            };
            put_varint(&mut out, shape_code as u64);
            put_varint(&mut out, donor_code as u64);
            previous_pair = Some(pair);
        }
        previous_span = row.span;
    }
    Ok(out)
}

pub(crate) fn decode_block(
    raw: &[u8],
    expected_archive: u32,
    expected_n_shapes: u32,
    expected_first: u32,
    expected_last: u32,
) -> Result<RouteBlock> {
    if raw.len() < ROUTE_HEADER_BYTES {
        bail!("shape-route block is truncated before its fixed header");
    }
    if &raw[..4] != ROUTE_MAGIC {
        bail!("shape-route block has bad magic");
    }
    let mut cursor = StrictCursor::at(raw, 4);
    let version = cursor.u32("shape-route version")?;
    if version != ROUTE_CODEC_VERSION {
        bail!("unsupported shape-route codec version {version}; expected {ROUTE_CODEC_VERSION}");
    }
    let archive_ordinal = cursor.u32("shape-route archive ordinal")?;
    let n_shapes = cursor.u32("shape-route source shape count")?;
    let first_span = cursor.u32("shape-route first span")?;
    let last_span = cursor.u32("shape-route last span")?;
    let n_spans = usize::try_from(cursor.u32("shape-route span count")?)
        .context("shape-route span count exceeds usize")?;
    if archive_ordinal != expected_archive {
        bail!(
            "shape-route block archive ordinal {archive_ordinal} differs from expected {expected_archive}"
        );
    }
    if n_shapes != expected_n_shapes {
        bail!("shape-route block shape count {n_shapes} differs from expected {expected_n_shapes}");
    }
    if first_span != expected_first || last_span != expected_last {
        bail!(
            "shape-route block span bounds {first_span}..{last_span} differ from expected {expected_first}..{expected_last}"
        );
    }
    if n_spans == 0 || n_spans > MAX_SPANS_PER_BLOCK {
        bail!("shape-route block span count {n_spans} is outside 1..={MAX_SPANS_PER_BLOCK}");
    }
    let mut spans = Vec::with_capacity(n_spans);
    let mut previous_span = 0u32;
    let mut total_pairs = 0usize;
    for span_index in 0..n_spans {
        let span_code = u32::try_from(cursor.varint("shape-route span")?)
            .context("shape-route span exceeds u32")?;
        let span = if span_index == 0 {
            span_code
        } else {
            if span_code == 0 {
                bail!("shape-route spans are duplicate or not strictly sorted");
            }
            previous_span
                .checked_add(span_code)
                .context("shape-route span delta overflows u32")?
        };
        let pair_count = usize::try_from(cursor.varint("shape-route pair count")?)
            .context("shape-route pair count exceeds usize")?;
        if pair_count == 0 {
            bail!("shape-route span {span} has no candidate pairs");
        }
        if pair_count > cursor.remaining() / 2 {
            bail!("shape-route pair count exceeds remaining block bytes");
        }
        total_pairs = total_pairs
            .checked_add(pair_count)
            .context("shape-route pair count overflow")?;
        if total_pairs > MAX_ROUTE_PAIRS {
            bail!("shape-route block exceeds the {MAX_ROUTE_PAIRS}-pair safety limit");
        }
        let mut pairs = Vec::new();
        pairs
            .try_reserve(pair_count.min(4_096))
            .context("allocating shape-route pairs")?;
        let mut previous_pair: Option<RoutePair> = None;
        for pair_index in 0..pair_count {
            let shape_code = u32::try_from(cursor.varint("shape-route shape id")?)
                .context("shape-route shape id exceeds u32")?;
            let donor_code = u32::try_from(cursor.varint("shape-route donor offset")?)
                .context("shape-route donor offset exceeds u32")?;
            let pair = match previous_pair {
                None => RoutePair {
                    shape_id: shape_code,
                    donor_offset: donor_code,
                },
                Some(previous) if shape_code == 0 => {
                    if donor_code == 0 {
                        bail!("shape-route pairs contain a duplicate donor offset");
                    }
                    RoutePair {
                        shape_id: previous.shape_id,
                        donor_offset: previous
                            .donor_offset
                            .checked_add(donor_code)
                            .context("shape-route donor delta overflows u32")?,
                    }
                }
                Some(previous) => RoutePair {
                    shape_id: previous
                        .shape_id
                        .checked_add(shape_code)
                        .context("shape-route shape-id delta overflows u32")?,
                    donor_offset: donor_code,
                },
            };
            if pair_index > 0 && previous_pair.is_some_and(|previous| pair <= previous) {
                bail!("shape-route pairs are duplicate or not strictly sorted");
            }
            if pair.shape_id >= n_shapes {
                bail!(
                    "shape-route shape id {} is outside the {n_shapes}-shape source dictionary",
                    pair.shape_id
                );
            }
            pair.acceptor_offset(span)?;
            pairs.push(pair);
            previous_pair = Some(pair);
        }
        spans.push(SpanRoute { span, pairs });
        previous_span = span;
    }
    if !cursor.is_empty() {
        bail!(
            "shape-route block has {} trailing byte(s)",
            cursor.remaining()
        );
    }
    let block = RouteBlock {
        archive_ordinal,
        n_shapes,
        spans,
    };
    validate_block(&block)?;
    if block.first_span() != Some(first_span) || block.last_span() != Some(last_span) {
        bail!("shape-route decoded rows disagree with the declared span bounds");
    }
    if encode_block(&block)? != raw {
        bail!("shape-route block is not canonically encoded");
    }
    Ok(block)
}

pub(crate) fn validate_blocks_against_binding(
    binding: &ShapeRouteBinding,
    blocks: &[RouteBlock],
) -> Result<()> {
    if blocks.len() != binding.blocks.len() {
        bail!("shape-route block count differs from its manifest directory");
    }
    for (block, descriptor) in blocks.iter().zip(&binding.blocks) {
        if block.archive_ordinal != binding.archive_ordinal || block.n_shapes != binding.n_shapes {
            bail!("shape-route block disagrees with its archive/shape-count binding");
        }
        if &block.descriptor()? != descriptor {
            bail!("shape-route block disagrees with its manifest span directory");
        }
    }
    Ok(())
}

pub(crate) fn verify_reconstruction(
    shapes_raw: &[u8],
    binding: &ShapeRouteBinding,
    blocks: &[RouteBlock],
) -> Result<()> {
    validate_blocks_against_binding(binding, blocks)?;
    let derived = derive_from_shapes(shapes_raw, binding.archive_ordinal)?;
    if derived.n_shapes != binding.n_shapes {
        bail!("reconstructed source shape count differs from the route binding");
    }
    if derived.blocks != blocks {
        bail!("shape-route state differs from exact reconstruction of the bound shapes section");
    }
    Ok(())
}

struct StrictCursor<'a> {
    raw: &'a [u8],
    position: usize,
}

impl<'a> StrictCursor<'a> {
    fn new(raw: &'a [u8]) -> Self {
        Self { raw, position: 0 }
    }

    fn at(raw: &'a [u8], position: usize) -> Self {
        Self { raw, position }
    }

    fn remaining(&self) -> usize {
        self.raw.len().saturating_sub(self.position)
    }

    fn is_empty(&self) -> bool {
        self.position == self.raw.len()
    }

    fn u32(&mut self, label: &str) -> Result<u32> {
        let end = self
            .position
            .checked_add(4)
            .context("shape-route cursor position overflow")?;
        let bytes: [u8; 4] = self
            .raw
            .get(self.position..end)
            .with_context(|| format!("{label} is truncated"))?
            .try_into()
            .unwrap();
        self.position = end;
        Ok(u32::from_le_bytes(bytes))
    }

    fn varint(&mut self, label: &str) -> Result<u64> {
        let start = self.position;
        let mut value = 0u64;
        let mut shift = 0u32;
        loop {
            let byte = *self
                .raw
                .get(self.position)
                .with_context(|| format!("{label} varint is truncated"))?;
            self.position += 1;
            if shift == 63 && byte & 0x7e != 0 {
                bail!("{label} varint overflows u64");
            }
            value |= ((byte & 0x7f) as u64) << shift;
            if byte & 0x80 == 0 {
                break;
            }
            if shift >= 63 {
                bail!("{label} varint exceeds ten bytes");
            }
            shift += 7;
        }
        if self.position - start != varint_len(value) {
            bail!("{label} uses a noncanonical varint");
        }
        Ok(value)
    }
}

fn varint_len(mut value: u64) -> usize {
    let mut len = 1;
    while value >= 0x80 {
        value >>= 7;
        len += 1;
    }
    len
}

fn put_varint(out: &mut Vec<u8>, mut value: u64) {
    loop {
        let byte = (value & 0x7f) as u8;
        value >>= 7;
        if value == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_shapes(shapes: &[&[(u32, u32)]]) -> Vec<u8> {
        let mut raw = Vec::new();
        for shape in shapes {
            put_varint(&mut raw, shape.len() as u64);
            for &(offset, length) in *shape {
                put_varint(&mut raw, offset as u64);
                put_varint(&mut raw, length as u64);
            }
        }
        raw
    }

    fn fixture_block() -> RouteBlock {
        RouteBlock {
            archive_ordinal: 3,
            n_shapes: 4,
            spans: vec![
                SpanRoute {
                    span: 10,
                    pairs: vec![
                        RoutePair {
                            shape_id: 0,
                            donor_offset: 10,
                        },
                        RoutePair {
                            shape_id: 0,
                            donor_offset: 25,
                        },
                        RoutePair {
                            shape_id: 1,
                            donor_offset: 8,
                        },
                    ],
                },
                SpanRoute {
                    span: 200,
                    pairs: vec![RoutePair {
                        shape_id: 3,
                        donor_offset: 40,
                    }],
                },
            ],
        }
    }

    #[test]
    fn streaming_derivation_preserves_repeated_spans_offsets_and_shape_ids() {
        let raw = encode_shapes(&[&[(0, 10), (20, 5), (35, 5)], &[(0, 8), (18, 4)], &[(0, 50)]]);
        let encoded = derive_encoded_from_shapes(&raw, 7).unwrap();
        let derived = derive_from_shapes(&raw, 7).unwrap();
        assert_eq!(encoded.blocks.len(), derived.blocks.len());
        for (encoded_block, decoded_block) in encoded.blocks.iter().zip(&derived.blocks) {
            assert_eq!(encoded_block.raw, encode_block(decoded_block).unwrap());
        }
        assert_eq!(derived.archive_ordinal, 7);
        assert_eq!(derived.n_shapes, 3);
        assert_eq!(derived.total_spans(), 1);
        assert_eq!(derived.total_pairs(), 3);
        assert_eq!(
            derived.blocks[0].spans[0].pairs,
            vec![
                RoutePair {
                    shape_id: 0,
                    donor_offset: 10,
                },
                RoutePair {
                    shape_id: 0,
                    donor_offset: 25,
                },
                RoutePair {
                    shape_id: 1,
                    donor_offset: 8,
                },
            ]
        );
        let route = &derived.blocks[0].spans[0];
        assert!(route.matches(0, 100, 110, 120).unwrap());
        assert!(route.matches(0, 100, 125, 135).unwrap());
        assert!(route.matches(1, 200, 208, 218).unwrap());
        assert!(!route.matches(2, 100, 110, 120).unwrap());
        assert!(!route.matches(0, 100, 111, 121).unwrap());
    }

    #[test]
    fn single_block_shape_has_no_route_and_distinct_spans_partition_at_256() {
        let single = encode_shapes(&[&[(0, 75)]]);
        let derived = derive_from_shapes(&single, 0).unwrap();
        assert_eq!(derived.n_shapes, 1);
        assert!(derived.blocks.is_empty());

        let shapes: Vec<Vec<(u32, u32)>> =
            (1..=257).map(|span| vec![(0, 1), (1 + span, 1)]).collect();
        let refs: Vec<&[(u32, u32)]> = shapes.iter().map(Vec::as_slice).collect();
        let derived = derive_from_shapes(&encode_shapes(&refs), 12).unwrap();
        assert_eq!(derived.blocks.len(), 2);
        assert_eq!(derived.blocks[0].spans.len(), 256);
        assert_eq!(derived.blocks[0].first_span(), Some(1));
        assert_eq!(derived.blocks[0].last_span(), Some(256));
        assert_eq!(derived.blocks[1].spans.len(), 1);
        assert_eq!(derived.blocks[1].first_span(), Some(257));
        assert_eq!(
            derived.blocks[1].descriptor().unwrap().section_name,
            "s.12.257"
        );
    }

    #[test]
    fn block_codec_is_canonical_bounded_and_context_bound() {
        let block = fixture_block();
        let raw = encode_block(&block).unwrap();
        assert_eq!(decode_block(&raw, 3, 4, 10, 200).unwrap(), fixture_block());
        assert!(decode_block(&raw, 2, 4, 10, 200).is_err());
        assert!(decode_block(&raw, 3, 5, 10, 200).is_err());
        assert!(decode_block(&raw, 3, 4, 11, 200).is_err());

        let mut unknown = raw.clone();
        unknown[4..8].copy_from_slice(&2u32.to_le_bytes());
        assert!(decode_block(&unknown, 3, 4, 10, 200)
            .unwrap_err()
            .to_string()
            .contains("unsupported shape-route codec version"));

        let mut trailing = raw.clone();
        trailing.push(0);
        assert!(decode_block(&trailing, 3, 4, 10, 200).is_err());
        for end in 0..raw.len() {
            assert!(decode_block(&raw[..end], 3, 4, 10, 200).is_err());
        }

        // First row's span occupies byte 28 and its pair count byte 29 in this fixture. Replace
        // the canonical one-byte count with the two-byte spelling of the same value.
        let mut noncanonical = raw.clone();
        noncanonical.splice(29..30, [0x83, 0x00]);
        assert!(decode_block(&noncanonical, 3, 4, 10, 200)
            .unwrap_err()
            .to_string()
            .contains("noncanonical varint"));

        // The first pair begins at byte 30 for this fixture. The source dictionary contains
        // shapes 0..4, so shape id 4 is out of range even though its one-byte spelling is valid.
        let mut out_of_range = raw.clone();
        out_of_range[30] = 4;
        assert!(decode_block(&out_of_range, 3, 4, 10, 200)
            .unwrap_err()
            .to_string()
            .contains("outside the 4-shape"));

        // The second same-shape pair's donor-delta byte is 33. Zero repeats the preceding tuple.
        let mut duplicate = raw.clone();
        duplicate[33] = 0;
        assert!(decode_block(&duplicate, 3, 4, 10, 200)
            .unwrap_err()
            .to_string()
            .contains("duplicate donor"));
    }

    #[test]
    fn malformed_pairs_shapes_and_checked_additions_are_rejected() {
        let mut duplicate = fixture_block();
        duplicate.spans[0].pairs[1] = duplicate.spans[0].pairs[0];
        assert!(encode_block(&duplicate).is_err());

        let mut unsorted = fixture_block();
        unsorted.spans[0].pairs.swap(0, 1);
        assert!(encode_block(&unsorted)
            .unwrap_err()
            .to_string()
            .contains("duplicate or not strictly sorted"));

        let too_many_spans = RouteBlock {
            archive_ordinal: 0,
            n_shapes: 1,
            spans: (1..=MAX_SPANS_PER_BLOCK as u32 + 1)
                .map(|span| SpanRoute {
                    span,
                    pairs: vec![RoutePair {
                        shape_id: 0,
                        donor_offset: 0,
                    }],
                })
                .collect(),
        };
        assert!(encode_block(&too_many_spans)
            .unwrap_err()
            .to_string()
            .contains("expected 1..=256"));

        let mut out_of_range = fixture_block();
        out_of_range.spans[0].pairs[0].shape_id = out_of_range.n_shapes;
        assert!(encode_block(&out_of_range).is_err());

        let mut overflow = fixture_block();
        overflow.spans[0].pairs[0].donor_offset = u32::MAX;
        assert!(encode_block(&overflow).is_err());

        let route = SpanRoute {
            span: 10,
            pairs: vec![RoutePair {
                shape_id: 0,
                donor_offset: 20,
            }],
        };
        assert!(route
            .matches(0, u32::MAX - 5, 14, 24)
            .unwrap_err()
            .to_string()
            .contains("genomic donor overflow"));

        let reversed = encode_shapes(&[&[(0, 10), (9, 5)]]);
        assert!(derive_from_shapes(&reversed, 0)
            .unwrap_err()
            .to_string()
            .contains("reversed or overlapping"));
        let zero_span = encode_shapes(&[&[(0, 10), (10, 5)]]);
        assert!(derive_from_shapes(&zero_span, 0)
            .unwrap_err()
            .to_string()
            .contains("zero-span"));
        assert!(derive_from_shapes(&[100], 0)
            .unwrap_err()
            .to_string()
            .contains("exceeds the remaining"));
        assert!(derive_from_shapes(&[0x80, 0], 0)
            .unwrap_err()
            .to_string()
            .contains("noncanonical varint"));
    }

    #[test]
    fn binding_names_directory_and_full_reconstruction_are_exact() {
        let raw = encode_shapes(&[&[(0, 10), (20, 5), (35, 5)], &[(0, 8), (18, 4)]]);
        let derived = derive_from_shapes(&raw, 5).unwrap();
        let root = [1u8; 32];
        let shapes_digest = [2u8; 32];
        let binding = derived.binding(root, shapes_digest).unwrap();
        validate_binding(&binding, 5, root, shapes_digest).unwrap();
        verify_reconstruction(&raw, &binding, &derived.blocks).unwrap();
        assert_eq!(parse_route_section_name("s.5.10").unwrap(), (5, 10));
        assert!(parse_route_section_name("s.05.10").is_err());
        assert!(parse_route_section_name("s.5.010").is_err());

        let mut wrong_root = binding.clone();
        wrong_root.source_root = [9u8; 32];
        assert!(validate_binding(&wrong_root, 5, root, shapes_digest).is_err());

        let mut wrong_name = binding.clone();
        wrong_name.blocks[0].section_name = "s.6.10".into();
        assert!(validate_binding(&wrong_name, 5, root, shapes_digest).is_err());

        let mut overlapping = binding.clone();
        overlapping.blocks.push(overlapping.blocks[0].clone());
        assert!(validate_binding(&overlapping, 5, root, shapes_digest)
            .unwrap_err()
            .to_string()
            .contains("overlaps"));

        let mut wrong_routes = derived.blocks.clone();
        wrong_routes[0].spans[0].pairs.pop();
        assert!(verify_reconstruction(&raw, &binding, &wrong_routes).is_err());
    }
}
