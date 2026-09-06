//! Placement-local geometry matching. No archive-format or index changes.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "snake_case")]
pub(super) enum MatchWithin {
    Record,
    AnyPlacement,
    AllPlacements,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(super) enum Engine {
    Auto,
    Scalar,
    Compiled,
}

pub(super) fn parse_atom(
    name: &str,
    kind: &str,
    locus: &str,
    chrom_names: &[String],
) -> Result<Predicate> {
    let kind = match kind {
        "region" => PredicateKind::Region,
        "junction" => PredicateKind::Junction,
        "terminal" | "terminal-tail" => PredicateKind::Terminal,
        "overlap" => PredicateKind::Overlap,
        "start" => PredicateKind::Start,
        "end" => PredicateKind::End,
        "junction-near" => PredicateKind::JunctionNear,
        "path" => PredicateKind::Path,
        "subpath" => PredicateKind::Subpath,
        _ => bail!("predicate {name} has unknown kind {kind:?}"),
    };
    let (coordinates, threshold) =
        if matches!(kind, PredicateKind::Overlap | PredicateKind::JunctionNear) {
            let (coordinates, value) = locus.rsplit_once('/').context(
                "overlap and junction-near require an explicit /N threshold after the locus/strand",
            )?;
            let threshold: u32 = value
                .parse()
                .context("threshold must be a nonnegative integer number of bases")?;
            if kind == PredicateKind::Overlap && threshold == 0 {
                bail!("overlap threshold must be at least one base");
            }
            (coordinates, threshold)
        } else {
            (locus, 0)
        };
    let mut path = Vec::new();
    let (chrom, start, end, strand_rev) = if matches!(
        kind,
        PredicateKind::Path | PredicateKind::Subpath
    ) {
        let (coordinates, strand) = if let Some(s) = coordinates.strip_suffix(":+") {
            (s, Some(false))
        } else if let Some(s) = coordinates.strip_suffix(":-") {
            (s, Some(true))
        } else {
            (coordinates, None)
        };
        let (chrom, junctions) = coordinates
            .split_once(':')
            .context("path must be chrom:donor-acceptor[,donor-acceptor...][:strand]")?;
        if junctions.split(',').count() > 64 {
            bail!("path exceeds 64 junctions");
        }
        for junction in junctions.split(',') {
            let (_, donor, acceptor) = super::super::parse_locus(&format!("{chrom}:{junction}"))?;
            if path
                .last()
                .is_some_and(|&(_, previous_end)| previous_end >= donor)
            {
                bail!("path junctions must be disjoint and in increasing genomic order on either strand");
            }
            path.push((donor, acceptor));
        }
        (chrom.to_owned(), path[0].0, path.last().unwrap().1, strand)
    } else {
        parse_stranded_locus(coordinates)?
    };
    let chrom_id = chrom_names
        .iter()
        .position(|candidate| candidate == &chrom)
        .with_context(|| format!("predicate {name} names unknown chromosome {chrom}"))?
        as u32;
    Ok(Predicate {
        name: name.to_owned(),
        kind,
        locus: locus.to_owned(),
        chrom,
        chrom_id,
        start,
        end,
        strand_rev,
        threshold,
        path,
    })
}

pub(super) fn is_geometry(predicate: &Predicate, region_match: RegionMatchArg) -> bool {
    predicate.kind != PredicateKind::Terminal
        && (predicate.kind != PredicateKind::Region || region_match == RegionMatchArg::AlignedBlock)
}

pub(super) fn needs_full_geometry(predicate: &Predicate, region_match: RegionMatchArg) -> bool {
    matches!(
        predicate.kind,
        PredicateKind::Overlap | PredicateKind::Start | PredicateKind::End
    ) || (predicate.kind == PredicateKind::Region && region_match == RegionMatchArg::AlignedBlock)
}

fn absolute(position: u32, offset: u32, length: u32) -> Result<u32> {
    position
        .checked_add(offset)
        .and_then(|value| value.checked_add(length))
        .context("shape overflows genomic coordinates")
}

/// Placement-local adapter shared with GQ. Scope, strand/chromosome checks,
/// omitted-geometry proofs and logical work accounting remain with the caller.
pub(crate) struct ShapePredicate(Predicate);

impl ShapePredicate {
    pub(crate) fn new(
        kind: &str,
        start: u32,
        end: u32,
        threshold: u32,
        path: Vec<(u32, u32)>,
    ) -> Self {
        Self(Predicate {
            name: String::new(),
            kind: match kind {
                "overlaps" => PredicateKind::Overlap,
                "start_in" => PredicateKind::Start,
                "end_in" => PredicateKind::End,
                "junction" => PredicateKind::JunctionNear,
                "path" => PredicateKind::Path,
                "subpath" => PredicateKind::Subpath,
                _ => unreachable!("unsupported shared shape predicate"),
            },
            locus: String::new(),
            chrom: String::new(),
            chrom_id: 0,
            start,
            end,
            strand_rev: None,
            threshold,
            path,
        })
    }

    pub(crate) fn matches(&self, position: u32, shape: &Shape) -> Result<bool> {
        matches_shape(&self.0, position, shape)
    }
}

pub(super) fn matches_shape(predicate: &Predicate, position: u32, shape: &Shape) -> Result<bool> {
    match predicate.kind {
        PredicateKind::Region | PredicateKind::Overlap => {
            let minimum = predicate.threshold.max(1);
            let mut matched = false;
            for &(offset, length) in &shape.blocks {
                let start = absolute(position, offset, 0)?;
                let end = absolute(position, offset, length)?;
                matched |= end
                    .min(predicate.end)
                    .saturating_sub(start.max(predicate.start))
                    >= minimum;
            }
            Ok(matched)
        }
        PredicateKind::Start | PredicateKind::End => {
            let block = if predicate.kind == PredicateKind::Start {
                shape.blocks.first()
            } else {
                shape.blocks.last()
            };
            match block {
                Some(&(offset, length)) => Ok(interval_contains(
                    predicate,
                    absolute(
                        position,
                        offset,
                        if predicate.kind == PredicateKind::End {
                            length
                        } else {
                            0
                        },
                    )?,
                )),
                None => Ok(false),
            }
        }
        PredicateKind::Junction
        | PredicateKind::JunctionNear
        | PredicateKind::Path
        | PredicateKind::Subpath => {
            let mut next = 0;
            let mut matched = false;
            for blocks in shape.blocks.windows(2) {
                let pair = (
                    absolute(position, blocks[0].0, blocks[0].1)?,
                    absolute(position, blocks[1].0, 0)?,
                );
                if matches!(
                    predicate.kind,
                    PredicateKind::Junction | PredicateKind::JunctionNear
                ) {
                    matched |= pair.0.abs_diff(predicate.start) <= predicate.threshold
                        && pair.1.abs_diff(predicate.end) <= predicate.threshold;
                } else if !matched {
                    if pair == predicate.path[next] {
                        next += 1;
                    } else if predicate.kind == PredicateKind::Path {
                        next = usize::from(pair == predicate.path[0]);
                    }
                    matched = next == predicate.path.len();
                }
            }
            Ok(matched)
        }
        PredicateKind::Terminal => bail!("terminal evidence is not placement-linked"),
    }
}

pub(super) fn validate_local_expression(
    expression: &Expression,
    predicates: &[Predicate],
    region_match: RegionMatchArg,
) -> Result<()> {
    match expression {
        Expression::Predicate(index) => {
            if !is_geometry(&predicates[*index], region_match) {
                bail!("predicate {} is record-level, not placement-local; terminal events lack placement linkage and anchor regions need --region-match aligned-block", predicates[*index].name);
            }
        }
        Expression::Not(value) => validate_local_expression(value, predicates, region_match)?,
        Expression::And(left, right) | Expression::Or(left, right) => {
            validate_local_expression(left, predicates, region_match)?;
            validate_local_expression(right, predicates, region_match)?;
        }
    }
    Ok(())
}

type Geometry = (u32, u32, bool, u32);
const MAX_CACHE_ENTRIES: usize = 65_536;

pub(super) struct CompiledMatcher<'a, 'b> {
    context: &'a MatchContext<'b>,
    junctions: FxHashMap<(u32, u32, u32), [u64; 2]>,
    other_geometry: Vec<usize>,
    cache: FxHashMap<Geometry, u64>,
    geometry_predicates: usize,
}

impl<'a, 'b> CompiledMatcher<'a, 'b> {
    pub(super) fn new(context: &'a MatchContext<'b>) -> Self {
        let mut junctions: FxHashMap<_, [u64; 2]> = FxHashMap::default();
        let mut other_geometry = Vec::new();
        for (index, predicate) in context.predicates.iter().enumerate() {
            if predicate.kind == PredicateKind::Junction {
                let masks = junctions
                    .entry((predicate.chrom_id, predicate.start, predicate.end))
                    .or_default();
                for (strand, mask) in masks.iter_mut().enumerate() {
                    if strand_matches(predicate.strand_rev, strand == 1) {
                        *mask |= 1 << index;
                    }
                }
            } else if is_geometry(predicate, context.region_match) {
                other_geometry.push(index);
            }
        }
        Self {
            context,
            junctions,
            other_geometry,
            cache: FxHashMap::default(),
            geometry_predicates: context
                .predicates
                .iter()
                .filter(|p| is_geometry(p, context.region_match))
                .count(),
        }
    }

    fn geometry_mask(&mut self, key: Geometry, scalar: bool) -> Result<u64> {
        if !scalar {
            if let Some(&mask) = self.cache.get(&key) {
                return Ok(mask);
            }
        }
        let (chrom, position, strand, shape_id) = key;
        let shape = self
            .context
            .shapes
            .get(shape_id as usize)
            .context("molecule references missing shape")?;
        let mut mask = 0;
        if scalar {
            for (index, predicate) in self.context.predicates.iter().enumerate() {
                if is_geometry(predicate, self.context.region_match)
                    && chrom == predicate.chrom_id
                    && strand_matches(predicate.strand_rev, strand)
                    && matches_shape(predicate, position, shape)?
                {
                    mask |= 1 << index;
                }
            }
            return Ok(mask);
        }
        if !self.junctions.is_empty() {
            for blocks in shape.blocks.windows(2) {
                let donor = absolute(position, blocks[0].0, blocks[0].1)?;
                let acceptor = absolute(position, blocks[1].0, 0)?;
                if let Some(masks) = self.junctions.get(&(chrom, donor, acceptor)) {
                    mask |= masks[usize::from(strand)];
                }
            }
        }
        for &index in &self.other_geometry {
            let predicate = &self.context.predicates[index];
            if chrom == predicate.chrom_id
                && strand_matches(predicate.strand_rev, strand)
                && matches_shape(predicate, position, shape)?
            {
                mask |= 1 << index;
            }
        }
        // Stop inserting at the cap; do not thrash or grow with the archive.
        if self.cache.len() < MAX_CACHE_ENTRIES {
            self.cache.insert(key, mask);
        }
        Ok(mask)
    }

    pub(super) fn masks(
        &mut self,
        molecule: &MolRec,
        tails: &[TailObservation],
        expression: &Expression,
        within: MatchWithin,
        engine: Engine,
    ) -> Result<(u64, u64, Option<TruthValue>)> {
        let context = self.context;
        let scalar = engine == Engine::Scalar
            || (engine == Engine::Auto
                && (context.predicates.len() < 8 || self.geometry_predicates < 7));
        if scalar && within == MatchWithin::Record {
            let (observed, complete) = context.masks(molecule, tails)?;
            return Ok((observed, complete, None));
        }
        let mut observed = 0;
        let mut complete = 0;
        for (index, predicate) in context.predicates.iter().enumerate() {
            if context.absence_is_complete(predicate, molecule) {
                complete |= 1 << index;
            }
            if !is_geometry(predicate, context.region_match)
                && context.matches(predicate, molecule, tails)?
            {
                observed |= 1 << index;
            }
        }
        let mut count = 0;
        let mut selected = if within == MatchWithin::AllPlacements {
            TruthValue::True
        } else {
            TruthValue::False
        };
        let mut inspect = |key| -> Result<()> {
            let mask = self.geometry_mask(key, scalar)?;
            observed |= mask;
            count += 1;
            if within != MatchWithin::Record {
                let value = expression.evaluate(mask, u64::MAX);
                selected = if within == MatchWithin::AllPlacements {
                    selected.and(value)
                } else {
                    selected.or(value)
                };
            }
            Ok(())
        };
        for chain in &molecule.chains {
            for &(position, shape) in &chain.reps {
                inspect((molecule.chrom, position, molecule.strand_rev, shape))?;
            }
        }
        if context.placements != PlacementScopeArg::Unique {
            for &(position, shape, pattern, _) in &molecule.mms {
                if context.placements == PlacementScopeArg::Direct {
                    inspect((molecule.chrom, position, molecule.strand_rev, shape))?;
                } else {
                    let alternatives = context
                        .patterns
                        .context("all-alternative matching requires patterns")?
                        .get(pattern as usize)
                        .context("molecule references missing pattern")?;
                    for alternative in alternatives {
                        let position = u32::try_from(
                            i64::from(position)
                                .checked_add(alternative.offset)
                                .context("alternative coordinate overflow")?,
                        )
                        .context("alternative coordinate out of range")?;
                        inspect((
                            alternative.chrom,
                            position,
                            molecule.strand_rev != alternative.strand_flip,
                            if alternative.shape == SAME_SHAPE {
                                shape
                            } else {
                                alternative.shape
                            },
                        ))?;
                    }
                }
            }
        }
        if count == 0 {
            selected = TruthValue::False;
        }
        // Every representative and omitted member of a unique chain shares its
        // exact junction path. Compact endpoints do not make junction-only
        // quantifiers uncertain; only referenced geometry-sensitive atoms can.
        let omitted =
            within != MatchWithin::Record && expression.referenced_mask() & !complete != 0;
        if omitted
            && ((within == MatchWithin::AnyPlacement && selected != TruthValue::True)
                || (within == MatchWithin::AllPlacements && selected != TruthValue::False))
        {
            selected = TruthValue::Unknown;
        }
        Ok((
            observed,
            complete,
            (within != MatchWithin::Record).then_some(selected),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rows::MolChain;

    fn predicates(values: &[&str]) -> Vec<Predicate> {
        parse_predicates(
            &values.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
            &["chr1".into(), "chr2".into()],
        )
        .unwrap()
    }

    fn molecule(reps: &[(u32, u32)], weight: u32) -> MolRec {
        MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec::smallvec![MolChain {
                weight,
                reps: reps.iter().copied().collect()
            }],
            mms: Default::default(),
        }
    }

    #[test]
    fn compiled_matches_scalar_across_strands_modes_and_geometry_retention() {
        let predicates = predicates(&[
            "u=region:chr1:0-1000",
            "r=region:chr1:100-111:+",
            "j=junction:chr1:110-130",
            "k=junction:chr2:115-135:-",
            "t=terminal:chr1:109-111",
            "o=overlap:chr1:100-108/8",
            "s=start:chr1:100-101",
            "e=end:chr1:140-141",
            "n=junction-near:chr1:112-128/2",
            "p=path:chr1:110-130",
            "q=subpath:chr1:110-130",
            "r2=region:chr2:114-116:-",
        ]);
        let shapes = [
            Shape {
                blocks: vec![(0, 10), (30, 10)],
            },
            Shape {
                blocks: vec![(0, 8)],
            },
        ];
        let patterns = [vec![
            PatAlt {
                chrom: 0,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            },
            PatAlt {
                chrom: 1,
                offset: 5,
                strand_flip: true,
                shape: SAME_SHAPE,
            },
        ]];
        let tails = [TailObservation {
            local_record: 0,
            chrom_id: 0,
            cleavage_anchor: 110,
            strand_rev: false,
        }];
        for region_match in [RegionMatchArg::Anchor, RegionMatchArg::AlignedBlock] {
            for placements in [
                PlacementScopeArg::Unique,
                PlacementScopeArg::Direct,
                PlacementScopeArg::All,
            ] {
                let context = MatchContext {
                    predicates: &predicates,
                    shapes: &shapes,
                    patterns: Some(&patterns),
                    region_match,
                    placements,
                };
                let mut compiled = CompiledMatcher::new(&context);
                for position in 90..125 {
                    for reverse in [false, true] {
                        for weight in [2, 3] {
                            let mut record = molecule(&[(position, 0), (position + 10, 1)], weight);
                            record.strand_rev = reverse;
                            record.mms.push((100, 0, 0, 1));
                            let expected = context.masks(&record, &tails).unwrap();
                            for _ in 0..2 {
                                // Recheck the cached path, too.
                                let (observed, complete, selection) = compiled
                                    .masks(
                                        &record,
                                        &tails,
                                        &Expression::Predicate(0),
                                        MatchWithin::Record,
                                        Engine::Compiled,
                                    )
                                    .unwrap();
                                assert_eq!((observed, complete), expected);
                                assert_eq!(selection, None);
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn placement_quantifiers_do_not_stitch_distinct_geometries() {
        let predicates = predicates(&["a=region:chr1:100-110", "b=region:chr1:200-210"]);
        let shapes = [Shape {
            blocks: vec![(0, 10)],
        }];
        let context = MatchContext {
            predicates: &predicates,
            shapes: &shapes,
            patterns: None,
            region_match: RegionMatchArg::AlignedBlock,
            placements: PlacementScopeArg::Unique,
        };
        let mut matcher = CompiledMatcher::new(&context);
        let expression = Expression::And(
            Box::new(Expression::Predicate(0)),
            Box::new(Expression::Predicate(1)),
        );
        let record = molecule(&[(100, 0), (200, 0)], 2);
        assert_eq!(context.masks(&record, &[]).unwrap(), (3, 3));
        assert_eq!(expression.evaluate(3, 3), TruthValue::True);
        for engine in [Engine::Scalar, Engine::Compiled] {
            assert_eq!(
                matcher
                    .masks(&record, &[], &expression, MatchWithin::AnyPlacement, engine)
                    .unwrap()
                    .2,
                Some(TruthValue::False)
            );
            let incomplete = molecule(&[(100, 0), (200, 0)], 3);
            assert_eq!(
                matcher
                    .masks(
                        &incomplete,
                        &[],
                        &expression,
                        MatchWithin::AnyPlacement,
                        engine
                    )
                    .unwrap()
                    .2,
                Some(TruthValue::Unknown)
            );
            let either = Expression::Or(
                Box::new(Expression::Predicate(0)),
                Box::new(Expression::Predicate(1)),
            );
            assert_eq!(
                matcher
                    .masks(&record, &[], &either, MatchWithin::AllPlacements, engine)
                    .unwrap()
                    .2,
                Some(TruthValue::True)
            );
            assert_eq!(
                matcher
                    .masks(
                        &incomplete,
                        &[],
                        &either,
                        MatchWithin::AllPlacements,
                        engine
                    )
                    .unwrap()
                    .2,
                Some(TruthValue::Unknown)
            );
            let empty = MolRec {
                chains: Default::default(),
                ..record.clone()
            };
            assert_eq!(
                matcher
                    .masks(&empty, &[], &either, MatchWithin::AllPlacements, engine)
                    .unwrap()
                    .2,
                Some(TruthValue::False)
            );
        }
    }

    #[test]
    fn paths_are_observed_ordered_and_optionally_consecutive() {
        let predicates = predicates(&[
            "p=path:chr1:110-130,170-190",
            "s=subpath:chr1:110-130,170-190",
            "q=path:chr1:110-130,140-160",
            "one=path:chr1:140-160",
        ]);
        let shape = Shape {
            blocks: vec![(0, 10), (30, 10), (60, 10), (90, 10)],
        };
        assert!(!matches_shape(&predicates[0], 100, &shape).unwrap());
        for predicate in &predicates[1..] {
            assert!(matches_shape(predicate, 100, &shape).unwrap());
        }
        for bad in [
            "p=path:chr1:170-190,110-130",
            "p=path:chr1:",
            "p=overlap:chr1:1-2/0",
            "p=junction-near:chr1:1-2",
        ] {
            assert!(
                parse_predicates(&[bad.into()], &["chr1".into()]).is_err(),
                "{bad}"
            );
        }
    }

    #[test]
    fn overlap_endpoints_and_tolerance_have_explicit_boundaries() {
        let predicates = predicates(&[
            "a=overlap:chr1:105-115/5",
            "b=overlap:chr1:105-115/6",
            "s=start:chr1:100-101",
            "e=end:chr1:140-141",
            "notend=end:chr1:139-140",
            "n=junction-near:chr1:112-128/2",
            "notnear=junction-near:chr1:112-128/1",
        ]);
        let shape = Shape {
            blocks: vec![(0, 10), (30, 10)],
        };
        for (predicate, expected) in predicates
            .iter()
            .zip([true, false, true, true, false, true, false])
        {
            assert_eq!(matches_shape(predicate, 100, &shape).unwrap(), expected);
        }
        let p = &predicates[3];
        assert!(matches_shape(p, u32::MAX, &shape).is_err());
    }

    #[test]
    fn local_expressions_reject_unlinked_terminal_and_record_anchor() {
        let predicates = predicates(&["a=region:chr1:0-1", "t=terminal:chr1:0-1"]);
        assert!(validate_local_expression(
            &Expression::Predicate(0),
            &predicates,
            RegionMatchArg::Anchor
        )
        .is_err());
        assert!(validate_local_expression(
            &Expression::Predicate(1),
            &predicates,
            RegionMatchArg::AlignedBlock
        )
        .is_err());
    }

    #[test]
    fn compact_chain_preserves_junction_only_quantifier_truth() {
        let predicates = predicates(&["j=junction:chr1:110-130"]);
        let shapes = [
            Shape {
                blocks: vec![(0, 10), (30, 10)],
            },
            Shape {
                blocks: vec![(0, 5), (25, 10)],
            },
        ];
        let context = MatchContext {
            predicates: &predicates,
            shapes: &shapes,
            patterns: None,
            region_match: RegionMatchArg::AlignedBlock,
            placements: PlacementScopeArg::Unique,
        };
        let mut matcher = CompiledMatcher::new(&context);
        let record = molecule(&[(100, 0), (105, 1)], 3);
        assert_eq!(
            matcher
                .masks(
                    &record,
                    &[],
                    &Expression::Predicate(0),
                    MatchWithin::AllPlacements,
                    Engine::Compiled
                )
                .unwrap()
                .2,
            Some(TruthValue::True)
        );
        let absent = Expression::Not(Box::new(Expression::Predicate(0)));
        assert_eq!(
            matcher
                .masks(
                    &record,
                    &[],
                    &absent,
                    MatchWithin::AnyPlacement,
                    Engine::Compiled
                )
                .unwrap()
                .2,
            Some(TruthValue::False)
        );
    }

    #[test]
    fn alternative_placements_are_quantified_independently() {
        let predicates = predicates(&["a=region:chr1:100-110", "b=region:chr2:100-110"]);
        let shapes = [Shape {
            blocks: vec![(0, 10)],
        }];
        let patterns = [vec![
            PatAlt {
                chrom: 0,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            },
            PatAlt {
                chrom: 1,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            },
        ]];
        let context = MatchContext {
            predicates: &predicates,
            shapes: &shapes,
            patterns: Some(&patterns),
            region_match: RegionMatchArg::AlignedBlock,
            placements: PlacementScopeArg::All,
        };
        let mut matcher = CompiledMatcher::new(&context);
        let record = MolRec {
            chains: Default::default(),
            mms: smallvec::smallvec![(100, 0, 0, 1)],
            ..molecule(&[], 0)
        };
        let both = Expression::And(
            Box::new(Expression::Predicate(0)),
            Box::new(Expression::Predicate(1)),
        );
        for engine in [Engine::Scalar, Engine::Compiled] {
            let (mask, _, truth) = matcher
                .masks(&record, &[], &both, MatchWithin::AnyPlacement, engine)
                .unwrap();
            assert_eq!(mask, 3);
            assert_eq!(truth, Some(TruthValue::False));
            assert_eq!(
                matcher
                    .masks(
                        &record,
                        &[],
                        &Expression::Predicate(0),
                        MatchWithin::AnyPlacement,
                        engine
                    )
                    .unwrap()
                    .2,
                Some(TruthValue::True)
            );
            assert_eq!(
                matcher
                    .masks(
                        &record,
                        &[],
                        &Expression::Predicate(0),
                        MatchWithin::AllPlacements,
                        engine
                    )
                    .unwrap()
                    .2,
                Some(TruthValue::False)
            );
        }
    }

    #[test]
    fn expression_limits_prevent_unbounded_recursive_work() {
        let names = FxHashMap::from_iter([("a".to_owned(), 0)]);
        assert!(
            ExpressionParser::new(&format!("{}a", "!".repeat(129)), &names)
                .parse()
                .is_err()
        );
        assert!(
            ExpressionParser::new(&" ".repeat(MAX_EXPRESSION_BYTES + 1), &names)
                .parse()
                .is_err()
        );
    }
}
