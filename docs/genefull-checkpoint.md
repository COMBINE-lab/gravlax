# GeneFull replay: unfinished handoff checkpoint

This branch is a work in progress, not a release candidate. No version bump,
release tag, or publication has been performed.

## Location and base

- Branch: `codex/genefull-replay`
- Base: `ee00b23` (main, version 0.2.1)
- Source checkout: `/tmp/gravlax-v016` on `nomad01.umiacs.umd.edu`
- The biological data belong on nomad00. The current task was executing on
  nomad01 even though its configured project path names the nomad00 project.
- Task ID: `01a053f1-8244-7820-b21d-bee1104f645f`. An older entry with the same
  ID appears on nomad00; preserve the latest conversation when handing off.

## Implemented so far

- `anno::assign::GeneFullIndex`: builds chromosome-specific gene-span indexes
  from annotation transcript/exon bounds. Tests aligned blocks against those
  spans, including introns, with the selected strand relationship. Genes inside
  a skipped alignment intron are not counted solely from the outer read span.
- Candidate assignment accepts the optional GeneFull index for primary and
  alternative placements, feeding the existing global UMI-collapse machinery.
- `replay-rows --gene-full` selects this path for streaming and eager replay,
  including the shared BAM-input path. The flag conflicts with `--velocity`
  and `--audit-multigene`; default Gene behavior is intended to remain unchanged.
- Uniform-report parameters record `gene_full` and `counting_model`.
- Archive and compiled-annotation encodings are unchanged.

## Checks performed

- `cargo check -p gravlax --locked` succeeded after the implementation edits.
- It reported an unused `replay_rows_stranded` import in
  `crates/aie/src/archivecmd.rs`; cleanup is still pending.
- `git diff --check` passed before this checkpoint.
- No new unit, integration, STARsolo comparison, or biological benchmark has
  been run. Compilation is not evidence of scientific correctness.

## Reference semantics and unfinished validation

STAR 2.7.11b source was fetched to `/tmp/STAR-genefull-reference` on nomad01,
commit `b1edc1208d91a53bf40ebae8669f71d50b994851`:

- `source/Transcriptome.cpp` derives each gene's full bounds from exon records.
- `source/Transcriptome_geneFullAlignOverlap.cpp` unions genes whose spans
  overlap any aligned block on the accepted strand, without requiring
  transcript-junction concordance.
- `source/SoloReadFeature_record.cpp` consumes the resulting gene sets.

Next steps:

1. Review annotation edge cases, including genes represented on multiple
   chromosomes/strands, half-open boundaries, disjoint isoforms, nested genes,
   multi-gene and multi-placement evidence, and unusual splice junctions.
2. Test candidate sets against a simple independent overlap implementation and
   synthetic STARsolo GeneFull runs. Test strand options and UMI-collapse cases.
3. Verify streaming/eager, archive/BAM, and GTF/compiled-annotation equivalence;
   ensure default Gene results remain unchanged and incompatible flags fail.
4. Validate report/metadata model identity. Python convenience support and user
   documentation remain unimplemented. Run the full applicable Rust/Python
   tests, formatting checks, and lint checks.
5. Audit assignment-statistic units: the existing accumulator increments its
   assigned counter per representative row, whereas its denominator is the
   number of molecule records. The historical 34% ratio is therefore not yet
   validated as a fraction of distinct archived molecules. Report unique
   molecule/class and final UMI quantities with explicit denominators.
6. Run the matched brain-nucleus analysis below before interpreting biological
   consequences or preparing the next release.

## Biological rerun to complete on nomad00

The original study root is
`/scratch5/rob/recomb-2027-annotation-independent-evidence`; the user notes
the underlying storage is `/mnt/scratch5` on nomad00. Verify actual paths and
repository status before using them. Read applicable repository instructions
and preserve existing work.

Known relative inputs from the saved analysis scripts:

- `runs/archive/d2p.aie`
- `annotations/gencode.v49.annotation.gtf`
- `runs/oracle/d2p/v49/Solo.out/Gene/raw/`
- `runs/oracle/d2p/v49/Solo.out/Gene/filtered/barcodes.tsv`
- `runs/oracle/d2p/v49/Solo.out/GeneFull/raw/`

The original STARsolo run requested both Gene and GeneFull, but
`scripts/97_d2p_campaign.sh` selected Gene output and Gene-called nuclei for
the replay comparison. Its log reports 130,270,701 molecule records,
44,325,139 assigned representatives, and 28,630,629 collapsed UMIs across the
requested raw barcode output. These are not all counts restricted to the
6,460 called nuclei. Saved scripts/logs are available on nomad01 in
`/tmp/gravlax-paper-scripts-postv1` if needed for recovery.

First compare Gene and GeneFull on the same called nuclei, using matched
STARsolo reference outputs. Then separately assess nucleus calling under the
intron-inclusive model. Reevaluate replay fidelity, annotation-change effects,
clustering/markers, and annotation-extension gains as appropriate. Preserve
original results and write new outputs to separate directories. Do not label
GeneFull results as measured until the runs and comparisons have completed.

The manuscript checkout on nomad01 is
`/tmp/moleceular-evidence-store-paper-audit`; its latest editorial changes were
pushed as `c897a26` on `recomb-refinement`. No GeneFull numerical results have
been substituted into the manuscript.
