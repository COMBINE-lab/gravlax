# Post-correction molecule BAM interchange

`aie export-molecule-bam` writes the archive's post-barcode-correction,
post-UMI-classing molecule abstraction to a sequence-free BAM. Its purpose is twofold:

1. provide an explicit, inspectable interchange form for Gravlax molecule evidence; and
2. measure BAM/CRAM container overhead against `.aie` after matching the retained function.

This is not an ordinary read-alignment BAM. The archive does not retain nucleotide UMI values,
and the archive used for the container comparison has more global UMI classes than can be embedded
as distinct 12-mers. The export therefore does not invent `UB` values. Generic SAM tools can read
and compress the file, but reconstructing the molecule abstraction requires the tags below.

## Placement records

Each chain representative or multimapper alternative is one mapped, sequence-free record. The
reference, 1-based start, strand flag, and `M`/`N` CIGAR encode its genomic placement. Records are
grouped by increasing `XM`; chain and multimapper group ids are independently dense from zero.

| Tag | Type | Meaning |
|---|---|---|
| `CB` | string | corrected 16-base cell barcode |
| `XC` | integer | dense cell id |
| `XI` | integer | opaque global UMI-class id |
| `XM` | integer | dense molecule id |
| `XW` | integer | signature read weight |
| `XK` | character | `C` for chain representative, `M` for multimapper alternative |
| `XG` | integer | group id within the record kind |
| `XA` | integer | alternative/representative index within the group |
| `XP` | integer | anchor flag; exactly one multimapper alternative has value 1 |
| `NH` | integer | 1 for a chain, number of alternatives for a multimapper group |

Non-anchor multimapper alternatives have the SAM secondary flag. All records for a molecule agree
on cell, UMI class, anchor chromosome, and anchor strand. A chain has one or two representatives.

## UMI adjacency records

After all placements, every retained undirected 1-mismatch edge is an unmapped record with
`XK:E`, smaller endpoint `XI`, and larger endpoint `XJ`. Edges are strictly sorted, unique, and
cell-scoped. This explicit suffix is necessary because class ids deliberately carry no nucleotide
sequence semantics.

## Round trip and validation

```sh
aie export-molecule-bam sample.aie --fai reference.fa.fai --out molecules.bam
aie replay-rows molecules.bam --from-molecule-bam --gtf genes.gtf \
  --barcodes barcodes.tsv --out-dir replay-from-bam
```

The importer rejects missing or inconsistent tags, nondense molecule/cell/class/group ids,
cross-cell or unsorted edges, nonempty sequence/qualities, unsupported CIGAR operations, and
invalid anchor/secondary/NH relationships. The CLI golden test exercises chain reduction,
multimapping, and a one-mismatch UMI edge and requires byte-identical replay matrices before and
after this interchange round trip.

For a storage comparison, convert the BAM to CRAM with a fixed CRAM version, reference digest,
`samtools` version, command line, and thread count recorded alongside the result. Label it
"post-correction molecule CRAM" rather than a generic BAM/CRAM baseline: its custom tags are
Gravlax-specific even though the container is standards-compliant.
