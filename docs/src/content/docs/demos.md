---
title: Reproducible demonstrations
description: Colab notebooks that execute Gravlax claims from an immutable, verified data capsule.
---

Three Google Colab notebooks exercise the preprint's central workflows. They
contain no saved scientific output and no fallback data address. Their
checked-in defaults pin the published `demo-data-v1`
[`demo-manifest.json`](https://github.com/COMBINE-lab/gravlax/releases/download/demo-data-v1/demo-manifest.json)
and its SHA-256,
`82c34aad442d478f1cb1243a6ccfe8ad9f937b81d9e1f946a8eb2cfc498214fd`. Each notebook verifies the
CLI, Python wheel, archives, annotations, and design files before executing
them. Both installed programs must report the exact version declared by the manifest. A collection
embeds canonical source paths,
so the multi-archive notebooks download each rooted source archive and build a
fresh collection inside Colab; they do not treat a detached `.aicollection` as
a portable data bundle.

- [Annotation reinterpretation](https://colab.research.google.com/github/COMBINE-lab/gravlax/blob/main/notebooks/01_annotation_reinterpretation.ipynb)
  compares two annotations against one rooted evidence archive and displays
  exact count deltas and a signed gene-level delta chart.
- [Multi-donor event discovery](https://colab.research.google.com/github/COMBINE-lab/gravlax/blob/main/notebooks/02_multi_donor_event_discovery.ipynb)
  reverse-searches a collection for recurrent junction and splice-event
  evidence, applies unique-chain raw-UMI-class sample, donor, group, and
  annotation-gap filters, and charts the strongest event's support across
  donor/group combinations.
- [Federated junction and co-occurrence](https://colab.research.google.com/github/COMBINE-lab/gravlax/blob/main/notebooks/03_federated_junction_cooccurrence.ipynb)
  queries one exact junction across a collection and drills into Boolean
  same-record region, junction, and terminal-tail predicates. It charts exact
  junction support by source archive and the resulting evidence patterns.

Every chart is generated from typed output produced during that notebook run.
The notebooks use a small deterministic SVG renderer built from the Python
standard library and Colab's IPython display API, so there is no plotting
package to resolve and no cached chart data.

The checked-in
[`demo-manifest.template.json`](https://github.com/COMBINE-lab/gravlax/blob/main/notebooks/demo-manifest.template.json)
lists every required asset and story parameter. It is deliberately incomplete:
null fields stop execution rather than silently selecting mutable or invented
data. A publishable manifest must validate against
[`demo-manifest.schema.json`](https://github.com/COMBINE-lab/gravlax/blob/main/notebooks/demo-manifest.schema.json),
use immutable release or repository-record URLs, and contain lowercase SHA-256
digests for every transport object. Archive assets should additionally declare
their rooted `aie-directory-root-v2` identities.

The notebooks are now one-click demonstrations pinned to the immutable
`demo-data-v1` capsule and the v0.1.5 software it declares. Clearing a locator,
using a mutable URL, or changing any downloaded byte still fails closed rather
than selecting fallback or cached results.
