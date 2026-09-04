---
title: Reproducible demonstrations
description: Colab notebooks that execute Gravlax claims from an immutable, verified data capsule.
---

Three Google Colab notebooks exercise the preprint's central workflows. They
contain no saved scientific output and no fallback data address. Each requires
the HTTPS URL and SHA-256 of a published `gravlax.demo-capsule.v1` manifest,
then verifies the CLI, Python wheel, archives, annotations, and design files
before executing them. Both installed programs must report the exact version
declared by the manifest. A collection embeds canonical source paths,
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

The notebooks become one-click demonstrations only after the data capsule and
the release containing these command schemas have immutable public URLs. Until
then, they are executable, fail-closed launchers for configured capsules—not
simulated demonstrations.
