---
name: retinal-mosaic-data
description: Work with ISETBio cone and midget-RGC mosaic data, including loading, generating, migrating, caching, or reviewing lattices, serialized cMosaics, pre-baked mRGCMosaics, and compute-ready mosaics. Use when a task concerns files under data/datafiles/cones, data/datafiles/rgc/lattices, ganglioncells/mosaics, mosaic-loading APIs, or cone-to-RGC connectivity.
---

# Retinal Mosaic Data

## Overview

Distinguish geometric source lattices from completed retinal circuit
models. Preserve the established `cMosaic` and `mRGCMosaic` APIs while
using “cone-to-RGC retinal circuit” for a model that includes both cells and
their connectivity.

## Identify the Data Type Before Editing a Loader

| Data | Location | Contents | Role |
| --- | --- | --- | --- |
| Cone lattice | SDR `cone_lattices/`; local gallery `data/datafiles/cones/lattices/` | Retina-wide cone RF-center positions | Crop geometry used to construct a `cMosaic`. Resolved by `retinalattice.import.sourceLatticeFile`. |
| RGC lattice | SDR `mrgc_lattices/`; local gallery `data/datafiles/rgc/lattices/` | Retina-wide midget-RGC RF-center positions | Crop destination geometry used to generate an RGC mosaic. Same resolver. |
| Serialized cone mosaic | SDR `cone_mosaics/` | A local, complete `cMosaic` | Load through `mosaicLoad` for a precomputed cone realization. |
| Pre-baked RGC mosaic | SDR `mrgc_on_circuits/`; local gallery `ganglioncells/mosaics/ONmRGC/` | A local `mRGCMosaic`, its input `cMosaic`, and computed connectivity/model data | Load through `mRGCMosaic.loadPrebakedMosaic`, which resolves via `sdrPrebakedCircuitFile`. |
| Compute-ready RGC mosaic | External resource directory selected by `getpref('isetbio').rgcResources` | A further prepared `mRGCMosaic` used by optimization workflows | Load through `loadComputeReadyRGCMosaic`. Not part of the SDR deposit. |

The four SDR collections are cached under `data/sdr/isetbio-mosaics/`, which
Git ignores. The local galleries above are empty in a fresh checkout; a file
generated locally into one of them takes precedence over the deposited copy.
See `sdr-mosaic-data` for the deposit, cache, and legacy-filename contract, and
`docs/sdr-mosaic-data.md` for the user-facing account.

Do not confuse “mosaic” in a filename with a completed circuit. A lattice
contains positions, not cone types, cone–RGC wiring, receptive-field weights,
or response gains. A pre-baked `mRGCMosaic` contains these circuit elements
and its input cone mosaic.

## Follow the Actual Construction Path

Use this dependency order when reasoning about synthesis:

```text
cone lattice → cMosaic
RGC lattice  ─┘
                → cone-to-mRGC connectivity → mRGCMosaic circuit model
```

Key implementation points:

- `cMosaic` normally crops a cone lattice in
  `retinalattice.import.finalConePositions`.
- `retinalattice.import.finalMRGCPositions` crops an RGC lattice.
- `RGCMosaicConstructor.compute.componentsForRFcenterConnectivity` accepts a
  supplied `cMosaic`, builds cone-to-RGC center connectivity with
  `coneToMidgetRGCConnector`, and returns the components used to instantiate
  `mRGCMosaic('withComponents', components)`.
- `RGCMosaicConstructor.compute.centerConnectedMosaic` is the usual
  high-level pipeline, but it creates its own compatible input cone mosaic;
  it is not a general public wrapper for a caller-supplied `cMosaic`.
- `mRGCMosaic` is the established class name. Do not rename it merely because
  its instances model a circuit.

When building a workflow around an arbitrary `cMosaic`, validate compatible
eye, retinal coordinate transforms, eccentricity/extent, and extra cone
support for RGC surrounds. Do not reuse precomputed connectivity after
changing the input cone positions or cone types.

## Choose the Smallest Correct Loader Change

- Use `mosaicLoad` for the serialized cone-mosaic library. Its legacy
  parameter and filename calls select a manifest record and fetch it on
  demand; use `mosaicLoad('list')` to discover available cone mosaics. Do
  not overload the `cMosaic` constructor’s existing `name` parameter: it is
  a human-readable label, not a resource identifier.
- Use `mRGCMosaic.loadPrebakedMosaic` for the pre-baked ON-mRGC circuit
  models. Preserve its parameter-to-filename convention and its cropping
  behavior. The filename it builds is the lookup key into the deposit,
  because the published circuit metadata lost the sign of the x
  eccentricity; do not switch to selecting by eccentricity and size.
- Use `mRGCMosaic.loadComputeReadyRGCMosaic` only for the configured external
  compute-ready resource set; it currently relies on ISETBio preferences.
- Keep the lattices available for any workflow that regenerates a mosaic. A
  pre-baked file cannot replace them as a general source of geometry.

For remote or on-demand data work, use ISETCam’s `ieWebGet` rather than adding
another downloader. Keep public loaders and their cache paths stable; resolve
or fetch an absent resource beneath them. Require a stable deposit path,
filename, byte size, checksum, and data-format/version metadata before
removing a bundled asset.

## Inspect and Validate

- Search all call sites before relocating data or changing a loader:
  `rg -n "mosaicLoad|loadPrebakedMosaic|loadComputeReadyRGCMosaic|finalConePositions|finalMRGCPositions"`.
- Inspect representative MAT-file variables before writing migration logic;
  do not infer content from filenames.
- Update data-path tests if their contract deliberately changes. All four
  collections have been migrated, so a test must not assert that a bundled
  lattice or circuit exists on disk.
- Validate both a first fetch/load and a cache-hit load, plus a clear error
  when the resource is unavailable. For circuit generation, verify dimensions
  and cone/RGC counts of the connectivity matrices as well as successful
  construction.
