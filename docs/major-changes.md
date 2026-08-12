# Major Changes

This document summarizes changes affecting public APIs, compatibility, or major
implementation paths.

## 2026-08: Retinal Mosaic Data Moved to the Stanford Digital Repository

The large retinal data files are no longer stored in the repository. They are
published in the SDR deposit `isetbio-mosaics`,
<https://purl.stanford.edu/vm447vf0975>, and downloaded on demand into an
untracked cache under `data/sdr/isetbio-mosaics/`. All four collections have
migrated: 5 cone lattices, 6 midget-RGC lattices, 27 serialized cone mosaics,
and 21 pre-baked ON-mRGC circuits, about 1.2 GB in total.

Public entry points are unchanged. `mosaicLoad`,
`retinalattice.import.finalConePositions`,
`retinalattice.import.finalMRGCPositions`, and
`mRGCMosaic.loadPrebakedMosaic` take the same arguments and return the same
results; remote resolution sits below them. Every download is verified against
the manifest's SHA-256, and a corrupt cache entry is replaced rather than
loaded.

Behavior that did change:

- A first call for a given file downloads it. Sizes run from 200 KB to 551 MB.
- `data/datafiles/cones/lattices/`, `data/datafiles/rgc/lattices/`, and
  `ganglioncells/mosaics/ONmRGC/` are empty in a fresh checkout. A file you
  generate locally into one of them still takes precedence over the deposited
  copy.
- `mRGCMosaic.listPrebakedMosaics` and
  `retinalattice.listPrecomputedPatches` list what the deposit publishes as
  well as what is present locally, rather than listing a directory.
- `mRGCMosaic.loadPrebakedMosaic` returns the directory that actually holds
  the loaded file, which is the cache directory for a downloaded circuit. It
  no longer necessarily joins with the returned filename.

Usage, caching, and troubleshooting are documented in
[sdr-mosaic-data.md](sdr-mosaic-data.md).

## 2026-06: Validation Infrastructure Reorganized

ISETBio and ISETCam now keep their public test runners and validation
infrastructure in top-level `validate/` directories. Runner names are singular
and consistent:

- ISETBio: `isetbioUnitTest`, `isetbioTutorialTest`, `isetbioExampleTest`
- ISETCam: `ieUnitTest`, `ieTutorialTest`, `ieExampleTest`

The tutorial and example runners use `'selection'` to run one script and
`'start'` to run from one script through the remainder of the deterministic
plan. The earlier plural names, positional selector, and `'select'` aliases
were removed during this early-development reorganization.

ISETCam owns the shared `ieRunTutorialExampleTests` engine and `ieTestReport`.
The engine provides cross-repository discovery, isolation, skipping, durable
checkpoints, and consistent reporting. Usage is documented in
[testing.md](testing.md); the shared architecture is documented in ISETCam's
[tutorial/example test architecture](../../isetcam/docs/tutorial-example-test-architecture.md).

## 2026-06: RemoteDataToolbox References Retired

RemoteDataToolbox is deprecated. ISETBio removed its RDT configuration and
obsolete cone-mosaic artifact-publishing examples. The ToolboxToolbox local
hook still configures local validation data and retains inert RDT preference
fields currently required by UnitTestToolbox.

Compatibility references to current `isetvalidate` directory and listing names
remain temporarily while independently maintained repositories migrate.

## 2026-06: Bundled Data Consolidated

The former `dataiset/` and `data/datafiles/` trees were merged into:

`data/datafiles/`

Use `isetbioDataPath` for ISETBio-owned data and `isetRootPath` for ISETCam-owned
data. The move also corrected `isetbioRootPath`, which still reflected its old
source location. Tests cover representative paths and golden values from major
data collections.

## 2026-06: Legacy Eye-Movement Struct API Removed

The `emCreate`, `emGet`, and `emSet` struct API was removed. Cone-mosaic eye
movements now use `fixationalEM`; both `cMosaic.emGenSequence` and
`coneMosaicRect.emGenSequence` use this class.

Mosaic compute methods consume generated numeric paths from the attached
`fixEMobj`, `emPositions`, or an explicit `emPath`. The independent analysis
utility `emFitParameters` remains available.

The external `isetvalidate` script
`isetbioRDT/eyemovements/v_ibioRDT_eyeMovementsPhysio.m` must migrate to
`fixationalEM`.

## Rectangular Cone-Mosaic Window

`coneMosaicRect.window` and `coneRectWindow` use the App Designer
`coneRectWindow_App`. Maintenance improved launch behavior, protected the
window from absorption-computation failures, and added GUI smoke tests.

The unused GUIDE implementation remains in
`cones/rectangular/coneMosaicWindow.m` and `.fig` pending confirmation that the
App Designer window covers all required workflows.
