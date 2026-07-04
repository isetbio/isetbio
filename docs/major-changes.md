# Major Changes

This document summarizes changes affecting public APIs, compatibility, or major
implementation paths.

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
`validate/README.md`; the shared architecture is documented in ISETCam's
`docs/tutorial-example-test-architecture.md`.

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
