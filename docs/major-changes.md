# Major Changes

This document records repository changes that affect public APIs, compatibility,
or major implementation paths. It also lists significant deprecations that are
under consideration.

## 2026-06: Legacy Eye-Movement Struct API Removed

The legacy eye-movement struct API has been removed:

- `emCreate`
- `emGet`
- `emSet`

Current cone-mosaic eye movements are generated with the `fixationalEM` class.
Both `cMosaic.emGenSequence` and `coneMosaicRect.emGenSequence` use
`fixationalEM`.

The mosaic compute methods consume previously generated numeric eye-movement
paths:

- `cMosaic.compute` reads paths from the attached `fixEMobj`.
- `coneMosaicRect.compute` reads paths from `emPositions` or an explicit
  `emPath` argument.

Code that previously constructed or modified an eye-movement struct must migrate
to `fixationalEM`, a mosaic's `emGenSequence` method, or an explicitly supplied
numeric eye-movement path.

The historical analysis utility `emFitParameters` remains available and is
independent of the removed struct API.

### External Compatibility

The separate `isetvalidate` repository contains at least one validation that
uses the removed API:

`isetbioRDT/eyemovements/v_ibioRDT_eyeMovementsPhysio.m`

That validation must be migrated to `fixationalEM` before it can run against
this version of ISETBio.

## Rectangular Cone-Mosaic Window

`coneMosaicRect.window` and `coneRectWindow` use the App Designer
`coneRectWindow_App`. Recent maintenance improved the default launch behavior,
protected cone-absorption computation failures from closing or damaging the
window, and added GUI smoke tests.

This work does not introduce a new public API. It strengthens the existing
App Designer path in preparation for a possible legacy-window removal.

## Pending Decision: GUIDE Cone-Mosaic Window

The legacy GUIDE implementation remains in:

- `cones/rectangular/coneMosaicWindow.m`
- `cones/rectangular/coneMosaicWindow.fig`

Repository code no longer directly calls `coneMosaicWindow`; rectangular mosaic
windows route through `coneRectWindow_App`. Before removing the GUIDE files,
confirm that the App Designer window covers the required workflows and migrate
or document any remaining external callers.

