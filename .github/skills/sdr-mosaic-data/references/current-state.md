# Current staging state

## Staging root

`/Volumes/Wandell/MATLAB/isetbio` is ready for the first SDR deposit. It
contains:

- `README.md`: plain-language guide for SDR visitors.
- `manifest.json`: top-level index.
- `cone_lattices/manifest.json`: 5 source lattices.
- `cone_mosaics/manifest.json`: 27 serialized local `cMosaic` objects.
- `mrgc_lattices/manifest.json`: 6 source lattices.
- `mrgc_on_circuits/manifest.json`: 21 pre-baked ON-mRGC circuit models.

All 59 MAT files have normalized names and a per-file byte count, SHA-256,
MATLAB format, eye, and collection-specific metadata. Every manifest checksum
was verified after creation.

`data_sabesan/` is deliberately outside the initial deposit. It contains
Professor Ram Sabesan's adaptive-optics measurements and its own README.

## Published deposit and remaining work

The first public release is available at
`https://purl.stanford.edu/vm447vf0975` (DRUID `vm447vf0975`). Its verified
Stacks download root is `https://stacks.stanford.edu/file/vm447vf0975`.
`manifest.json` and a nested cone-mosaic MAT file were both verified from this
root on 2026-08-10.

ISETCam now supports manifest-relative `isetbio-mosaics` asset download and a
verified full-deposit mirror. The manifest records the current ISETBio commit
as the manifest-generation provenance; it does not claim that commit generated
the individual historical data files.

All four collections are now migrated. Fetch, checksum verification, and cache
repair are shared by `sdrMosaicCacheRoot`, `sdrMosaicRecords`,
`sdrMosaicFetch`, and `sdrMosaicVerify` in `utility/sdr/`.

- Cone mosaics resolve through `mosaicLoad`. The 27 bundled files were removed
  in commit `757f1e26d`.
- Cone and mRGC lattices resolve through
  `retinalattice.import.sourceLatticeFile`, which both `finalConePositions`
  and `finalMRGCPositions` call. A locally generated lattice wins over the
  deposited copy.
- Pre-baked ON-mRGC circuits resolve through `sdrPrebakedCircuitFile`, called
  by `mRGCMosaic.loadPrebakedMosaic`.

The 5 cone lattices, 6 mRGC lattices, and 16 ON-mRGC circuits were removed
from the working tree after each file was verified byte-identical to its
deposited copy. The three gallery directories remain, each holding a README,
because lattice generation still writes into them.

## Eccentricity sign defect in mrgc_on_circuits

Staging dropped the sign of the x eccentricity from both the deposited path
and the `eccentricity_degrees` field. The circuits at x = -2 and x = +2
degrees are both `ecc-x2-y0` with `eccentricity_degrees` `[2, 0]`; they differ
only in size, so nothing currently collides, but a circuit cannot be selected
from published metadata alone.

Circuits are therefore selected by their historical ISETBio filename, using
`data/datafiles/sdr/isetbio-mosaics-legacy-filenames.json`. Every pairing in
that map was established by comparing SHA-256 checksums, and the cone-mosaic
entries were recovered from the tree of commit `757f1e26d^`. Five deposited
circuits never had a repository copy, so they carry an unverified candidate
name and cannot be loaded by name.

Fold `legacy_filename`, and a signed eccentricity for the circuits, into the
manifests when a later deposit version is cut. Add fields; do not rename
deposited paths, which would break published URLs.
`test_sdrLegacyFilenameMapFullOnly` fails once the sign is fixed, which is the
signal to revisit filename-keyed selection.

## Naming contract

Directories describe asset type. `left` or `right` gives eye identity.

- Lattices: `nominal-fov-Ndeg.mat`.
- Local cone mosaics: `ecc-xX-yY__size-wW-hH.mat`.
- ON-mRGC circuits: `ecc-xX-yY__size-wW-hH__optics-N__surround-N.mat`.

The detailed rationale, data model, metadata contract, and migration sequence
are in `local/MOSAIC-LOAD.md`.
