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

An ISETCam registry entry exists for browsing and for recording the verified
download root. Manifest-aware asset download and the remote ISETBio loader are
still deferred. The manifest records the current ISETBio commit as the
manifest-generation provenance; it does not claim that commit generated the
individual historical data files.

## Naming contract

Directories describe asset type. `left` or `right` gives eye identity.

- Lattices: `nominal-fov-Ndeg.mat`.
- Local cone mosaics: `ecc-xX-yY__size-wW-hH.mat`.
- ON-mRGC circuits: `ecc-xX-yY__size-wW-hH__optics-N__surround-N.mat`.

The detailed rationale, data model, metadata contract, and migration sequence
are in `local/MOSAIC-LOAD.md`.
