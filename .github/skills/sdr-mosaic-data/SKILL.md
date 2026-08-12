---
name: sdr-mosaic-data
description: Create, publish, organize, or fetch the ISETBio retinal mosaic data deposit in the Stanford Digital Repository (SDR). Use for SDR collection/deposit setup, PURL or Stacks download-root capture, staged manifests, ISETCam ieWebGet registry and download support, or the ISETBio manifest-aware untracked mosaic-data cache.
---

# SDR Mosaic Data

Use this skill for the ISETBio `isetbio-mosaics` SDR workflow. Read
`references/current-state.md` first. Then read `local/MOSAIC-LOAD.md` when it
is available; it is the detailed design record and remains authoritative for
the loader contract and deferred decisions.

## Published deposit

`isetbio-mosaics` is published at
`https://purl.stanford.edu/vm447vf0975` (DRUID `vm447vf0975`). Its verified
Digital Stacks download root is
`https://stacks.stanford.edu/file/vm447vf0975`. Both `manifest.json` and a
nested mosaic asset have been verified at that root. Keep the deposited
directory hierarchy intact when resolving relative manifest paths.

## Deposit workflow

1. Inspect the staged root and every manifest before uploading. Preserve the
   four collection directories and their nested `left`/`right` paths. Do not
   include `data_sabesan` unless its separate collection and metadata contract
   have been approved.
2. Create the SDR collection/deposit named `isetbio-mosaics`, retaining the
   hierarchy. Use the top-level `README.md` for the plain-language file guide
   and the SDR description for the short catalog-level summary.
3. The first public release is PURL `https://purl.stanford.edu/vm447vf0975`,
   DRUID `vm447vf0975`, with Stacks root
   `https://stacks.stanford.edu/file/vm447vf0975`. For later versions, record
   the public version and re-test `manifest.json` plus one nested MAT-file URL.
   Do not derive a download root from a guessed URL pattern.
4. Compare uploaded file sizes/checksums with the staged manifests before
   treating the deposit as ready for loader work.

## Fetch and cache design

Use ISETCam's `ieWebGet`; do not add a second HTTP downloader in ISETBio.
The current implementation is `../isetcam/utility/file/ieWebGet.m`. Its SDR
registry is private (`urlDeposit`) and its cases are resource-specific, so
adding a registry row alone is not enough for arbitrary manifest paths.

- The `isetbio-mosaics` ISETCam resource records the PURL for browsing and
  the verified Stacks root for files. `ieWebGet` accepts one SDR-relative
  asset path plus a caller-supplied destination, and `deposit file = 'all'`
  mirrors every manifest-listed asset.
- Set the default ISETBio cache root to
  `fullfile(isetbioRootPath, 'data', 'sdr')`. `data/sdr/.gitignore` excludes
  downloaded files from Git. The resulting cache must be
  `<cache-root>/isetbio-mosaics/<collection>/<eye>/<file>.mat`.
- `mosaicLoad` fetches `manifest.json` before resolving a cone mosaic, looks
  up an exact record, and does not reconstruct a remote filename from MATLAB
  arguments. It downloads to a temporary sibling, verifies SHA-256 and byte
  count, then renames into the mirrored cache path. It verifies cache hits
  and replaces corrupt entries.

See `references/iewebget-and-cache.md` for the interface and validation
details.

## API and test boundaries

Keep `mosaicLoad`, cone/mRGC lattice importers, and
`mRGCMosaic.loadPrebakedMosaic` as public entry points. Place remote resolution
below them. Do not repurpose `cMosaic`'s `name` parameter, modify the staged
file layout casually, or rewrite Git history.

All four collections are migrated. The shared fetch layer is in `utility/sdr/`:
`sdrMosaicCacheRoot`, `sdrMosaicRecords`, `sdrMosaicFetch`, `sdrMosaicVerify`,
`sdrLegacyFilenameMap`, and `sdrPrebakedCircuitFile`. Lattices resolve through
`retinalattice.import.sourceLatticeFile`. A locally generated file always wins
over the deposited copy, so the three gallery directories remain, each with a
README. Add a new collection by resolving through these helpers rather than
writing a second fetch path.

Circuits are selected by their historical ISETBio filename, not by manifest
metadata. Staging dropped the sign of the x eccentricity for
`mrgc_on_circuits`, so x = -2 and x = +2 both read `[2 0]`. The legacy map at
`data/datafiles/sdr/isetbio-mosaics-legacy-filenames.json` is the lookup table;
every pairing in it was established by SHA-256 comparison. Before deleting any
bundled asset, confirm it is byte-identical to its deposited copy.

Keep core unit tests offline and deterministic. Live checks are named
`FullOnly`: `test_sdrMosaicCacheFullOnly` covers first fetch, cache hit,
corrupt-cache replacement, missing records, and golden crop parity;
`test_sdrLegacyFilenameMapFullOnly` checks the map against the manifests and
fails once the eccentricity sign is fixed. Offline fake-transport tests are
still outstanding and need an injection seam, since `sdrMosaicFetch` calls
`ieWebGet` directly and `sdrMosaicCacheRoot` is fixed.

## Documentation follow-up

`docs/sdr-mosaic-data.md` is the user-facing account of the deposit: loading,
the cache, full-deposit mirroring, locally generated data, and troubleshooting.
`docs/major-changes.md` carries the compatibility entry. Update both when
loader behavior changes.

The wiki (`isetbio/wiki`) is still to be written, and `docs/sdr-mosaic-data.md`
is its intended source. The repository `docs/` directory may later become a set
of pointers to wiki pages; do not restructure it as part of loader work.

## Deferred Git-history cleanup

Removing a migrated asset from the current tree is a normal commit and makes
future repository checkouts smaller only when paired with a later history
rewrite. Do not rewrite history during the initial SDR migration. Revisit that
separate project after roughly one to two years, once the SDR loader and
deposit have proven stable and downstream users have migrated.

At that time, preserve the SDR PURL, manifests, checksums, and a verified bare
repository backup before using a scoped `git filter-repo` rewrite to remove
only approved migrated paths. Inventory release tags, branches, forks, and
downstream repositories; announce the force-push/re-clone cutover; then verify
a fresh clone, the affected releases, and representative SDR loads. Never
rewrite history merely to reduce repository size without that coordination.
