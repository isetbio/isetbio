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

- Add the `isetbio-mosaics` resource to ISETCam, with the PURL for browsing
  and the verified Stacks download root for files.
- Extend `ieWebGet` with a narrow manifest-asset fetch path that accepts an
  SDR-relative file path and a local destination. Preserve existing resource
  behavior.
- Set the default ISETBio cache root to
  `fullfile(isetbioRootPath, 'data', 'sdr')`. `data/sdr/.gitignore` excludes
  downloaded files from Git. The resulting cache must be
  `<cache-root>/isetbio-mosaics/<collection>/<eye>/<file>.mat`.
- Fetch `manifest.json` before resolving an asset. Look up an exact manifest
  record; do not reconstruct a remote filename from MATLAB arguments.
- Download to a temporary file in the destination directory, verify its
  SHA-256 and byte count, then rename it into the mirrored cache path. For a
  cache hit, verify the checksum before returning it; replace a corrupt cache
  entry. Give missing-resource and checksum failures clear errors.

See `references/iewebget-and-cache.md` for the interface and validation
details.

## API and test boundaries

Keep `mosaicLoad`, cone/mRGC lattice importers, and
`mRGCMosaic.loadPrebakedMosaic` as public entry points. Place remote resolution
below them. Do not repurpose `cMosaic`'s `name` parameter, modify the staged
file layout casually, remove bundled assets, or rewrite Git history as part of
the initial loader work.

Keep core unit tests offline and deterministic. Add fake-transport tests for
manifest lookup, mirrored paths, first fetch, cache hit, corrupt cache, and
missing resource. Put a real SDR smoke test behind an opt-in `FullOnly` mode.
