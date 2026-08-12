# ieWebGet and cache reference

## Existing ISETCam behavior

`ieWebGet` stores resource names, PURLs, and Stacks roots in the private
`urlDeposit` helper inside `../isetcam/utility/file/ieWebGet.m`. Its
`isetbio-mosaics` case accepts a caller-provided `downloaddir` and one
manifest-relative `depositfile`, preserving that path below the destination.
It also accepts `depositfile = 'all'` to mirror every asset named by the
top-level and collection manifests, verifying each MAT file's byte count and
SHA-256.

Keep URL construction and transport in ISETCam. ISETBio resolves a public
loader request through its collection manifest, then calls that ISETCam
capability.

## Required fetch sequence

1. Fetch the top-level manifest into
   `data/sdr/isetbio-mosaics/manifest.json`.
2. Fetch a collection manifest when needed and select one exact file record.
3. Construct the local path by appending the manifest's collection and `path`
   fields to `data/sdr/isetbio-mosaics/`.
4. Create parent directories. Download to a unique temporary sibling file.
5. Check byte count and SHA-256 against the manifest. Rename only after both
   checks succeed. Remove the temporary file on failure.
6. On an existing local file, recheck its SHA-256. A bad cache file must not
   be returned as a hit.

Use the Stacks download root plus the manifest-relative path as the remote
asset path only after confirming the deployed SDR hierarchy and URL encoding.

## Integration and verification

Register the resource after the PURL and download root are known. Keep the
PURL and Stacks root separate: the PURL is for people, while the Stacks root is
for `ieWebGet` transfers.

All four collections are migrated and share one fetch path, so these cases are
verified once rather than per loader. `test_sdrMosaicCacheFullOnly` covers a
first download, a cache hit, corrupt-cache replacement, missing records, and
the deposited lattices reproducing the golden crop values.

Still outstanding: the offline equivalents. Cache-hit and corrupt-cache tests
need either a network or a seam, because `sdrMosaicFetch` calls `ieWebGet`
directly and `sdrMosaicCacheRoot` is fixed. Adding a fake transport means
introducing an injection point; decide whether that testing seam is worth the
API surface before writing those tests. The offline suite currently covers
manifest-independent behavior only, through `test_sdrLegacyFilenameMap` and
`test_dataPaths`.

Do not make the repository's core test suite depend on the public SDR service.
Live checks are named `FullOnly` so they stay opt-in.
