# ieWebGet and cache reference

## Existing ISETCam behavior

`ieWebGet` stores resource names, PURLs, and Stacks roots in the private
`urlDeposit` helper inside `../isetcam/utility/file/ieWebGet.m`. It downloads
only known resource types through switch cases. It already supports a caller
provided `downloaddir`, but it does not currently provide a general,
manifest-relative single-file fetch API.

Add the SDR resource and a minimal generic asset case in ISETCam. Keep URL
construction and transport there. Let ISETBio resolve its resource identifier
through its collection manifest and call that ISETCam capability.

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

Validate these cases before changing public loaders:

- first manifest and asset download;
- cache hit with no network request;
- corrupt cached asset replacement;
- missing manifest record;
- failed checksum or byte-count validation; and
- one opt-in live SDR download.

Use a temporary cache plus fake transport for ordinary unit tests. Do not make
the repository's core test suite depend on the public SDR service.
