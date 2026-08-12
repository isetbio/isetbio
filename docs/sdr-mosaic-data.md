# Retinal mosaic data and the Stanford Digital Repository

ISETBio's large retinal data files are no longer stored in the Git repository.
They are published in the Stanford Digital Repository (SDR) and downloaded on
demand the first time your code needs one.

- Deposit: **isetbio-mosaics**, <https://purl.stanford.edu/vm447vf0975>
- 59 files, about 1.2 GB in total

For most people this changes nothing. Loading a mosaic works as it always did;
the first call for a given file pauses to download it, and later calls read a
verified local copy.

## Why

The deposit holds roughly 1.2 GB of lattices, cone mosaics, and pre-baked
retinal circuits. Keeping them in Git meant every clone paid for all of them,
whether or not a given project used any. It also meant the data had no citable
identity of its own. The SDR gives the data a permanent public URL, a version
history, and published checksums, and it lets you download only what you use.

## Using it

Nothing special is required. Load data the way you always have.

```matlab
% A precomputed cone mosaic: 0.5 x 0.5 degrees at the fovea
cm = mosaicLoad([0.5 0.5], [0 0]);

% See what the deposit publishes
mosaicLoad('list');                  % cone mosaics
mosaicLoad('list', 'cone_lattices'); % or any other collection
```

Source lattices are fetched automatically when you build a mosaic that needs
one, so there is no separate step:

```matlab
cm = cMosaic('sizeDegs', [1 1], 'eccentricityDegs', [2 0]);
```

Pre-baked ON-midget-RGC circuits load by their usual specifiers:

```matlab
mRGCMosaic.listPrebakedMosaics();     % what is available
[theMRGCMosaic, theOI] = mRGCMosaic.loadPrebakedMosaic( ...
    'human', 'JCNpaperDefaultSubject', 'JCNpaperNasal2DegsTinyMosaic', 'default');
```

### What to expect on a first call

The first use of a file downloads it. Sizes vary a great deal:

| Collection | Files | Typical size | What it is |
|---|---|---|---|
| `cone_lattices` | 5 | 18-21 MB | Retina-wide cone RF-center positions, cropped to build a `cMosaic` |
| `mrgc_lattices` | 6 | ~2.9 MB | Retina-wide midget-RGC RF-center positions |
| `cone_mosaics` | 27 | 0.2-30 MB | Serialized `cMosaic` objects at various sizes and eccentricities |
| `mrgc_on_circuits` | 21 | 2.7-551 MB | Complete ON-midget-RGC circuit models with connectivity |

A small cone mosaic arrives in a second or two. The largest circuit is 551 MB,
so plan for a wait the first time you load one.

## The local cache

Downloads are cached under:

```
isetbio/data/sdr/isetbio-mosaics/<collection>/<eye>/<file>.mat
```

The directory mirrors the deposit, so you can see what you have by looking at
it. Git ignores it, so cached files are never committed.

Every file is checked against the SHA-256 checksum published in the deposit
manifest, both when downloaded and on every later read. A truncated or
corrupted cache entry is detected and replaced automatically rather than
silently loaded. To force a fresh download, delete the file, or the whole
`data/sdr` directory.

### Downloading everything at once

If you are heading offline, or want to prime a shared machine, mirror the whole
deposit:

```matlab
ieWebGet('deposit name', 'isetbio-mosaics', ...
         'deposit file', 'all', ...
         'downloaddir', fullfile(isetbioRootPath, 'data', 'sdr', 'isetbio-mosaics'));
```

This downloads all 59 files, about 1.2 GB, verifying each one.

## Data you generate yourself

A file you generate locally always takes precedence over the deposited copy.
The loaders look in the repository directories first:

- `data/datafiles/cones/lattices/` — cone source lattices
- `data/datafiles/rgc/lattices/` — midget-RGC source lattices
- `ganglioncells/mosaics/ONmRGC/` — pre-baked ON-mRGC circuits

These are empty in a fresh checkout, and each holds a README explaining its
role. If you run `retinalattice.generatePatch` or synthesize your own circuit,
the result lands in one of them and is used in preference to anything
downloaded. Lattice generation progress files (`*_progress.mat`) are local
artifacts; they are neither committed nor deposited.

## Troubleshooting

**A download fails or the checksum does not match.** The loader raises an error
rather than returning a suspect file, and removes the partial download. Re-run
the call. If it fails repeatedly, confirm you can reach
<https://stacks.stanford.edu/file/vm447vf0975/manifest.json>.

**"No SDR cone mosaic has size ... at eccentricity ..."** The deposit has no
file with exactly those parameters. Run `mosaicLoad('list')` to see what
exists; sizes and eccentricities must match a published record exactly.

**A circuit will not load by name.** Five deposited circuits were contributed
without a matching repository copy, so their historical filename could not be
established. `mRGCMosaic.listPrebakedMosaics` marks them
`(candidate name, unverified)`. They can be downloaded directly by their
deposit path but not selected by name.

**Tests that need the network.** Ordinary test runs stay offline. Tests whose
names end in `FullOnly` reach the live deposit and are excluded from the
default suite.

## Known characteristic of the deposited circuits

A small number of RGCs in the `mrgc_on_circuits` files have all-NaN surround
weights and behave as center-only cells: 15 of 102,619 RGCs, in 6 of the 16
circuits surveyed. They respond about seven times too strongly because nothing
subtracts a surround. The behavior is pinned by tests and is described, with
the open question behind it, in [todo-science.md](todo-science.md).

## Notes for maintainers

The mapping between each deposited file and the filename ISETBio used before
the migration is recorded in
`data/datafiles/sdr/isetbio-mosaics-legacy-filenames.json`. Every pairing was
established by comparing SHA-256 checksums against a repository copy, so the
map is evidence rather than a naming convention.

That map is not merely historical. Staging dropped the sign of the x
eccentricity for the `mrgc_on_circuits` collection: the circuits at x = -2 and
x = +2 degrees both appear as `ecc-x2-y0` with `eccentricity_degrees` `[2 0]`.
A circuit therefore cannot be selected from the published metadata alone, and
`sdrPrebakedCircuitFile` selects by legacy filename instead. A future deposit
version should add `legacy_filename` to every record and a signed eccentricity
to the circuits — adding fields, never renaming deposited paths, which would
break published URLs.

The download and verification layer lives in `utility/sdr/`:
`sdrMosaicCacheRoot`, `sdrMosaicRecords`, `sdrMosaicFetch`, `sdrMosaicVerify`,
`sdrLegacyFilenameMap`, and `sdrPrebakedCircuitFile`. HTTP transport stays in
ISETCam's `ieWebGet`; ISETBio resolves a request through the manifest and then
calls it.

## See also

- [major-changes.md](major-changes.md) — the compatibility record for this migration
- ISETCam `ieWebGet` — transport and the deposit registry
