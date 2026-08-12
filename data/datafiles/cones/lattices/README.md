# Cone source lattices

This directory is the local cone lattice gallery. It is empty in a fresh
checkout.

The retina-wide cone lattices ISETBio ships with are published in the
`isetbio-mosaics` Stanford Digital Repository deposit,
<https://purl.stanford.edu/vm447vf0975>, rather than bundled here. They are
downloaded on demand, verified against the manifest checksum, and cached
under `data/sdr/isetbio-mosaics/cone_lattices/`, which Git ignores.

`retinalattice.import.sourceLatticeFile` searches this directory first, so a
lattice you generate locally with `retinalattice.generatePatch` is used in
preference to the deposited copy. Generation also writes its
`*_progress.mat` history here; those files are large and are neither
committed nor deposited.

See also: `sdrLegacyFilenameMap`, which pairs the historical filenames once
kept here with the deposited assets.
