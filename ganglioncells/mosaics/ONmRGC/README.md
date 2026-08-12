# Pre-baked ON-midget-RGC circuits

This directory is the local gallery for pre-baked ON-mRGC circuit models. It
is empty in a fresh checkout.

The circuits ISETBio ships with are published in the `isetbio-mosaics`
Stanford Digital Repository deposit,
<https://purl.stanford.edu/vm447vf0975>, rather than bundled here. They are
downloaded on demand, verified against the manifest checksum, and cached
under `data/sdr/isetbio-mosaics/mrgc_on_circuits/`, which Git ignores. The
circuits are large, from 2.7 MB to 551 MB, so a first load can take a while.

`mRGCMosaic.loadPrebakedMosaic` searches this directory first, so a circuit
you synthesize locally is used in preference to the deposited copy, and
`mRGCMosaic.listPrebakedMosaics` lists both.

Circuits are selected by their historical ISETBio filename rather than by
the manifest's eccentricity and size fields, because staging dropped the
sign of the x eccentricity. See `sdrLegacyFilenameMap` for the map and the
evidence behind it.
