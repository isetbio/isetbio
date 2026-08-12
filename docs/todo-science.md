# Science TODO

Open scientific questions and validation gaps. These are known, deliberately
left in place, and pinned by tests where possible so they cannot drift
unnoticed. Nothing here blocks use of ISETBio.

## RGCs with all-NaN surround weights

Some RGCs in the pre-baked ON-mRGC circuits have all-NaN surround weights, so
their surround sums to zero and they behave as center-only cells.

**Scope.** 15 of 102,619 RGCs, in 6 of the 16 circuits surveyed. RF center
matrices are clean everywhere. In the smallest circuit the affected cells are
RGCs 156, 264, 282, and 303 of 557.

**Consequence.** Under a uniform drive of 0.1, those four return exactly 0.1,
while the other 553 average 0.0146 — roughly seven times too strong, because
nothing subtracts a surround.

**Why it went unnoticed.** Computed responses contain no NaN at all. The
compute path silently drops NaN weights, so the affected cells look entirely
plausible and are simply missing their antagonism. Only assertions on values
expose it; a smoke test never would.

**Status.** Nicolas is not concerned by these cells, so they are left as they
are. The open question is *why* they arise during surround optimization —
whether it is a convergence failure for a small number of cells, or an expected
outcome of the fitting for cells in particular positions. Worth understanding
before the circuits are regenerated.

**Where it is pinned.**
`ganglioncells/_tests_/test_mRGCMosaicPrebakedGoldenValues.m`:

- `testSurroundConnectivityKnownDegenerateCells` asserts the exact set
- `testDegenerateCellsRespondWithoutSurround` demonstrates the consequence

Both are commented as documenting a known defect. If the circuits are ever
regenerated without these cells, both tests fail, which is the intended signal
to tighten them to require zero such cells. Regeneration would also mean
publishing a new SDR deposit version, so it is worth bundling with the manifest
fixes noted in [sdr-mosaic-data.md](sdr-mosaic-data.md).

## Missing golden coverage of the optics-dependent path

Nothing pins the optics-dependent mRGC path with golden values: visually
projected RF properties and spatial transfer functions. That is where the
model's scientific claims live, so it is the most valuable coverage still
missing.

It was left out deliberately on cost. Each case runs about 20 seconds even at
301 wavefront samples, because generating the native optics searches 29 defocus
values for the one maximizing the Strehl ratio. That belongs behind an opt-in
`FullOnly` test rather than in the core suite.

`ganglioncells/tests/xOLD/` holds the earlier PLOS2024 and in-vitro validation
scripts, which are the natural basis for STF-level regression.

A second, smaller gap: the pre-baked circuit tests pin a single small mosaic.
The NaN surround behavior above varies per circuit, so coverage across more of
them would be worth having.

## Watson A5/A6 round-trip error

Not a defect, and not to be "fixed" — Watson's Equations A5 and A6 are
separately fitted polynomials rather than an analytic inverse pair, so a
degrees → mm → degrees round trip is off by about 4% near the fovea. Both
implementations match the publication. Documented in
[rgc-watson.md](rgc-watson.md) and pinned by
`testWatsonConversionsAreNotExactInverses`. Listed here only so it is not
rediscovered as a bug.
