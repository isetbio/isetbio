
## **From Claude, July 28, 2027**

I reviewed the ganglion cell tests and added quantified ones — and doing so surfaced a defect in the shipped model you should know about before release.

## The existing tests weren't validating the model

`ganglioncells/_tests_` contained exactly one test file. It drives `optimizedStrehlRatioPSF` against `FakeConeMosaicForStrehlTest`, a stub returning a synthetic 3×3×3 PSF. Its golden values (0.375, 0.4375, 0.5) verify that the search picks the maximum — useful, but it's a test of arithmetic. Nothing asserted anything about the analytic formulae the mosaics are built from, or about the prebaked mosaics scientists will actually load.

## A real defect in the shipped mosaics

**Four of the 557 RGCs in the small mosaic (indices 156, 264, 282, 303) have all-NaN surround weights.** Their surround sums to zero, so they are center-only cells.

This is not cosmetic. Under a uniform drive of 0.1, those four return exactly **0.1** — the raw input, unattenuated — while the other 553 average  **0.0146** . They respond roughly seven times too strongly because nothing is subtracting a surround.

Across all 16 shipped mosaics it affects  **15 of 102,619 RGCs in 6 of the 16 files** . The RF center matrices are clean everywhere. The reason this went unnoticed is that computed responses contain no NaN at all — the compute path silently drops the NaN weights, so the cells look perfectly plausible and are simply missing their antagonism. A smoke test would never catch this; it only shows up when you assert on values.

I pinned it in two tests rather than papering over it: one asserts the exact set of degenerate cells, the other demonstrates the behavioral consequence. Both are commented as documenting a known defect, to be tightened to require zero such cells once the mosaics are regenerated.

## A second finding worth knowing

Watson's Equations A5 and A6 are separately fitted polynomials, not an analytic inverse pair. Their linear terms disagree: `1/0.268 = 3.731` versus the `3.556` used in A6. A degs→mm→degs round trip is therefore off by  **−4.4% at 0.5°, −4.2% at 1°, −3.7% at 2°** , converging only past 20°. That is precisely the eccentricity band most prebaked mosaics occupy, and 14 files call both conversions. The implementations match the publication, so I pinned the behavior with a test that documents it rather than "fixing" either function.

## What the tests cover

`test_RGCmodelsGoldenValues` pins Watson A5/A6 and the Croner & Kaplan relationships — the 6.7 surround/center Rc ratio, surround radius versus eccentricity, and integrated sensitivity ratio. I checked these against the literature rather than just recording output: 0.2683 mm/deg and the 6.7 ratio are the published values.

`test_mRGCMosaicPrebakedGoldenValues` runs against the real shipped mosaic (loads in 0.26s) and pins geometry, cone composition (L:M of 2270:1118 = 2.03, S at 10%), RF center connectivity (697 unit-valued connections, 1–2 cones per center, mean 1.2513 — the expected midget arrangement at 2°), and a deterministic noise-free response. The response drive is an analytic function of cone position, not a random draw, so expected values don't depend on the RNG implementation across MATLAB versions; the drive is separately checksummed so a change there fails distinctly from a change in compute.

Every golden number is paired with a qualitative assertion — surround radius must increase with eccentricity, sensitivity ratio must stay in (0,1), mean response to a zero-mean grating must be small relative to its spread — so a failure tells you what broke, not just that a digit moved.

Both runners auto-discover the files, so no runner changes were needed. `isetbioUnitTest` goes **84 → 96** passing, `ganglioncellsUnitTest`  **2 → 14** , zero failures.

## Coverage still missing

The largest remaining gap is that nothing exercises the **optics-dependent** path with golden values — visually-projected RF properties and STFs, which is where the model's scientific claims live. I left that out deliberately: it costs ~20s per case even at 301 wavefront samples, so it wants a `'full'`-mode test rather than `'core'`. The `xOLD/` directory holds the old PLOS2024 and in-vitro validation scripts that once covered this ground; those are the natural basis if you want STF-level regression before release. I'd also suggest extending the mosaic tests across more than the one small mosaic — the current file pins a single artifact, and the NaN defect shows per-mosaic variation matters.
