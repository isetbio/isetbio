# Watson (2014) retinal distance conversions

ISETBio uses Watson's formulae to convert between visual degrees and retinal
millimeters:

- `RGCmodels.Watson.convert.rhoDegsToMMs` — Appendix Equation A5
- `RGCmodels.Watson.convert.rhoMMsToDegs` — Appendix Equation A6

Near the fovea the A5 linear coefficient is the familiar 0.268 mm/deg.

## A5 and A6 are not exact inverses

The two equations are **separately fitted polynomials**, not an analytic
inverse pair, and their linear terms disagree:

| Source | Linear term |
|---|---|
| A5, inverted | 1 / 0.268 = 3.731 |
| A6, as published | 3.556 |

A degrees → mm → degrees round trip therefore does not return its input:

| Eccentricity | Round-trip error |
|---|---|
| 0.5° | −4.4% |
| 1° | −4.2% |
| 2° | −3.7% |
| > 20° | converges |

The discrepancy is largest exactly where it matters most. Most pre-baked
mosaics sit within a few degrees of the fovea, and 14 files in the repository
call both conversions.

## Why this is not "fixed"

Both implementations match the publication. Changing either one would make
ISETBio disagree with Watson (2014), which is worse than the round-trip error:
results would no longer be comparable with the published model or with other
implementations of it. The behavior is pinned by a test rather than corrected.

## What to do about it

- **Do not rely on a degrees → mm → degrees round trip** to recover an input
  value near the fovea. Keep the original value instead of recovering it.
- Convert once, in one direction, wherever possible. Chained conversions
  compound the error.
- When comparing a quantity in degrees against one derived through millimeters,
  expect a few percent of disagreement at low eccentricity, and check whether
  that is within the tolerance your comparison needs.

## Where this is verified

`ganglioncells/_tests_/test_RGCmodelsGoldenValues.m`:

- `testWatsonRhoDegsToMMs` — A5 golden values, including 0.2683 mm at 1°
- `testWatsonRhoMMsToDegs` — A6 golden values
- `testWatsonConversionsAreNotExactInverses` — pins the round-trip error at
  0.5°, 1°, and 2°, so a change in either equation fails distinctly

The same file pins the Croner & Kaplan relationships used alongside these
conversions: the 6.7 surround-to-center Rc ratio, surround radius versus
eccentricity, and the integrated sensitivity ratio. Those golden values were
checked against the literature rather than recorded from output.

## Reference

Watson, A. B. (2014). A formula for human retinal ganglion cell receptive field
density as a function of visual field location. *Journal of Vision*, 14(7):15.
<https://doi.org/10.1167/14.7.15>

See also `tutorials/data/t_data_rgcWatson.m`, which plots the density formula.
