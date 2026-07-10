# TODO: sceneEye Tutorial and Example Cleanup

This note tracks the remaining work needed to turn the skipped `sceneEye`
scripts into maintained teaching material.  The current tutorial smoke suite
intentionally skips these files because they depend on legacy ISET3d/PBRT,
Docker/cloud helpers, obsolete `sceneEye` API patterns, external scene assets,
or heavyweight rendering workflows.

## Goals

- Keep `tutorials/sceneEye` small, reliable, and tutorial-like.
- Maintain longer ray-tracing workflows in `examples/sceneEye`.
- Keep all `underDevelopment` paths skipped by default.
- Restore only scripts that can run deterministically in routine smoke tests.

## Current State

- Top-level `tutorials/sceneEye` now contains only four candidate tutorials:
  `t_eyeIntro.m`, `t_eyeFOV.m`, `t_eyeFocalDistance.m`, and `t_eyeDoF.m`.
- The other top-level sceneEye workflows were moved to `examples/sceneEye` as
  `s_*.m` files.
- The moved examples still carry `% SkipFile`; they are preserved for repair
  and interactive use, not yet active example smoke-test cases.
- Under-development and cloud sceneEye scripts remain under their existing
  `underDevelopment` paths and are skipped by the tutorial runner.

## Tutorial vs Example Rule

Keep a script in `tutorials/sceneEye` only if it is short, linear, and teaches
one core object/API concept with a compact render or a mocked/precomputed
result.  A tutorial should avoid cloud execution, large PBRT render sweeps,
manual Docker setup, and high-resolution comparisons.

Move a script to `examples/sceneEye` if it is still useful but mainly shows an
applied workflow: multiple eye models, multiple render conditions, crop-window
rendering, depth-of-field comparisons, accommodation sweeps, retina shape
experiments, stereo, cloud rendering, PSF/MTF analysis, or publication-style
figures.

## Next Pass

- Review the four remaining tutorial candidates and decide whether each can be
  made into a short, reliable tutorial.
- Leave `% SkipFile` on moved examples until each one has been repaired and
  tested with
  `isetbioExampleTest('selection', ...)`.
- Replace legacy dependency checks such as `piCamBio`, `piDockerExists`,
  `piDockerConfig`, `mcDockerExists`, and `mcDockerConfig` with current,
  documented setup checks or an explicit skip reason.
- Remove or comment high-resolution sections that repeat the same render with
  only parameter changes.  Keep one representative condition runnable.
- Where a PBRT scene is not in the repository or standard data download path,
  either document the required asset clearly or convert the script to an
  example that remains skipped.

## Candidate Tutorials to Rebuild

These are the best candidates for short, active tutorials after repair:

- `t_eyeIntro.m`: keep as the conceptual entry point, but make it runnable
  without obsolete helper functions.  Prefer a tiny render or a precomputed
  sceneEye object inspection.
- `t_eyeFOV.m`: keep if it can teach field-of-view and retinal geometry without
  a full heavyweight render.
- `t_eyeFocalDistance.m`: keep if reduced to one focal-distance comparison at
  low resolution.
- `t_eyeDoF.m`: keep only if the chess-set dependency is resolved and the render
  count is reduced to one compact depth-of-field demonstration.

## Moved Example Candidates

These files have been moved to `examples/sceneEye` rather than kept in the
active tutorial set.  They remain skipped until each one is repaired and tested
as an example:

- `s_cropWindowExample.m`
- `s_eyeAccommodate.m`
- `s_eyeArizona.m`
- `s_eyeDiffraction.m`
- `s_eyeLeGrand.m`
- `s_eyeMovement3D.m`
- `s_eyeNavarro.m`
- `s_eyeRetinaDistance.m`
- `s_eyeStereo.m`
- `s_retinaShape.m`
- `s_schematicEyeModels.m`
- `s_vergenceAccomm.m`

## Under-Development and Cloud Scripts

The following paths should remain skipped until someone deliberately promotes
individual scripts:

- `tutorials/sceneEye/analysis/underDevelopment/`
- `tutorials/sceneEye/cloud/underDevelopment/`

When promoting one of these, first decide whether it is a compact tutorial or
an applied example.  Most cloud, PSF-over-condition, MTF, eccentricity, and
defocus sweeps should become skipped examples in `examples/sceneEye` unless
they can run quickly and locally.

## Acceptance Checks

- `isetbioTutorialTest` stays at `0` failures.
- All paths containing `underDevelopment` are skipped by the runner.
- Any restored sceneEye tutorial passes
  `isetbioTutorialTest('selection','<tutorialName>')`.
- Any moved sceneEye example passes, or deliberately skips with a clear reason,
  under `isetbioExampleTest('selection','<exampleName>')`.
- A restored tutorial should generally run in under 10 seconds on a normal
  local MATLAB session; longer workflows belong in `examples/sceneEye`.

## Later: Examples Review

Moving useful but non-tutorial scripts out of `tutorials/` will increase the
number of files in `examples/`.  This is less urgent than stabilizing the
tutorial suite, but it should still be reviewed.

- Inventory `examples/` by topic and identify duplicates or near-duplicates.
- Keep examples that show reusable analysis patterns or realistic workflows.
- Remove, merge, or skip examples that are obsolete, overly narrow, broken, or
  mainly publication/data-generation scripts.
- Prefer one clear example per workflow family, with optional parameter
  variants commented rather than executed repeatedly.
- Keep `isetbioExampleTest` passing, or mark heavyweight/resource-dependent
  examples with `% SkipFile` and a specific reason.
