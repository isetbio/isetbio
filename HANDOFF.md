# Handoff: Tutorial Smoke Test Stabilization

Date: 2026-07-07
Branch: `main-unit-merge`
Current anchor commit: `f791b2202 Make tutorials faster and parfor optional`

## Context

This work was driven by repeated crashes and very long runtimes while running
`isetbioTutorialTest`. The immediate user goal is to make the tutorial smoke
suite suitable for routine automated execution: quick, deterministic, and not
dependent on heavyweight publication/data-generation workflows.

Repository instructions live in `.github/copilot-instructions.md`. Read that
file before continuing.

## What Was Found

The tutorial failures were not simply user misremembering. There were several
different issues:

- Some `t_*` files are really data-generation, publication, or prebaked-resource
  workflows rather than short tutorials.
- Some tutorial computations had grown too large for smoke tests.
- `t_cMosaicGenerate.m` had a severe size regression: a from-scratch cMosaic
  path was using a large default source lattice (`58 deg`) and spent tens of
  minutes processing roughly 100 million candidate RFs.
- `t_cMosaicComputeWithOIEnsemble.m` was doing a much larger OI ensemble
  calculation than a tutorial smoke test should require.
- Several cMosaic operations were starting `parfor`/parpool paths even for
  small tutorial-sized work.
- `visualStimulusGenerator.presentationDisplay` called a missing helper:
  `generateConventionalxyYDisplayDefaultParams`.

## Main Changes Already Made

These changes are included at the anchor commit above.

### Runtime reductions

- `tutorials/cmosaic/t_cMosaicGenerate.m`
  - Reduced tutorial-scale mosaic sizes.
  - Added explicit small source lattice for the from-scratch case.
  - Disabled `parfor` for tutorial examples.

- `tutorials/cmosaic/t_cMosaicComputeWithOIEnsemble.m`
  - Reduced stimulus size, image pixels, mosaic dimensions, OI sampling grid,
    and wavefront sample count.
  - Skips expensive visualizations while running under
    `ieRunTutorialExampleTests`.

- `tutorials/cmosaic/t_cMosaicCustomConeData.m`
  - Disabled `parfor` for the small custom-data tutorial mosaic.

- `tutorials/cmosaicrect/t_cmosaicRectBigArray.m`
  - Reduced scene FOV from `3` to `1`.
  - Skips interactive `coneMosaicRect.window` calls under the tutorial runner.

- `tutorials/eyemovement/t_fixationalEyeMovementsCharacterize.m`
  - Reduced eye movement duration/trials from `3 sec`, `1 ms`, `512 trials` to
    `1 sec`, `2 ms`, `32 trials`.
  - Made log-axis plot limits data-driven so the compact run still plots
    safely.

### Core cMosaic behavior

- `cones/cmosaic/@cMosaic/computeConeApertures.m`
  - Respects `obj.useParfor`; uses a plain `for` loop when disabled.

- `cones/cmosaic/@cMosaic/updateStateGivenKeptConeIndices.m`
  - Respects `obj.useParfor`; uses a plain `for` loop when disabled.

### Missing helper fix

- `isettools/ganglioncells/+visualStimulusGenerator/presentationDisplay.m`
  - Removed dependency on missing `generateConventionalxyYDisplayDefaultParams`
    and `generateConventionalxyYDisplay`.
  - Uses ISETCam display constructors/setters instead.

### Smoke-test skips

The following were marked with `% SkipFile` because they are heavyweight,
resource-dependent, prebaked mRGC, publication, or data-generation workflows:

- `tutorials/cmosaic/t_customHumanConeMosaicsBasedOnAOmeasurements.m`
- `tutorials/mrgc/advanced/t_mRGCMosaicEnableSurroundSconeInputs.m`
- `tutorials/mrgc/t_mRGCMosaicBasicComputation.m`
- `tutorials/mrgc/t_mRGCMosaicInspect.m`
- `tutorials/mrgc/t_mRGCMosaicSynthesizeAtStage2.m`
- `tutorials/mrgc/t_mRGCMosaicSynthesizeAtStage3.m`
- `tutorials/mrgc/t_mRGCMosaicVisualizeWithOptics.m`

Several other mRGC tutorials were already skipped before this work.

## Validation Already Run

MATLAB executable used:

```sh
/Applications/MATLAB_R2025b.app/bin/matlab
```

Batch MATLAB sessions need paths initialized explicitly:

```matlab
addpath(genpath(fullfile(getenv('HOME'),'Documents','MATLAB','isetcam')));
addpath(genpath(fullfile(getenv('HOME'),'Documents','MATLAB','isetbio')));
```

Focused checks that passed:

- `isetbioTutorialTest('selection','t_cMosaicGenerate')`
- `isetbioTutorialTest('selection','t_cMosaicComputeWithOIEnsemble')`
- `isetbioTutorialTest('selection','t_cmosaicRectBigArray')`
- `isetbioTutorialTest('selection','t_fixationalEyeMovementsCharacterize')`

Focused skip checks reported `SKIPPED` as intended:

- `t_customHumanConeMosaicsBasedOnAOmeasurements`
- `t_mRGCMosaicEnableSurroundSconeInputs`
- `t_mRGCMosaicBasicComputation`
- `t_mRGCMosaicInspect`
- `t_mRGCMosaicSynthesizeAtStage2`
- `t_mRGCMosaicSynthesizeAtStage3`
- `t_mRGCMosaicVisualizeWithOptics`

Useful observed runtimes:

- `t_cmosaicRectBigArray`: about 4 seconds after reduction.
- `t_fixationalEyeMovementsCharacterize`: about 16-17 seconds after reduction
  (previously about 85 seconds in the user's full-run log).
- `t_cMosaicGenerate`: no longer hits the huge 58 deg source lattice path.

Repository check:

```sh
git diff --check
```

This was clean before writing this handoff file.

## Most Recent Full-Run Evidence From User

The user's most recent crash report before these final skips/reductions:

```text
--- isetbioTutorialTest Run Summary ---
Run State:        Running
Total Planned:    83
Total Completed:  47
Total Passed:     37
Total Failed:     1
Total Skipped:    9
Total Unfinished: 36
Last active file: mrgc/advanced/t_mRGCMosaicEnableSurroundSconeInputs.m
Last checkpoint:  2026-07-07 08:04:38
```

The failure before the crash was:

- `tutorials/cmosaic/t_customHumanConeMosaicsBasedOnAOmeasurements.m`
  - Failed in local data parsing (`cell2mat` in `parseSingleMosaicData`).
  - This file is now skipped for smoke tests.

The crash-active file was:

- `tutorials/mrgc/advanced/t_mRGCMosaicEnableSurroundSconeInputs.m`
  - This file is now skipped for smoke tests.

## Recommended Next Steps On The Other Computer

1. Pull/check out the same branch:

   ```sh
   git checkout main-unit-merge
   git pull
   ```

2. Confirm the anchor commit is present:

   ```sh
   git log -1 --oneline
   ```

   Expected top commit at this handoff: `f791b2202 Make tutorials faster and parfor optional`.

3. Run a focused post-crash continuation cluster first:

   ```matlab
   addpath(genpath(fullfile(getenv('HOME'),'Documents','MATLAB','isetcam')));
   addpath(genpath(fullfile(getenv('HOME'),'Documents','MATLAB','isetbio')));
   isetbioTutorialTest('start','t_customHumanConeMosaicsBasedOnAOmeasurements')
   ```

   Watch that the AO/Sabesan and mRGC heavy files report skipped and the run
   advances past the old crash point.

4. If that looks good, run the full tutorial smoke suite:

   ```matlab
   isetbioTutorialTest
   ```

5. If MATLAB crashes again, inspect the newest progress log under:

   ```text
   local/*_isetbioTutorialTest/progress.log
   ```

   The `Last active file` in the checkpoint/run summary is the next file to
   triage. Prefer reducing tutorial computation sizes first. If a file is really
   a publication/data-generation/resource workflow, mark it with `% SkipFile`
   and leave a nearby comment explaining why.

## Current Handoff File Status

`HANDOFF.md` itself was written after the clean working tree check. If this
file is uncommitted on the source machine, commit or copy it before switching
computers.

## 2026-07-07 Follow-Up On This Machine

The post-crash continuation run completed successfully:

```matlab
isetbioTutorialTest('start','t_customHumanConeMosaicsBasedOnAOmeasurements')
```

Result:

- Total planned: 63
- Passed: 39
- Failed: 0
- Skipped: 24
- Unfinished: 0

Run log:

```text
local/2026-07-07_094305_isetbioTutorialTest/progress.log
```

Additional smoke-mode runtime reductions were then made for passing tutorials
that were still relatively slow:

- `tutorials/recipes/t_dynamicStimulusFixationalEM.m`
  - Uses a shorter duration, lower frame rate, smaller matched stimulus/mosaic
    FOV, and one trial for both direct Desktop use and tutorial smoke runs.
    This avoids the previous warning flood and keeps the script student-facing
    without top-level test-runner branching.
- `tutorials/oisequences/t_oisCreate.m`
  - Uses shorter temporal weights and skips movie generation under the tutorial
    runner.
- `tutorials/oisequences/t_oisHarmonic.m`
  - Uses smaller image dimensions, shorter temporal weights, and skips
    movie/montage visualization under the tutorial runner.
- `tutorials/optics/t_wvfWatsonJOV.m`
  - Uses fewer scale/Zernike demonstrations and suppresses GUI-heavy plots
    under the tutorial runner.
- `tutorials/outersegment/advancedTutorials_os/t_osTimeStep.m`
  - Runs fewer conditions, uses a shorter simulation time axis, and skips the
    large plotting routine under the tutorial runner.
- `tutorials/cones/t_conesPerifovea.m`
  - Uses a smaller rings/rays stimulus and FOV under the tutorial runner.

Focused checks passed for all six edited tutorials.

The full tutorial smoke suite then completed successfully:

```matlab
isetbioTutorialTest
```

Result:

- Total planned: 83
- Passed: 58
- Failed: 0
- Skipped: 25
- Unfinished: 0
- Started: 2026-07-07 09:55:12
- Finished: 2026-07-07 10:02:52

Run log:

```text
local/2026-07-07_095512_isetbioTutorialTest/progress.log
```

The remaining longer passing tutorials in that full run were mostly around
20-30 seconds:

- `t_cMosaicHarmonic`
- `t_cMosaicBasic`
- `t_RGBImageToExcitations`
- `t_dynamicStimulusToPhotocurrent`
- `t_linearFilters`

These are good next candidates if the goal is to push the full smoke suite
below its current roughly 8-minute runtime.
