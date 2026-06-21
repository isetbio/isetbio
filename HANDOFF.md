# Handoff: tutorial runner and MATLAB Desktop crashes

Date: 2026-06-21  
Repository: `/Users/wandell/Documents/MATLAB/isetbio`  
Dependency: `/Users/wandell/Documents/MATLAB/isetcam`

## Start here

Read `.github/copilot-instructions.md` before making changes. The shared
tutorial/example runner architecture is documented in:

`../isetcam/docs/tutorial-example-test-architecture.md`

The immediate unresolved problem is a repeatable native MATLAB Desktop crash
associated with `cMosaic` plotting. Native crashes terminate MATLAB and bypass
MATLAB `try/catch`. The same tutorials pass under `matlab -batch`, so batch
success is not sufficient evidence that the Desktop problem is fixed.

Do not replace `cMosaic.plot` with `scatter` as the final solution. The user
explicitly wants to preserve `cMosaic.plot`; it is a long-standing workhorse
that has worked repeatedly in the past.

## Current working tree

At handoff, these files have uncommitted changes:

- `cones/cmosaic/@cMosaic/plot.m`
- `tutorials/adaptiveOptics/t_AODisplay.m`

The current `t_AODisplay` once again ends with the desired call:

```matlab
theExcitations = theMosaic.compute(theOI);
theMosaic.plot('excitations',theExcitations);
```

Commit `dea8b0b9e` temporarily replaced this with `scatter`. The working-tree
change restores `cMosaic.plot`. Do not accidentally restore the scatter
workaround.

`t_AODisplay.m` has a CRLF/LF line-ending diff in addition to the substantive
change. Review diffs with `git diff --ignore-space-at-eol`.

## Crash chronology and evidence

### First full-suite crash

Run directory:

`local/2026-06-21_135331_isetbioTutorialsTest`

The progress log contained:

1. `t_AOBleachingKinetics` passed.
2. `t_AODisplay` was marked `ScriptStarted`.
3. MATLAB crashed before any pass/fail event was written.

The runner records `ScriptStarted` before its reset and execution, so that
checkpoint alone did not initially prove which operation crashed. The user
then ran `t_AODisplay` manually and found the native exit at its final line:

```matlab
theMosaic.plot('excitations',theExcitations);
```

### Numerical diagnostics for `t_AODisplay`

The tutorial was run in a disposable MATLAB R2025b batch process with the
rendering call changed to `'data only',true`. Results:

- cones: `12085`
- excitation size: `[1 1 12085]`
- all excitation values finite
- excitation range: approximately `[0.0104546, 286791]`

Thus there was no obvious size mismatch, NaN, Inf, or malformed response.
The complete original tutorial also ran successfully under `matlab -batch`.
Again, Desktop behavior differs from batch behavior.

### Failed hypotheses and experiments

1. **Colormap size**
   
   `cMosaic.visualize` creates `gray(numberOfCones)`. This was temporarily
   changed to `gray(256)` and tested. MATLAB Desktop still crashed. The change
   was reverted. Do not repeat this experiment as the primary fix.

2. **Scene/OI web-app interaction**
   
   `sceneWindow` and `oiWindow` were temporarily replaced with
   `sceneShowImage`/`oiShowImage` and ordinary figures. The user reported an
   even earlier crash, before either preview appeared. Those changes were
   reverted. The original `sceneWindow` and `oiWindow` calls are present now.

3. **Scatter workaround**
   
   Replacing the final detailed plot with a lightweight scatter map allowed
   `t_AODisplay` to complete after residual MATLAB processes were killed. The
   user does not accept abandoning the workhorse plot method, correctly. The
   working tree restores `cMosaic.plot`.

4. **Residual MATLAB processes**
   
   Killing residual MATLAB processes allowed MATLAB Desktop to start and the
   scatter-based tutorial to complete, but did not establish that residual
   processes caused the original plot crash. Cleanup commands used were:

   ```sh
   pkill -TERM -f '/Applications/MATLAB'
   sleep 3
   pkill -KILL -f '/Applications/MATLAB'
   ```

   Warn the user that these commands also stop MATLAB used by VS Code and can
   destroy unsaved MATLAB work.

### Second full-suite crash

Run directory:

`local/2026-06-21_143856_isetbioTutorialsTest`

The run passed the scatter-based `t_AODisplay`, then passed through
`t_cMosaicArbitraryOpticalImagePositioning`. It crashed after recording:

```text
ScriptStarted | [8/72] .../tutorials/cmosaic/t_cMosaicBasic.m
```

No later event was written. `t_cMosaicBasic` contains many `cMosaic.plot`
calls. Its first is:

```matlab
cm.plot('mosaic');
```

Therefore the suite checkpoint implicates `t_cMosaicBasic`, but it does not
prove which line inside that script caused the native exit.

### Latest manual crash

After restoring the historical invocation path described below, the focused
tests and batch tutorial passed, but the user reported another MATLAB Desktop
crash. This appears to have been a manual Desktop run, because no newer full
suite checkpoint was present at handoff. Ask exactly which tutorial/call was
running before drawing a new conclusion.

## Recent `cMosaic.plot` refactor

Commit `97b52c37f` (2026-06-12, "Refactor cMosaic.plot and update tests/examples")
substantially refactored `cones/cmosaic/@cMosaic/plot.m`.

The low-level polygon renderer in `cMosaic.visualize` and its
`renderPatchArray` helper have been essentially unchanged since 2023. This
makes the recent `plot.m` refactor a plausible regression boundary, although
the latest Desktop crash means the current hypothesis is not confirmed.

The current uncommitted `plot.m` changes restore the historical invocation
style for the two workhorse cases:

- `cm.plot('mosaic')` directly calls `cmosaic.visualize`.
- `cm.plot('excitations',...)` obtains `cmosaic.visualize('params')`, fills
  the parameter struct, and calls `cmosaic.visualize(params)`.

The newer trial/time-point selection and supplied figure/axes support are
retained for excitation plots. A diagnostic comment notes that native exits
bypass `try/catch` and should be treated as a graphics/resource regression,
not automatically as invalid excitation data.

This historical-invocation change passed tests and batch execution but did
not, according to the user's latest message, solve the Desktop crash.

## Validation already completed

Focused tests passed after restoring the historical invocation path:

```matlab
results = runtests( ...
    '/Users/wandell/Documents/MATLAB/isetbio/cones/cmosaic/_tests_/test_cMosaic.m');
assertSuccess(results);
```

The real focused tutorial runner also passed in MATLAB R2025b batch mode:

```matlab
run = isetbioTutorialsTest('t_cMosaicBasic');
```

Checkpoint:

`local/2026-06-21_144406_isetbioTutorialsTest`

Result: one planned, one passed, run completed.

## Crash reports

No new macOS diagnostic report was found for the afternoon crashes. The most
recent file found was:

`~/Library/Logs/DiagnosticReports/MATLAB-2026-06-21-104201.ips`

It predates these particular runs. It reports R2025b, macOS 26.5,
`EXC_BAD_ACCESS`/`SIGSEGV` on the Qt GUI thread, with Chromium/font-rendering
symbols. This is evidence of general Desktop graphics instability on this
machine, but it must not be presented as the stack trace for the later
`cMosaic.plot` crashes.

Older MATLAB crash dumps showed native graphics/hardcopy frames, but they also
predate this investigation.

## Suggested next investigation

Proceed conservatively; repeated Desktop crashes are costly for the user.

1. Start with a fresh MATLAB Desktop after confirming residual MATLAB
   processes are gone.
2. Ask the user which exact script/call produced the latest crash. Do not infer
   it from the completed batch checkpoint.
3. Reproduce with a very small standalone `cMosaic` and `cm.plot('mosaic')`,
   then increase mosaic size. This separates basic mosaic rendering from
   excitation rendering and from tutorial state.
4. Compare R2025b Desktop with another installed MATLAB release if practical.
   A release-specific graphics regression is plausible because batch passes
   and the plotting code historically worked.
5. Inspect MathWorks release notes/bug reports for R2025b on macOS 26.5,
   especially Qt/Chromium/patch graphics native exits. Use official MathWorks
   sources when browsing.
6. If the problem is isolated to `renderPatchArray`, preserve aperture
   rendering while reducing graphics payload correctly. A promising internal
   optimization is to use one `FaceVertexCData` value per face rather than
   repeating the same value for every vertex, and to preallocate `Faces`
   instead of growing it with `cat`. This has not yet been tried and should be
   tested for visual equivalence; do not silently replace apertures with
   scatter points.
7. Consider a focused visual regression test that checks patch face count,
   color mapping, and supplied axes. Batch tests cannot detect the Desktop
   native exit, but they can protect rendering semantics during optimization.

## Tutorial/example runner context

The shared engine is:

`../isetcam/utility/ieRunTutorialExampleTests.m`

ISETBio wrappers are:

- `tutorials/isetbioTutorialsTest.m`
- `examples/isetbioExamplesTest.m`

ISETCam has corresponding wrappers. Runs save an atomic `checkpoint.mat`, a
`progress.log`, and `planned-files.txt`. `ieTestReport` in ISETCam reads either
completed returns or crash-surviving checkpoints and can list passed, failed,
or skipped files.

An incomplete checkpoint remains in state `Running` because a dead MATLAB
process cannot update it to `Interrupted`. `currentFile` means the file that
had been marked started; it does not always prove the script body began,
because the runner currently records `ScriptStarted` before its reset call.

Also note that the runner's pre/post tutorial reset calls are outside the
inner script `try/catch`. Ordinary reset errors should eventually be included
in catchable per-file handling, but no `try/catch` can contain a native MATLAB
process crash.

## User priorities

- Preserve `cMosaic.plot`; do not give up on this established API.
- Prevent tutorial runs from accumulating figures and session data.
- Tutorials should not require a parallel pool.
- Long tutorials should be shortened or skipped when justified.
- Crash checkpoints must remain readable because native crashes return no
  MATLAB value.

