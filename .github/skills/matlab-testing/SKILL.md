---
name: matlab-testing
description: Develop, run, diagnose, or report ISETBio MATLAB tests. Use when changing production MATLAB code, adding function-based tests or local _tests_ runners, running isetbioUnitTest, isetbioTutorialTest, or isetbioExampleTest, investigating SkipFile behavior, checkpoints, or ieTestReport output.
---

# ISETBio MATLAB Testing

## Instructions & Guidelines

### Choose the Smallest Relevant Check

ISETBio uses ISETCam's shared test engine and reporting tools. Begin with the
nearest subsystem runner, then run directly related tutorial or example smoke
tests, and run the repository core suite before sharing a substantial change.

```matlab
results = isetbioUnitTest;          % core unit suite
results = isetbioUnitTest('full');  % includes FullOnly tests

results = conesUnitTest;            % focused subsystem examples
results = eyemovementUnitTest;
results = opticalimageUnitTest;
```

Run all teaching scripts only for broad changes:

```matlab
tutorialRun = isetbioTutorialTest;
exampleRun = isetbioExampleTest;
```

Run one script with `selection`, or resume from a named script with `start`:

```matlab
run = isetbioTutorialTest('selection','t_cMosaicBasic');
run = isetbioExampleTest('selection','s_matlabConv2');
run = isetbioTutorialTest('start','t_cMosaicBasic');
```

### Test Authoring

- Place tests with the subsystem or behavior they protect, in a colocated
  `_tests_` directory.
- Name function-based test files `test_<subject>.m` and begin them with
  `tests = functiontests(localfunctions)`.
- Test accessors, computations, dimensions/shapes, invariants, validation
  errors, and stable golden values using explicit named tolerances.
- Keep core tests deterministic and non-interactive. Control randomness and
  classify GUI, smoke, slow, or resource-heavy cases outside the core suite.
- Provide each `_tests_` directory a local `<area>UnitTest.m` runner based on
  `TestSuite.fromFolder`, `TestRunner.withTextOutput`, and `ieTestReport`.
  Default to `core` and accept `full` for all cases.
- Runners must close figures created by the test run while preserving figures
  that existed before it started.
- Do not duplicate an ISETCam test when ISETCam already owns the behavior.
  Migrate legacy `isetvalidate` checks to the subsystem they protect.

### Script Runs, Skips, and Reports

Tutorial and example runners discover `t_*` and `s_*` files. To deliberately
discover but skip a resource-dependent or interactive script, add this exact
comment anywhere in the file:

```matlab
% SkipFile
```

Use it sparingly and remove it when routine execution becomes practical. The
legacy `% UTTBSkip` marker remains compatible but new work uses `% SkipFile`.

Runs save `checkpoint.mat`, `progress.log`, and `planned-files.txt` below a
timestamped directory in `local/`. Report a completed or interrupted run with:

```matlab
ieTestReport(run,'List',{'failed','skipped'});
ieTestReport('/path/to/run/directory','List','all');
```

When a sandboxed MATLAB batch run silently fails or exits with status 1,
retry with permission to write MATLAB preferences or caches outside the repo.

## Known coverage gap

Nothing pins the optics-dependent mRGC path with golden values: visually
projected RF properties and spatial transfer functions. It was left out on
cost, roughly 20 seconds per case, so it belongs behind an opt-in `FullOnly`
test rather than the core suite. See `docs/todo-science.md` for the details and
for the known NaN surround behavior that the pre-baked circuit tests pin.
