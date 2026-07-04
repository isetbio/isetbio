# Testing ISETBio

This guide is for developers who have changed ISETBio and want to check that
their work has not broken existing behavior. Most users do not need to run
these tests.

ISETBio uses ISETCam's shared test execution and reporting tools. Read the
[ISETCam testing guide](../../isetcam/docs/testing.md) for the common workflow,
the full explanation of `selection`, `start`, `SkipFile`, checkpoints, and
`ieTestReport`.

## Before testing

Start MATLAB with both ISETCam and ISETBio, including their subdirectories, on
the path. ISETCam is a required dependency. A clean MATLAB session is
recommended for tutorial and example runs because these runners reset ISET
state and close figures between scripts.

ISETBio provides three repository-wide runners:

| Runner | Purpose |
| --- | --- |
| `isetbioUnitTest` | Run automated function-based unit tests |
| `isetbioTutorialTest` | Run tutorial scripts as smoke tests |
| `isetbioExampleTest` | Run example scripts as smoke tests |

## A practical testing workflow

During development:

1. Run the nearest subsystem unit-test runner.
2. Run directly related tutorials or examples with `'selection'`.
3. Run the core repository unit tests before sharing a substantial change.
4. Use the full unit suite and complete script suites for broad changes.

## Unit tests

Run the core suite, which excludes test files whose names contain `FullOnly`:

```matlab
results = isetbioUnitTest;
```

Run the full suite:

```matlab
results = isetbioUnitTest('full');
```

Subsystem `_tests_` directories provide focused runners such as:

```matlab
results = conesUnitTest;
results = eyemovementUnitTest;
results = opticalimageUnitTest;
```

Focused ISETBio runners use the same default `core` mode and optional `'full'`
mode. Unit runners return MATLAB `TestResult` arrays suitable for
`ieTestReport`.

## Tutorial and example tests

Run every tutorial or example:

```matlab
tutorialRun = isetbioTutorialTest;
exampleRun = isetbioExampleTest;
```

Run one script:

```matlab
run = isetbioTutorialTest('selection','t_cMosaicBasic');
run = isetbioExampleTest('selection','s_matlabConv2');
```

Run from one script through the remainder of the sorted plan:

```matlab
run = isetbioTutorialTest('start','t_cMosaicBasic');
```

Use the canonical marker for a script that must be discovered but skipped:

```matlab
% SkipFile
```

Scripts whose purpose is to generate or refresh repository data files should
be named `data_*.m`, not `t_*.m` or `s_*.m`.

## Reports and interrupted runs

The tutorial and example runners return the same run structure used by
ISETCam. List failures or skips with:

```matlab
ieTestReport(run,'List',{'failed','skipped'});
```

Runs save `checkpoint.mat`, `progress.log`, and `planned-files.txt` in a
timestamped directory below the ISETBio `local/` directory. If MATLAB exits
before returning a variable, pass the checkpoint or run directory directly to
`ieTestReport`:

```matlab
ieTestReport('/path/to/run/directory','List','all');
```

For additional reporting and failure-investigation examples, see the
[ISETCam testing guide](../../isetcam/docs/testing.md).
