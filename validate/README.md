# ISETBio validation

The public ISETBio validation runners live in this directory. They use
ISETCam's shared execution and reporting infrastructure, so ISETCam must be on
the MATLAB path.

## Unit tests

Run the core unit-test suite:

```matlab
results = isetbioUnitTest;
```

Run the full suite, including files whose names contain `FullOnly`:

```matlab
results = isetbioUnitTest('full');
```

Subsystem runners in colocated `_tests_` directories support the same `core`
and `full` modes.

## Tutorials and examples

Run every tutorial or example:

```matlab
tutorialRun = isetbioTutorialTest;
exampleRun = isetbioExampleTest;
```

Run one file with `selection`. The value may be a script stem, file name,
path relative to the suite directory, or full path:

```matlab
run = isetbioTutorialTest('selection','t_cMosaicBasic');
run = isetbioExampleTest('selection','s_someExample.m');
```

Run one file and every file after it in the deterministic, path-sorted plan:

```matlab
run = isetbioTutorialTest('start','t_cMosaicBasic');
```

Add `% SkipFile` to a tutorial or example that should be discovered but not
executed. Skips should be rare and documented in the source file.

## Reports and checkpoints

Each tutorial or example run creates a timestamped directory below `local/`
with `checkpoint.mat`, `progress.log`, and `planned-files.txt`. The checkpoint
is updated before and after each script and remains useful if MATLAB exits.

The runner prints a summary automatically. Reports can also be regenerated or
expanded with ISETCam's `ieTestReport`:

```matlab
ieTestReport(run)
ieTestReport(run,'List',{'failed','skipped'})
ieTestReport('/path/to/run/directory','List','all')
```

For implementation details, see ISETCam's
`docs/tutorial-example-test-architecture.md`.
