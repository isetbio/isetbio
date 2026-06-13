# ISETBio AI Instructions

Use this file as the shared startup guidance for Copilot, Claude, Codex,
Gemini, and other AI coding assistants working in this repository.

## Repository Context

- MATLAB is the primary runtime.
- The main repository is `isetbio`. ISETCam (`../isetcam`) is a required
  dependency and is always expected to be on the MATLAB path when ISETBio is
  used or tested.
- ISETBio code and tests may directly use ISETCam utilities, including
  `ieTestReport`. Do not duplicate utilities already supplied by ISETCam.
- related local repositories may include `isetvalidate`, and `tools/UnitTestToolbox`.
- Many independently maintained repositories depend on ISETBio. Before removing
  or changing public APIs, paths, data locations, setup behavior, or integration
  hooks, search for likely external usage and prefer staged deprecation when an
  immediate change could disrupt collaborators. Allow time for dependent
  repositories to migrate unless coordinated cleanup is explicitly requested.
- For VS Code MATLAB setup, see `.vscode/matlab-setup.md`.
- For MATLAB Command Window path setup, use `.github/matlab-paths.md`.

## Tutorials and Examples

ISETBio keeps `tutorials/` and `examples/` as separate teaching surfaces for
different goals and audiences.

- **Tutorials (`tutorials/`)**

  - Audience: learners (including new students) who can program and are
    learning image systems engineering and ISETCam object fundamentals.
  - Purpose: short, heavily commented introductions to key objects and APIs.
  - Expected content:
    - object creation and setup
    - `*Get`/`*Set` usage for key properties
    - basic visualization (`*Window`, `*Plot`)
    - one simple quantitative computation/checkpoint
  - Expected behavior: runs relatively quickly and is easy to read linearly.
- **Examples (`examples/`)**

  - Audience: users looking for realistic analysis patterns to adapt.
  - Purpose: applied workflows and more advanced computations using ISETCam.
  - Expected content:
    - end-to-end numerical analyses or visualization workflows
    - realistic parameter choices and tradeoff exploration
    - code that users may copy/adapt as a starting point for their own work
  - Expected behavior: can be longer and more detailed than tutorials.

When adding or editing files, preserve this distinction. If content is mainly
onboarding and API orientation, place it in `tutorials/`. If content is mainly
applied workflow, analysis, or deeper exploration, place it in `examples/`.

You can convert these tutorials and examples into HTML documentation by running
the `s_publishTutorials` and `s_publishExamples` utilities (provided by ISETCam)
from the MATLAB command window.

For student contributors, prioritize clarity, reproducibility, and instructional
value: use clear comments, stable outputs, and explicit links to related wiki
pages, tests, and nearby tutorials/examples.

## ISETCam Pipeline

Prefer existing object-specific functions before writing new utilities.

1. Scene: `scene*` functions, accessed with `sceneGet` and `sceneSet`.
2. Optical image: `oi*` functions, accessed with `oiGet` and `oiSet`.
3. Sensor: `sensor*` functions, accessed with `sensorGet` and `sensorSet`.
4. Image processing: `ip*` functions, accessed with `ipGet` and `ipSet`.
5. Display: `display*` functions, accessed with `displayGet` and `displaySet`.

Common constructors and compute functions include `sceneCreate`,
`oiCreate`, `oiCompute`, `sensorCreate`, `sensorCompute`, `ipCreate`,
`ipCompute`, and `displayCreate`.

For object diagnostics, prefer existing plotting functions such as
`scenePlot`, `oiPlot`, `sensorPlot`, `ipPlot`, and `displayPlot` over ad hoc
plotting.

## Search Guidance

- Use `rg` for text search and `fd` for filename/path search when using a
  terminal.
- Before adding behavior, search for nearby examples with the relevant object
  prefix.
- For color transforms and color science utilities, search `color/` before
  implementing new code.
- For new scene patterns or chart behavior, check existing examples in
  `scene/` and especially related pattern/chart code.

## Coding Style

- Keep edits minimal and consistent with existing MATLAB style.
- Reuse established constructors, getters, setters, plotting helpers, and
  object naming conventions.
- Prefer vectorized MATLAB where it improves clarity or performance.
- Update function header comments when behavior changes, especially `Syntax`,
  `Inputs`, `Returns`, and `See also`.
- Do not add dependencies unless they are necessary and consistent with the
  repository.

## Validation

- Validate modified files with MATLAB diagnostics or focused test commands when
  practical.
- Place tests for major objects and computational areas in colocated `_tests_`
  directories. Use ISETCam's `_tests_` directories as the reference
  implementation when an ISETBio convention is not yet established.
- Write function-based MATLAB tests in files named `test_<subject>.m`, starting
  each file with `tests = functiontests(localfunctions)`.
- Prefer focused, descriptively named test functions that cover accessors,
  computations, dimensions and shapes, invariants, important validation
  errors, and stable golden-value fingerprints with explicit named tolerances.
- Keep core tests deterministic and non-interactive. Control random-number
  generation when randomness is required, and classify GUI, smoke, slow, or
  resource-heavy tests outside the core suite.
- Give each `_tests_` directory a local `<area>UnitTest.m` runner built with
  `TestSuite.fromFolder`, `TestRunner.withTextOutput`, and `ieTestReport`.
  Local runners should run the `core` suite by default and accept `full` to
  include all tests.
- Local and repository-wide runners must close figures created during testing
  while preserving figures that were open before the test run.
- Run the full ISETBio unit-test suite with `isetbioUnitTests` and render or
  summarize its output with ISETCam's `ieTestReport`. `isetbioUnitTests` is
  the ISETBio master runner; `ieUnitTests` remains the ISETCam master runner.
- When converting legacy `isetvalidate` scripts into built-in unit tests,
  place each test with the ISETBio subsystem or behavior it protects rather
  than copying the legacy validation directory layout. Do not duplicate a
  test already maintained by ISETCam merely because the validation script
  historically lived under an ISETBio validation directory.
- Treat `isetvalidate` as the broader system/regression validation suite when
  relevant to a change.
- MATLAB is available through the VS Code MATLAB extension.
- A local MATLAB executable is available at
  `/Applications/MATLAB_R2025b.app/bin/matlab` and can be used with `-batch`
  for non-interactive checks.
- If launching MATLAB from a sandboxed shell fails silently or exits with
  status 1, retry unsandboxed or escalated because MATLAB may need to write
  preferences or cache files outside the repository.

## When Uncertain

Choose the simplest implementation that matches existing `scene*`, `oi*`,
`sensor*`, `ip*`, and `display*` patterns. Ask the user only when the choice
would materially affect behavior, API shape, or test expectations.
