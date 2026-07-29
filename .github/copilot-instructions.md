# ISETBio AI Instructions

Use this file as the shared startup guidance for AI coding assistants working
in this repository. Load a skill in `.github/skills/` when its description
matches the requested work; skills contain task-specific procedures.

## Repository Context

- MATLAB is the primary runtime.
- ISETCam (`../isetcam`) is a required dependency and must be on the MATLAB
  path when ISETBio is used or tested. Reuse ISETCam utilities, including
  `ieTestReport`; do not duplicate them in ISETBio.
- Related local repositories can include `isetvalidate` and
  `tools/UnitTestToolbox`.
- Many independently maintained repositories depend on ISETBio. Before
  changing a public API, path, data location, setup behavior, or integration
  hook, search for likely external callers and prefer a staged deprecation
  unless coordinated cleanup is explicitly requested.

## Repository Layout

- Computational subsystems include `cones/`, `opticalimage/`, `eyemovement/`,
  `ganglioncells/`, `mouse/`, and `wrappers/`.
- Teaching material is in `tutorials/` and `examples/`; automated validation
  runners are in `validate/`.
- Repository-owned data is in `data/`; developer documentation is in `docs/`.

## Universal Engineering Rules

- Keep edits minimal and match nearby MATLAB conventions.
- Reuse established constructors, getters, setters, plotting helpers, and
  object naming conventions before adding utilities.
- Prefer vectorized MATLAB only when it improves clarity or performance.
- Update a function header's `Syntax`, `Inputs`, `Returns`, and `See also`
  comments when behavior changes.
- Do not add dependencies unless necessary and consistent with the repository.
- Use `rg` for content search and `fd` for filename search when available.
- Validate changed code with the narrowest relevant MATLAB diagnostic or test.
  Use `isetbioUnitTest` for the repository-wide unit suite.

## Skills

Activate the matching skill before following detailed workflow guidance:

- `matlab-environment` — MATLAB paths and VS Code setup.
- `matlab-testing` — unit, tutorial, and example validation.
- `tutorial-example-authoring` — creating or maintaining teaching scripts.
- `isetcam-pipeline` — Scene, OI, sensor, IP, and display work.
- `sceneeye-maintenance` — sceneEye tutorial/example cleanup and promotion.
- `matlab-script-review` — read-only review of scripts against nearby tests.

## When Uncertain

Choose the smallest implementation consistent with nearby code. Ask before a
choice that materially changes behavior, public API shape, or test expectations.
