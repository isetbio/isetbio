---
name: matlab-script-review
description: Review existing ISETBio MATLAB tutorial, example, or demo scripts without implementing behavior. Use when asked to assess runnability, comment quality, overlap, consolidation candidates, or teaching-script coverage against the nearest object-specific _tests_ directory, especially scripts under scripts/ or scene-related examples.
---

# MATLAB Script Review

## Instructions & Guidelines

### Scope and Operating Rules

Default to read-only analysis. Do not edit scripts, tests, or documentation
unless explicitly asked. Start from the named script directory and the nearest
`_tests_` directory for that object area. Keep the review local: inspect
unrelated areas only when they control the reviewed behavior.

When execution is requested, use the cheapest focused MATLAB command or the
closest test before proposing a broader run.

### Review Workflow

1. Inventory target scripts and group them by topic.
2. Inventory the nearest tests and map them to those topics.
3. For every script, assess purpose, APIs exercised, runtime dependencies or
   failure risks, comment quality, and overlap with neighboring scripts.
4. Compare demonstrations against tests: identify behavior tested but not
   taught, and repeated demonstrations with no distinct teaching value.
5. Recommend a minimal cleanup plan with concrete keep, merge, split, skip,
   or retire candidates.

Lead with findings, risks, and overlap candidates. Clearly distinguish what
was executed from static inference. Cite the closest tests; state when no
corresponding demonstration exists. For a merge recommendation, name the
overlapping scripts and the one teaching goal that should remain.

### Scene Review Heuristics

Prefer existing `scene*` constructors, accessors, and plotting helpers. Check
`s_sceneDemo.m` and `s_sceneExamples.m` as possible overview overlaps. Review
illuminant scripts (`s_sceneIlluminant.m`, `s_sceneIlluminantSpace.m`,
`s_sceneIlluminantMixtures.m`, `s_sceneChangeIlluminant.m`) as a cluster, and
do the same for reflectance scripts (`s_sceneReflectanceCharts.m`,
`s_sceneReflectanceChartBasisFunctions.m`, `s_sceneReflectanceSamples.m`).

Compare relevant scene examples with nearby tests including `test_scenedemo.m`,
`test_sceneexamples.m`, `test_sceneChangeIlluminant.m`,
`test_sceneIncreaseSize.m`, and `test_sceneHCCompress.m` when present.

Favor a small number of high-value scripts with distinct teaching goals. Keep a
script that uniquely explains a core concept even when a unit test covers the
same API surface; flag historical duplicates that add no teaching value.
