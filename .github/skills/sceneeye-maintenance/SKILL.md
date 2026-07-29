---
name: sceneeye-maintenance
description: Repair, promote, reorganize, or evaluate sceneEye tutorials and examples. Use when editing tutorials/sceneEye, examples/sceneEye, sceneEye underDevelopment paths, PBRT or Docker/cloud-dependent scripts, or deciding whether a sceneEye workflow belongs in a routine tutorial smoke test or a skipped example.
---

# sceneEye Maintenance

## Instructions & Guidelines

### Placement Rule

Keep a sceneEye file in `tutorials/sceneEye` only when it is short, linear, and
teaches one core object/API concept using a compact render or mocked/precomputed
result. Avoid cloud execution, manual Docker setup, large PBRT sweeps, and
high-resolution comparison grids.

Place applied workflows in `examples/sceneEye`: multiple eye models or render
conditions, crop-window rendering, depth-of-field or accommodation sweeps,
retina-shape experiments, stereo, cloud rendering, PSF/MTF analysis, and
publication-style figures.

Keep `tutorials/sceneEye/analysis/underDevelopment/` and
`tutorials/sceneEye/cloud/underDevelopment/` skipped unless a deliberate
promotion makes one script routine, local, and reliable.

### Repair and Promotion Workflow

1. Identify obsolete helpers, unavailable assets, external dependencies, and
   repeated high-resolution work.
2. Replace legacy dependency checks such as `piCamBio`, `piDockerExists`,
   `piDockerConfig`, `mcDockerExists`, and `mcDockerConfig` with a current
   documented setup check or an explicit skip reason.
3. Reduce repeated parameter sweeps to one representative runnable condition.
4. If a PBRT scene is not in the repository or standard data-download path,
   document the asset requirement clearly or retain the script as a skipped
   example.
5. Test a promoted file with its selected tutorial or example runner.

### Acceptance Criteria

- `isetbioTutorialTest` has zero failures.
- Under-development paths remain skipped.
- A restored tutorial passes `isetbioTutorialTest('selection','<name>')` and
  generally completes in under 10 seconds on a normal local MATLAB session.
- A moved example passes `isetbioExampleTest('selection','<name>')`, or keeps
  `% SkipFile` with a specific documented reason.
- Retain one clear example per workflow family; merge, remove, or skip
  obsolete, overly narrow, or duplicated scripts.
