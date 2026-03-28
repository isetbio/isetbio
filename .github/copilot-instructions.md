# Copilot Instructions for 'isetbio' Workspace

## Purpose
Use these instructions to navigate ISETBio-style code quickly and choose project-consistent functions before writing new code. 

## Matlab setup

For instructions on how to setup matlab for use with VS Code, see [this guide](../.vscode/matlab-setup.md).

## Workspace Context
This workspace heavily depends on another repository, **ISETCam**. 
You must treat `isetcam` and `isetbio` as a cohesive ecosystem.

- `isetcam`: Provides the foundation for generic scene radiance and simple optical formulation.
- `isetbio`: Provides precise biologically accurate models of human vision.
- `isetvalidate` / `ISETValidations`: The standalone tool for unit tests of both.

Assume MATLAB is the primary runtime.

## Core Naming Conventions
When searching for functionality, ISETBio mixes `classdef` style structures (newer code) with traditional structure-passing functions (older code). 

### ISETCam Base Objects
ISETBio routines will frequently create and consume the following base objects from ISETCam:
- **Scene**: `scene*` methods.
- **Optical Image**: `oi*` methods. 
  - *Exception*: ISETBio provides specialized human optics, frequently using `sceneEye*` or custom `oiCompute*` methods.

### ISETBio Core Objects (Biological Hierarchy)
1. **sceneEye**: Human wavefront optics for turning scene radiance into retinal irradiance.
2. **Macular**: `Macular` class methods for macular pigment transmittance.
3. **Cone Photoreceptors**: 
   - `cMosaic` (modern wrapper for cone mosaic construction and responses)
   - `coneMosaic` (older cone abstraction)
   - `outersegment*` or `os*`
4. **Bipolar Cells**: `bipolar*` or `bipolarLayer`
5. **Retinal Ganglion Cells**: 
   - `mRGCMosaic` (modern midget RGC mosaic)
   - `irgcMosaic` (inner retina ganglion cells) 
   - `rgc*` functions

## Plotting Conventions
Always prefer object-specific `plot()` or `visualize()` methods provided by the underlying classes rather than ad-hoc plotting routines. 

- For traditional structures from ISETCam, use: `scenePlot`, `oiPlot`, `sensorPlot`.
- For specific ISETBio classes (like `cMosaic` or `mRGCMosaic`), invoke their intrinsic method, e.g., `mosaic.visualize()`.

## MATLAB Workflow Expectations
- Keep edits minimal and consistent with existing MATLAB style.
- Prefer vectorized operations for performance-sensitive loops, especially when modeling millions of cells in retinal mosaics.
- Validate modified files using local testing folders (`_tests_`) and system-level validation when extending objects.
- Be extremely cautious when modifying the geometric algorithms parsing retinal lattices (`retinalattice`) or calculating response outputs, as downstream models expect precise multidimensional arrays.
- **Terminal Tooling**: If executing sweeping codebase searches via a terminal, assume a macOS environment with `rg` (ripgrep) and `find_by_name` installed. Prefer using `rg` over `grep` or slower MATLAB search functions.

## If Uncertain
If multiple plausible implementations exist, prefer the newer object-oriented approach (e.g., `cMosaic` / `mRGCMosaic`) over legacy procedural scripts. If still unsure, ask the user for refinement.
