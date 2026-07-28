---
name: isetcam-pipeline
description: Implement, modify, debug, or review ISETBio image-system pipeline code. Use when tasks involve scene, optical image, sensor, image processing, or display objects; scene patterns or charts; color transforms; object constructors, getters, setters, compute methods, or diagnostic plots.
---

# ISETCam Pipeline Development

## Instructions & Guidelines

### Use the Existing Object Families

Follow the pipeline and its established APIs before adding utilities:

1. Scene: `scene*`, `sceneGet`, `sceneSet`.
2. Optical image: `oi*`, `oiGet`, `oiSet`.
3. Sensor: `sensor*`, `sensorGet`, `sensorSet`.
4. Image processing: `ip*`, `ipGet`, `ipSet`.
5. Display: `display*`, `displayGet`, `displaySet`.

Typical lifecycle calls are `sceneCreate`, `oiCreate`, `oiCompute`,
`sensorCreate`, `sensorCompute`, `ipCreate`, `ipCompute`, and
`displayCreate`.

For diagnostics, prefer the corresponding `scenePlot`, `oiPlot`, `sensorPlot`,
`ipPlot`, or `displayPlot` function over ad hoc figures.

### Search Before Implementing

- Search for nearby implementations using the relevant object prefix.
- Search `color/` before adding color transforms or color-science utilities.
- For scene patterns or charts, inspect `scene/` and related existing pattern
  or chart code before designing new behavior.
- Respect existing property names, input/output conventions, data shapes, and
  public APIs. Avoid parallel helper implementations when an object-family API
  already expresses the operation.

Validate behavior with the nearest colocated tests and a compact object-level
workflow when it clarifies correctness.
