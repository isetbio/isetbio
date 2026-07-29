---
name: tutorial-example-authoring
description: Create, move, review, publish, or test ISETBio tutorials, examples, and data-generation scripts. Use when working in tutorials/ or examples/, naming t_*.m, s_*.m, or data_*.m files, choosing tutorial versus example scope, adding SkipFile, or publishing teaching scripts to HTML.
---

# Tutorial and Example Authoring

## Instructions & Guidelines

### Select the Correct Teaching Surface

Use `tutorials/` for short, heavily commented introductions for learners who
know programming but are learning image-systems engineering and ISETCam object
fundamentals. A tutorial should present object setup, key `*Get`/`*Set` calls,
basic `*Window` or `*Plot` output, and one simple quantitative checkpoint.
Keep it fast and readable linearly.

Use `examples/` for realistic, adaptable analysis patterns. Examples may be
longer and can include end-to-end numerical analysis, visualization workflows,
realistic parameters, and tradeoff exploration.

Use `data_*.m` for scripts that generate or refresh repository data. Do not
name those side-effecting scripts `t_*.m` or `s_*.m`.

### Authoring Standards

- Prioritize clarity, reproducibility, stable outputs, and instructional value.
- Use a clear header, readable section structure, and comments that explain the
  purpose of each section.
- Link to related tests, nearby teaching scripts, or wiki material when useful.
- Reuse existing constructors, accessors, and plotting helpers instead of
  embedding alternate APIs or ad hoc visualizations.
- Validate a changed teaching script with the corresponding selected runner:

```matlab
isetbioTutorialTest('selection','t_exampleName');
isetbioExampleTest('selection','s_exampleName');
```

### Skipping and Publishing

Only add `% SkipFile` when a script requires unavailable external data or
toolboxes, deliberate interaction, unusual expense, or a documented known
failure. State the reason nearby. New and repaired files should use the
canonical marker, not `% UTTBSkip`.

Publish tutorials or examples with the ISETCam-provided utilities:

```matlab
s_publishTutorials
s_publishExamples
iePublish('filename.m')
```

`iePublish` applies the expected HTML formatting and embeds figure styles.
