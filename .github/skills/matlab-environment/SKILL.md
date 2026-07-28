---
name: matlab-environment
description: Configure or troubleshoot a local MATLAB and VS Code environment for ISETBio. Use when setting MATLAB paths, configuring the MathWorks VS Code extension, starting MATLAB Desktop or Live Editor, or diagnosing MATLAB session startup and connection behavior.
---

# MATLAB Environment

## Instructions & Guidelines

Use a MATLAB session with both ISETCam and ISETBio, including subdirectories,
on its path. ISETCam is required for ISETBio development and testing.

Run the following in the MATLAB Command Window for the session you intend to
use:

```matlab
addpath(genpath(fullfile(getenv('HOME'),'Documents','MATLAB','isetcam')));
addpath(genpath(fullfile(getenv('HOME'),'Documents','MATLAB','isetbio')));

which ieInit
```

For repeated `parfor`-heavy tutorial or example work, optionally warm the
pool after path setup:

```matlab
ieParallelPoolWarmUp('config','conservative','runSilent',true);
```

### VS Code

Install the recommended `MathWorks.language-matlab` extension. The repository
uses `MATLAB.matlabConnectionTiming: "onDemand"`; retain that workspace default
unless the team intentionally changes it. Set `MATLAB.installPath` in personal
VS Code settings when automatic discovery does not find MATLAB, for example:

```jsonc
"MATLAB.installPath": "/Applications/MATLAB_R2025b.app",
"files.associations": { "*.m": "matlab" }
```

When using MATLAB Desktop or Live Editor, start the Desktop first, wait for it
to become responsive, then open VS Code and a `.m` file. The extension may run
a separate background MATLAB process for language features; this is expected.

If a VS Code-started session prevents Desktop MATLAB from opening later, quit
both applications and stop only the lingering MathWorks helper processes before
relaunching the Desktop.

### Personal startup.m

Keep desktop-specific personal initialization from running in headless or VS
Code sessions. A guard can test `usejava('desktop')`, `VSCODE_PID`, and
`VSCODE_IPC_HOOK_CLI`; in a VS Code-only session, print the path-setup reminder
and, if needed for debugging, call `matlab.engine.shareEngine`.

Do not commit machine-specific MATLAB installation paths or personal
`startup.m` changes to this repository.
