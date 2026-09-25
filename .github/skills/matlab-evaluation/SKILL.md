---
name: matlab-evaluation
description: Run MATLAB non-interactively from a shell on this machine — locating the MATLAB binary, setting ISETCam/ISETBio paths for a batch session, and invoking -batch evaluation. Use when a shell command must execute MATLAB code, a test runner, or iePublish; use matlab-testing to choose which tests to run.
---

# Matlab Evaluation Workflow

## Finding Matlab

MATLAB is not on `$PATH` on this machine, so call the executable by its full
path. Installed releases live in `/Applications` as `MATLAB_R20XXx.app`
bundles, with the executable at `bin/matlab`. List them rather than assuming a
version, and use the newest unless the user specifies otherwise or it fails
for a version-specific reason:

```bash
ls /Applications | rg -i '^MATLAB_R'
```

## Running non-interactively

Use `-batch`, not `-r`. `-batch` runs headless, still processes MATLAB's
startup files, exits when the statement finishes, and returns a non-zero exit
code with the error on stderr if the statement throws. `-r` needs an explicit
`exit` and a `try/catch` to avoid hanging or masking failures. Check `$?`
after the call; terminal text alone does not reliably signal failure.

```bash
cd ~/Documents/MATLAB/isetbio
/Applications/MATLAB_R2026a.app/bin/matlab -batch "\
addpath(genpath('~/Documents/MATLAB/isetcam')); \
addpath(genpath('~/Documents/MATLAB/isetbio')); \
results = isetbioUnitTest;"
echo "exit: $?"
```

Personal `startup.m` path setup is usually guarded by `usejava('desktop')` and
skipped in batch sessions, so add ISETCam and ISETBio explicitly as above. Add
only the toolboxes the code needs; never `addpath(genpath('~/Documents/MATLAB'))`,
which mixes unrelated projects. Use `which functionName` to confirm resolution.
`Unrecognized function or variable` usually means a missing path.

Scripts that call `pause`, `input`, `keyboard`, or `waitfor` block a batch
session with nobody to answer. Run unfamiliar code with a time limit — the
shell tool's own timeout; macOS has no `timeout` command unless coreutils is
installed.

## Working directory

`-batch` inherits the shell's working directory as MATLAB's `pwd`. `iePublish`
and some scripts resolve relative paths and write output relative to `pwd`, so
`cd` to the intended directory first.
