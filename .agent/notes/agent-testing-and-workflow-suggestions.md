# Agent-Based Coding Suggestions for ISETBio

Date: 2026-03-28

## Context
ISETBio is transitioning from mostly system-wide validation scripts (e.g., `v_ibio_*`) to supporting local agent workflows, mirroring the strategies adopted in the sister repository `ISETCam`. 

This note outlines how AI coding agents and human contributors should approach writing and testing new functionality safely.

## Recommended Testing Strategy (Hybrid)

### 1) System Level: `isetvalidate`
- Validations such as `v_ibio_mRGCMosaic` or similar should be maintained for long-term behavioral coverage.
- Treat validation scripts as the primary guarantee of mathematical properties passing end-to-end biological modeling.

### 2) Feature Level: Function-Local Tests in `isetbio`
- Sibling `_tests_` directories inside component folders. For instance, if you add complex geometry to a mosaic algorithm, colocate standard assertions near that function.
- Naming format: `ieTest_<functionname_or_class>.m`.

### 3) Test Structure per Function
Each `ieTest_` script in ISETBio should generally address:
1. Validating input arguments and properties.
2. Checking spatial/temporal output dimensions (very important for overlapping mosaics).
3. Using golden values for deterministic algorithms given a known seed. 

## Migration Guidance for `runtests` / `matlab.unittest`
- Do not migrate global scripts without explicit command. 
- Use standard MATLAB `.m` scripts containing raw assertions unless user requests `matlab.unittest`.

## Agent Workflow Suggestions for `isetbio/.agent`

Agents operating in ISETBio should establish and enforce policies through this folder structure:

- `.agent/notes/` — Core policies, design choices, and migration tasks.
- `.agent/workflows/` — Checklists for routine chores (upgrading cone models, building a new dataset).
- `.agent/templates/` — Snippets for consistent file initialization.

## Typical Workflows for AI Agents in ISETBio

1. **Bugfix/Refactor Workflow**
   - Ensure the problem is locally reproducible.
   - Run the corresponding `isetvalidate` script.
   - Isolate the feature bug and draft a failing case in `_tests_`.
   - Iterate to green, making sure the return formats (`cMosaic` struct signatures) do not change unexpectedly.

2. **Feature Workflow**
   - Use existing abstract classes (e.g. extending an abstract base RGC model) if applicable.
   - Maintain naming conventions aligned with `setup/copilot-instructions.md`.
   - Update `Contents.m` headers in relevant subdirectories to aid users' `help` queries. 

## Key Technical Details
- **External Dependencies:** Because ISETBio strictly depends upon `ISETCam`, assume agents might run code that expects paths already set. Agents should not blindly install missing functions if those functions conventionally reside in `ISETCam`. Instead, ensure the user has updated both repositories.
