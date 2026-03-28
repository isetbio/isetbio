# Connecting MATLAB to VS Code for ISETBio (macOS)

## 1) Configure VS Code settings

To ensure the VS Code MATLAB extension finds your specific MATLAB installation, point it at the MATLAB app bundle. 

1. Open **Settings** in VS Code (`Cmd + ,`).
2. Search for `matlab.installPath`.
3. Set it to your MATLAB installation, for example:
   - `/Applications/MATLAB_R2024b.app` (or your specific version)

Alternatively, this workspace already provides a `.vscode/settings.json` overriding this value.

## 2) Update `startup.m` for environment detection

**ISETBio critically requires ISETCam to be on your MATLAB path.** The VS Code MATLAB extension typically runs MATLAB in a terminal/headless mode, bypassing your desktop UI path configurations.

To initialize ISETBio and ISETCam paths properly, have `startup.m` detect when it is running under VS Code:

```matlab
% Detect if MATLAB is being launched by VS Code
isVSCode = ~usejava('desktop') || ...
    ~isempty(getenv('VSCODE_PID')) || ...
    ~isempty(getenv('VSCODE_IPC_HOOK_CLI'));

if isVSCode
    disp('Initializing ISETBio paths for VS Code');
    
    % --- ISETCAM REQUIRED PATH ---
    % Since ISETBio depends entirely on ISETCam functions (scene, oi, etc.), 
    % you must add ISETCam to your path before running ISETBio code.
    addpath(genpath('/Users/wandell/Documents/MATLAB/isetcam'));
    
    % --- ISETBIO PATH ---
    addpath(genpath('/Users/wandell/Documents/MATLAB/isetbio'));
    
    % Share the engine so the VS Code extension can connect for debugging
    matlab.engine.shareEngine;

    fprintf('MATLAB initialized for VS Code workspace.\n');
    return
elseif isdeployed
    % Skip initialization for compiled apps
else
    % Standard Desktop initialization
    reset(groot);
    % Your usual plotting/graphics defaults here
end
```

## 3) Verification

- **Start MATLAB:** Click the MATLAB icon in the Activity Bar or open a `.m` file. The extension should start a MATLAB session in the integrated terminal.
- **Path check:** In the VS Code MATLAB terminal, run:

  ```matlab
  path
  ```

  Verify that **both** the `isetcam` and `isetbio` directories are included. If `isetcam` is missing, you will face "Undefined function or variable" errors for core objects like `sceneCreate` or `oiCreate`.
