function sceneWBCreate(sceneAll, workDir)
% Create a directory of waveband scene images in separate files
%
% Syntax:
%	sceneWBCreate(sceneAll, workDir)
%
% Description:
%    By default, this routine produces a directory with the scene name in
%    the current working directory that contains a set of  Matlab (.mat)
%    scene files. Each scene file contains photon data for a particular
%    wavelength. The files are named sceneXXX.mat where XXX is the
%    wavelength center of that image in nanometers.
%
%    Using this routine, along with the other waveband computational
%    routines, the user can run a larger spatial image by processing the
%    data from scene to sensor one waveband at a time.
%
% Inputs:
%    sceneAll - The original (complete) scene
%    workDir  - The working directory
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Fix example. Error is:
%           Error using sceneInterpolateW (line 69)
%           Wavelength extrapolation not allowed. Only interpolation
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneWBCreate.m' in the command window.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/20/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip.  This example throws an error. Fix and then delete these
    % comment lines to add back to autorun.

    scene = sceneCreate;
    sceneWBCreate(scene);
    sceneWBCreate(scene, pwd);
%}

if notDefined('sceneAll'), error('You must define the scene'); end

name = sceneGet(sceneAll, 'name');
if notDefined('workDir'), workDir = fullfile(pwd, name);  end
if ~exist(workDir, 'dir')
    w = warndlg('Creating work directory.');
    [p, n] = fileparts(workDir);
    chdir(p); mkdir(n); close(w);
end

curDir = pwd;
chdir(workDir);

nWave = sceneGet(sceneAll, 'nwave');
wave = sceneGet(sceneAll, 'wave');

scene = sceneClearData(sceneAll);
for ii = 1:nWave
    photons = sceneGet(sceneAll, 'photons', wave(ii));
    scene = sceneSet(scene , 'wave', wave(ii));
    scene = sceneSet(scene, 'photons', photons);
    fname = sprintf('scene%.0d.mat', wave(ii));
    vcSaveObject(scene, fname);
end

chdir(curDir);

end