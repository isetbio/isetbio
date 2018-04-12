function scene = sceneSpatialResample(scene, dx, units, method, promptUser)
% Spatial resample all wavebands of a scene
%
% Syntax:
%	scene = spatialResample(scene, dx, [units], [method], [promptUser])
%
% Description:
%    Spatial resample all wavebands of a scene
%
%    N.B. The source contains executable examples of usage, which can be
%    accessed by typing 'edit sceneSpatialResample.m' in the command window
%
% Inputs:
%    scene      - ISET scene
%    dx         - New sample spacing
%    units      - (Optional) Sample spatial units. Default 'm' (for meters)
%    method     - (Optional) linear, cubic or spline interpolation.
%                 Default 'linear'.
%    promptUser - (Optional) Whether or not to prompt the user to begin.
%                 Default true. Set to false to avoid waiting for user to
%                 press a key.
%
% Outputs:
%    scene      - The ISET scene after modification
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    sceneSpatialSupport, oiSpatialResample
%

% History:
%    xx/xx/16  BW   Copyright ISETBIO Team, 2016
%    12/20/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    scene = sceneSet(scene, 'fov', 3);
    ieAddObject(scene);
    sceneWindow;
%}
%{
    % ETTBSkip. This example provides user input. Skip on autorun.
    scene = sceneCreate;
    scene = sceneSpatialResample(scene, 100, 'um');
    ieAddObject(scene);
    sceneWindow;
%}

%% Set up parameters
if notDefined('scene'), error('scene required'); end
if notDefined('dx'), error('Sample spacing required'); end
if notDefined('units'), units = 'm'; end
if notDefined('method'), method = 'linear'; end
if notDefined('promptUser'), promptUser = true; end
% Always work in meters
dx = dx / ieUnitScaleFactor(units);

mLum = sceneGet(scene, 'mean luminance');

% Find the spatial support of the current scene, and its max/min
ss = sceneSpatialSupport(scene, 'meters'); % x and y spatial support
xmin = min(ss.x(:));
xmax = max(ss.x(:));
ymin = min(ss.y(:));
ymax = max(ss.y(:));

% Set up the new spatial support
% We want height/rows = dx exactly, if possible
% We get to set the FOV to make this work out.
xN = xmin:dx:xmax;
yN = ymin:dx:ymax;

% fprintf('Current  dx = %f meters\n', ss.y(2) - ss.y(1));
% fprintf('Proposed dx = %f meters\n', dx);
% fprintf('New scene size %d (rows) %d (cols)\n', length(yN), length(xN));
if length(xN) > 1000 || length(yN) > 1000
    if (promptUser)
        fprintf('Very large scene. Any key to continue\n');
        pause
    end
end

%% Interpolate the image for each waveband
nWave = sceneGet(scene, 'nwave');
wave = sceneGet(scene, 'wave');

% Precompute meshgrid for speed outside of loop
[X, Y] = meshgrid(ss.x, ss.y);
[Xq, Yq] = meshgrid(xN, yN);

photonsN = zeros(length(yN), length(xN), nWave);
for ii = 1:nWave
    photons = sceneGet(scene, 'photons', wave(ii));
    photonsN(:, :, ii) = interp2(X, Y, photons, Xq, Yq, method);
end

% Change up the photons and thus the row/col
scene = sceneSet(scene, 'photons', photonsN);
n = sceneGet(scene, 'name');
scene = sceneSet(scene, 'name', sprintf('%s-%s', n, method));

% Now adjust the FOV so that the dx works out perfectly
sr = sceneGet(scene, 'spatial resolution');
fov = sceneGet(scene, 'fov');
scene = sceneSet(scene, 'fov', fov * dx / sr(2));

% The spatial resampling can have a small effect
scene = sceneSet(scene, 'mean luminance', mLum);

end
