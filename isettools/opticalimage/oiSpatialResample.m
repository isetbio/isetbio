function oi = oiSpatialResample(oi, dx, units, method, promptUser)
% Spatial resample all wavebands of a scene
%
% Syntax:
%   oi = oiSpatialResample(oi, dx, units, method, promptUser)
%
% Description:
%    The current oi spatial sampling can be measured
%
%    There are examples contained in the code. To access, typed 'edit
%    oiSpatialResample.m' into the Command Window.
%
% Inputs:
%    oi         - Struct. An optical image structure.
%    dx         - Scalar Numeric. New spatial sampling difference
%    units      - (Optional) String. Spatial units. Default 'm' (meters).
%    method     - (Optional) String. Interpolation method. Default linear.
%    promptUser - (Optional) Boolean. Whether to wait for a key press from
%                 the user. Set to false to avoid the routine waiting for
%                 the user to enter a keypress. Default true (wait).
%
% Outputs:
%    oi         - Struct. The modified optical image structure.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    v_sceneSpatialResample (tests oi, too)
%

% History:
%    xx/xx/16       Copyright Imageval Consulting, LLC 2016
%    03/07/18  jnm  Formatting
%    06/24/19  JNM  Minor formatting adjustments

% Examples:
%{
    scene = sceneCreate;
    scene = sceneSet(scene, 'fov', 10);
    oi = oiCreate('diffraction limited');
    oi = oiCompute(oi, scene);
    ieAddObject(oi);
    oiWindow;  % 7.0 um samples

    oi = oiSpatialResample(oi, 1e-6);  % 1.0 um samples
    ieAddObject(oi);
    oiWindow;
%}

%% Set up parameters
if notDefined('oi'), error('oi required'); end
if notDefined('dx'), error('dx required'); end
if notDefined('units'), units = 'm'; end
if notDefined('method'), method = 'linear'; end
if notDefined('promptUser'), promptUser = true; end
%
% Always work in meters
dx = dx / ieUnitScaleFactor(units);

% Find the spatial support of the current oi, and its max/min
ss = oiSpatialSupport(oi, 'meters');  % x and y spatial support
xmin = min(ss.x(:));
xmax = max(ss.x(:));
ymin = min(ss.y(:));
ymax = max(ss.y(:));

% Set up the new spatial support
xN = xmin:dx:xmax;
yN = ymin:dx:ymax;

% fprintf('Current  dx = %f meters\n', ss.y(2) - ss.y(1));
% fprintf('Proposed dx = %f meters\n', dx);
% fprintf('New scene size %d (rows) %d (cols)\n', length(yN), length(xN));
if length(xN) > 1000 || length(yN) > 1000
    if (promptUser)
        fprintf('Very large scene. Press any key to continue\n');
        pause
    end
end

%% Interpolate the image for each waveband
nWave = oiGet(oi, 'nwave');
wave = oiGet(oi, 'wave');

% Precompute meshgrid for speed outside of loop
[X, Y] = meshgrid(ss.x, ss.y);
[Xq, Yq] = meshgrid(xN, yN);

photonsN = zeros(length(yN), length(xN), nWave);
for ii = 1:nWave
    photons = oiGet(oi, 'photons', wave(ii));
    photonsN(:, :, ii) = interp2(X, Y, photons, Xq, Yq, method);
end

% Change up the photons and thus the row/col
oi = oiSet(oi, 'photons', photonsN);
n = oiGet(oi, 'name');
oi = oiSet(oi, 'name', sprintf('%s-%s', n, method));

% Now adjust the FOV so that the dx works out perfectly
sr = oiGet(oi, 'spatial resolution');
fov = oiGet(oi, 'fov');
oi = oiSet(oi, 'fov', fov * dx / sr(2));

oi = oiSet(oi, 'illuminance', oiCalculateIlluminance(oi));

end
