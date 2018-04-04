function scene = sceneGridLines(scene, sz, lineSpacing, spectralType)
% Create scene comprising an array of grid lines
%
% Syntax:
%	scene = sceneGridLines(scene, [sz], [lineSpacing], [spectralType])
%
% Description:
%    Create a scene comprising an array of grid lines.
%
%    The grid line scene is useful for visualizing the geometric distortion
%    of a lens. The spectral power distribution of the lines is set to
%    equal photons ('ep'), unless spectralType is set to one of 'd65' or
%    'ee' (equal energy).
%
% Inputs:
%    scene        - A scene structure
%    sz           - (Optional) The scene size. Default 128.
%    lineSpacing  - (Optional) The number of samples between the lines.
%                   Default 16.
%    spectralType - (Optional) The spectral power distribution. Default
%                   'ep'. Other options include 'ee' and 'd65'.
%
% Outputs:
%    scene        - The modified scene structure.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    02/02/18  jnm  Formatting. Fix descr. to mention 'correct' default.

% Examples:
%{
    scene = sceneCreate;
    scene = sceneGridLines(scene);
    scene = sceneGridLines(scene, 128, 16, 'ee');
    scene = sceneGridLines(scene, 128, 16, 'ep');
%}


if notDefined('scene'), error('Scene structure required'); end
if notDefined('sz'), sz = 128; end
if notDefined('lineSpacing'), lineSpacing = 16; end
if notDefined('spectralType'), spectralType = 'ep'; end

scene = sceneSet(scene, 'name', 'gridlines');

scene = initDefaultSpectrum(scene, 'hyperspectral');
wave = sceneGet(scene, 'wave');
nWave = sceneGet(scene, 'nwave');

d = zeros(sz);
d(round(lineSpacing / 2):lineSpacing:sz, :) = 1;
d(:, round(lineSpacing / 2):lineSpacing:sz) = 1;

% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero. This is due to the
% ieCompressData scheme.
d(d == 0) = 1e-4;

switch lower(spectralType)
    case {'d65'}
        spd = ieReadSpectra('D65', wave);
        illPhotons = Energy2Quanta(wave, spd);
        % spect = Energy2Quanta(wave, illPhotons);
    case {'ee', 'equalenergy'}
        illPhotons = Energy2Quanta(wave, ones(nWave, 1));
    case {'ep', 'equalphoton'}
        illPhotons = ones(nWave, 1);
    otherwise
        error('Unknown spectral type:%s\n', spectralType);
end

data = bsxfun(@times, d, reshape(illPhotons, [1 1 nWave]));

scene = sceneSet(scene, 'illuminantPhotons', illPhotons);

% Allocate space for the photons
scene = sceneSet(scene, 'photons', data);
scene = sceneSet(scene, 'fov', 40);

end