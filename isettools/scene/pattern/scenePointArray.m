function scene = scenePointArray(scene, sz, pointSpacing, spectralType)
% Make a point array stimulus for evaluating the optics
%
% Syntax:
%	scene = scenePointArray(scene, [sz], [pointSpacing], [spectralType])
%
% Description:
%    The point array scene clarifies the PSF at a variety of locations in
%    the optical image. It also provides a sense of the geometrical
%    distortions in the scene.
%
% Inputs:
%    scene        - A scene structure
%    sz           - (Optional) Row & column size of the image. Default 128.
%    pointSpacing - (Optional) The distance between points in pixels.
%                   Default is 16.
%    spectralType - (Optional) The spectrum. Default 'D65' Other options
%                   include spectrums such as 'ee' (equal energy) and 'ep'
%                   (equal photons).
%
% Outputs:
%    scene        - The modified scene structure.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    02/02/18  jnm  Formatting

% Examples:
%{
    scene = sceneCreate;
    scene = scenePointArray(scene);
    scene = scenePointArray(scene, 64, 8, 'd65');
    scene = scenePointArray(scene, 64, 8, 'ee');
    scene = scenePointArray(scene, 64, 8, 'ep');
%}

if notDefined('scene'), error('Scene structure required'); end
if notDefined('sz'), sz = 128; end
if notDefined('pointSpacing'), pointSpacing = 16; end
if notDefined('spectralType'), spectralType = 'd65'; end

scene = sceneSet(scene, 'name', 'pointarray');

scene = initDefaultSpectrum(scene, 'multispectral');
wave  = sceneGet(scene, 'wave');
nWave = sceneGet(scene, 'nwave');

d = zeros(sz);
idx = round(pointSpacing / 2):pointSpacing:sz;
d(idx, idx) = 1;

switch lower(spectralType)
    case {'d65'}
        illPhotons = Energy2Quanta(wave, ieReadSpectra('D65', wave));
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