function sceneW = sceneExtractWaveband(scene,waveList)
%sceneExtractWaveband - Extract wave bands from the scene
%
%   sceneW = sceneExtractWaveband(scene,waveList)
%
% The list of evenly-spaced wavelengths, waveList in nm, is extracted from
% the original scene. The output scene contains the photons in the
% wavelength bands. The new scene does not have a calculated luminance, and
% in fact its luminance differs from the original scene.
%
% If the waveList is a single value, the spectral bin width is set to 1.
% Otherwise it is set to the difference in the (evenly spaced!) wavelength
% list.
%
%Example
%   sceneMonochrome = sceneExtractWaveband(scene,500);
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('scene'), scene = vcGetObject('scene'); end
if notDefined('waveList'), error('Wave list must be defined'); end

sceneW = scene;
sceneW = sceneSet(sceneW,'cphotons',sceneGet(scene,'photons',waveList));
sceneW = sceneSet(sceneW,'wave',waveList);

end