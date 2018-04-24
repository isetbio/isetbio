function sceneW = sceneExtractWaveband(scene, waveList)
% Extract wave bands from the scene
%
% Syntax:
%	sceneW = sceneExtractWaveband(scene, waveList)
%
% Description:
%    The list of evenly-spaced wavelengths, waveList in nm, is extracted
%    from the original scene. The output scene contains the photons in the
%    wavelength bands. The new scene does not have a calculated luminance,
%    and in fact its luminance differs from the original scene.
%
%    If the waveList is a single value, the spectral bin width is set to 1.
%    Otherwise it is set to the difference in the (evenly spaced!)
%    wavelength list.
%
%    There are examples in the code. Type 'edit sceneExtractWaveband' into
%    the Command Window to access.
%
% Inputs:
%    scene - (Optional) The scene structure. Default is to select an
%            existing scene.
%    waveList - The list of wavelengths
%
% Outputs:
%    sceneW   - The modified scene
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - The example does not work as desired. Either a single
%      value, or with a range (ex. 400:10:500).]
%    * TODO: Fix example.
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneExtractWaveband.m' in MATLAB's
%      command window.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/22/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match the Wiki.

% Examples:
%{
    scene = sceneCreate('macbethD65');
	sceneMonochrome = sceneExtractWaveband(scene, 500);
%}

if notDefined('scene'), scene = vcGetObject('scene'); end
if notDefined('waveList'), error('Wave list must be defined'); end

sceneW = scene;
sceneW = sceneSet(sceneW, 'wave', waveList);

end